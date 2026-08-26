"""
adapters/myheritage_v5.py
Адаптер MyHeritage (.csv/.tsv) -> ParseResult.

Формат файла MyHeritage на практике не стабилен:
 - Разделитель может быть как TSV (табуляция), так и CSV (запятая).
 - Строка заголовка может отличаться регистром и точным набором слов
   (RSID/RS_ID/SNP_ID, CHROMOSOME/CHROM/CHR, POSITION/POS/BP,
   RESULT/GENOTYPE/GT/ALLELE), в кавычках или без.
 - Перед заголовком обычно идёт ~12 строк комментариев ('#').

Задача 9: автоопределение разделителя и гибкий (case-insensitive,
с поддержкой синонимов) разбор заголовка вместо жёсткого требования
точного TSV-заголовка ('RSID','CHROMOSOME','POSITION','RESULT').

Логика разрешения ориентации аллелей идентична FTDNA — используется
тот же ReferenceGenome и те же QC-счётчики. На выходе — тот же
ParseResult, что и у ftdna_v3.

=============================================================================
ИЗМЕНЕНИЯ (промт "HRC / TopMed", интеграция лифтовера в адаптеры):
=============================================================================

Тот же баг и то же исправление, что и в adapters/ftdna_v3.py (см. его
докстринг для подробностей): parse_myheritage_v5() не принимала параметр
liftover, хотя main.py/gui/app.py вызывали её именно так, что приводило к
    TypeError: parse_myheritage_v5() got an unexpected keyword argument 'liftover'
Добавлен параметр liftover: Optional[ChainLiftover] = None, перенос
координаты происходит сразу после парсинга pos/chrom, до broad_key и до
reference.base_at().
"""
from __future__ import annotations
import csv
import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

try:
    from pyfaidx import Fasta  # noqa: F401  (нужен для ReferenceGenome, импортируется ниже)
except ImportError as exc:
    raise ImportError("Требуется пакет pyfaidx: pip install pyfaidx") from exc

from .base import ParsedVariant, ParseResult
from .ftdna_v3 import FTDNAFormatError, StrandQualityError, ReferenceGenome
from core.liftover import ChainLiftover

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Константы
# ---------------------------------------------------------------------------
NO_CALL = ("--", "NA", "N/A", "", "0")

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
SELF_COMPLEMENTARY_PAIRS = (frozenset("AT"), frozenset("CG"))

DEFAULT_BOTH_NON_REF_THRESHOLD_PCT = 0.1

# Нормализация хромосом MyHeritage:
#   "23" -> "X", "24" -> "Y", "25" -> "MT", "M" -> "MT"
CHROM_NORMALIZE = {
    "23": "X", "24": "Y", "25": "MT",
    "X": "X", "Y": "Y", "MT": "MT", "M": "MT",
    "chrX": "X", "chrY": "Y", "chrM": "MT", "chrMT": "MT",
}

REJECT_INVALID = "invalid"
REJECT_SELF_COMPLEMENTARY = "self_complementary"
REJECT_BOTH_NON_REF = "both_non_ref"

# Синонимы названий колонок (все сравниваются в нижнем регистре,
# без кавычек и пробелов по краям) — Задача 9, пункт 2.
COLUMN_SYNONYMS: dict[str, set[str]] = {
    "RSID": {"rsid", "rs_id", "rs id", "snp_id", "snp id", "marker", "markername", "id"},
    "CHROMOSOME": {"chromosome", "chrom", "chr"},
    "POSITION": {"position", "pos", "bp", "basepair", "base_pair"},
    "RESULT": {"result", "genotype", "gt", "allele", "alleles", "call"},
}
REQUIRED_COLUMNS = ("RSID", "CHROMOSOME", "POSITION", "RESULT")
MIN_MATCHED_COLUMNS = 3  # Задача 9, пункт 3: 3 из 4 полей достаточно


# ---------------------------------------------------------------------------
# Исключения
# ---------------------------------------------------------------------------
class MyHeritageFormatError(ValueError):
    pass


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------
def _normalize_chrom(chrom: str) -> str:
    """Нормализует хромосому MyHeritage: 23->X, 24->Y, 25->MT."""
    c = chrom.strip().strip('"').strip("'")
    if c.lower().startswith("chr") and c not in CHROM_NORMALIZE:
        c = c[3:]
    return CHROM_NORMALIZE.get(c, c)


def _strip_quotes(value: str) -> str:
    """Удаляет кавычки вокруг значения (MyHeritage часто их ставит)."""
    v = value.strip()
    if len(v) >= 2 and v[0] == v[-1] and v[0] in ('"', "'"):
        return v[1:-1]
    return v


def _clean_header_token(token: str) -> str:
    return _strip_quotes(token).strip().lower()


@dataclass
class _Resolved:
    ref: str
    alt: str
    gt: str
    reject_reason: str | None


def _resolve_genotype(result: str, ref_base: str) -> _Resolved:
    """
    Разрешает генотип MyHeritage относительно референса.
    Логика идентична FTDNA: если обе буквы совпадают с референсом -> 0/0,
    если одна совпадает -> 0/1, если обе другие -> 1/1,
    если обе не-референсные и разные -> both_non_ref (отклоняем).
    """
    if len(result) != 2 or any(b not in COMPLEMENT for b in result):
        return _Resolved(ref="", alt="", gt="", reject_reason=REJECT_INVALID)

    a1, a2 = result[0], result[1]
    alleles = frozenset(result)

    if len(alleles) == 2 and alleles in SELF_COMPLEMENTARY_PAIRS:
        return _Resolved(ref="", alt="", gt="",
                         reject_reason=REJECT_SELF_COMPLEMENTARY)

    if a1 == a2:
        if a1 == ref_base:
            return _Resolved(ref=ref_base, alt=".", gt="0/0", reject_reason=None)
        return _Resolved(ref=ref_base, alt=a1, gt="1/1", reject_reason=None)

    non_ref = [a for a in (a1, a2) if a != ref_base]
    if len(non_ref) == 1:
        return _Resolved(ref=ref_base, alt=non_ref[0], gt="0/1", reject_reason=None)
    return _Resolved(ref="", alt="", gt="", reject_reason=REJECT_BOTH_NON_REF)


def _chip_signature(positions: list[tuple[str, int]]) -> str:
    h = hashlib.sha256()
    for chrom, pos in positions:
        h.update(f"{chrom}:{pos}\n".encode("utf-8"))
    return h.hexdigest()[:16]


def save_position_cache(cache_dir: Path, result: ParseResult) -> Path:
    """
    Сохраняет позиции чипа (chrom, pos) в
    cache_dir/<chip_signature>.positions.json — тем же самым результатом
    парсинга (result.variants), который уже прошёл через гибкое
    определение разделителя/заголовка MyHeritage, а не повторным
    (потенциально неверным) разбором сырого CSV. download_donors.py
    читает этот файл через --positions-json вместо create_chip_positions().
    Идентична save_position_cache() в ftdna_v3.py/vcf_source.py — единый
    интерфейс для main.py::SOURCES[...]["save_position_cache"].
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature}.positions.json"
    if not out_path.exists():
        payload = [(v.chrom, v.pos) for v in result.variants]
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


def save_position_cache_broad(cache_dir: Path, result: ParseResult) -> Path:
    """
    Задача D (опционально) — см. докстринг adapters/base.py. Аналог
    save_position_cache(), но по chip_signature_broad/
    signature_positions_broad (все позиции физического дизайна чипа —
    см. исправление в parse_myheritage_v5() ниже). НЕ использовать вместе
    со старым флагом bcftools merge -0 — см. предупреждение в
    adapters/base.py.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature_broad}.positions.json"
    if not out_path.exists():
        payload = result.signature_positions_broad
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# Автоопределение разделителя (Задача 9, пункт 1)
# ---------------------------------------------------------------------------
def _detect_delimiter(sample_line: str) -> str:
    """
    Считает запятые и табуляции в строке-кандидате на заголовок и
    выбирает более частый символ как разделитель. При равенстве
    (или отсутствии обоих) по умолчанию — таб, так как исторически
    MyHeritage отдаёт TSV.
    """
    tabs = sample_line.count("\t")
    commas = sample_line.count(",")
    if commas > tabs:
        return ","
    return "\t"


# ---------------------------------------------------------------------------
# Гибкий разбор заголовка (Задача 9, пункт 2-4)
# ---------------------------------------------------------------------------
def _match_header(tokens: list[str]) -> Optional[dict[str, int]]:
    """
    Пытается сопоставить токены строки заголовка с REQUIRED_COLUMNS
    через синонимы. Возвращает {каноническое_имя: индекс_колонки} если
    удалось найти минимум MIN_MATCHED_COLUMNS полей, иначе None.
    """
    cleaned = [_clean_header_token(t) for t in tokens]
    found: dict[str, int] = {}
    for canonical, synonyms in COLUMN_SYNONYMS.items():
        for idx, token in enumerate(cleaned):
            if token in synonyms:
                found[canonical] = idx
                break
    if len(found) >= MIN_MATCHED_COLUMNS:
        return found
    return None


def _find_header(path: Path) -> tuple[int, str, dict[str, int]]:
    """
    Ищет строку заголовка, автоматически определяя разделитель и
    сопоставляя колонки через синонимы.
    Возвращает (индекс_строки_заголовка, разделитель, {канон: индекс}).
    """
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for i, line in enumerate(f):
            raw = line.rstrip("\r\n")
            if not raw.strip():
                continue
            if raw.lstrip().startswith("#"):
                continue
            delimiter = _detect_delimiter(raw)
            tokens = raw.split(delimiter)
            matched = _match_header(tokens)
            if matched:
                logger.info(
                    "MyHeritage: заголовок найден на строке %d, разделитель=%r, поля=%s",
                    i, delimiter, matched,
                )
                return i, delimiter, matched

    supported = "; ".join(
        f"{canon}: {', '.join(sorted(syns))}" for canon, syns in COLUMN_SYNONYMS.items()
    )
    raise MyHeritageFormatError(
        f"Не удалось найти строку заголовка в файле {path}.\n"
        f"Нужно найти минимум {MIN_MATCHED_COLUMNS} из 4 полей "
        f"({', '.join(REQUIRED_COLUMNS)}).\n"
        f"Поддерживаемые варианты названий колонок:\n  {supported}"
    )


def _read_rows(path: Path, header_line: int, delimiter: str) -> Iterator[list[str]]:
    """Читает строки MyHeritage после заголовка с автоопределённым разделителем."""
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for _ in range(header_line + 1):
            next(f, None)
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if not row or all(cell.strip() == "" for cell in row):
                continue
            yield [_strip_quotes(cell) for cell in row]


# ---------------------------------------------------------------------------
# Главная функция парсинга
# ---------------------------------------------------------------------------
def parse_myheritage_v5(
    csv_path: Path,
    reference: ReferenceGenome,
    both_non_ref_threshold_pct: float = DEFAULT_BOTH_NON_REF_THRESHOLD_PCT,
    liftover: Optional[ChainLiftover] = None,
) -> ParseResult:
    """
    Парсит файл MyHeritage (CSV или TSV, гибкий заголовок) и возвращает
    ParseResult. Интерфейс идентичен parse_ftdna_v3, включая параметр
    liftover (промт "HRC / TopMed") — см. докстринг parse_ftdna_v3()/
    adapters/ftdna_v3.py про место и причину переноса координаты.
    """
    csv_path = Path(csv_path)
    header_line, delimiter, columns = _find_header(csv_path)

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing_cols:
        logger.warning(
            "MyHeritage: не найдены колонки %s — строки, где они использовались бы, "
            "будут считаться некорректными (malformed_rows).",
            missing_cols,
        )

    result = ParseResult()
    positions_for_signature: list[tuple[str, int]] = []
    # Задача D: широкая сигнатура — отпечаток ФИЗИЧЕСКОГО ДИЗАЙНА чипа.
    #
    # ⚠ ИСПРАВЛЕНИЕ БАГА (то же самое, что и в adapters/ftdna_v3.py —
    # см. подробный комментарий там): раньше broad_key регистрировался
    # ПОСЛЕ проверки `if genotype in NO_CALL: continue`, из-за чего
    # широкая сигнатура зависела от того, у какого конкретно человека
    # сколько позиций не прочиталось ("--"/"NA"/...) — а это отличается
    # от образца к образцу даже на одном и том же чипе. В результате
    # chip_signature_broad почти никогда не совпадала у разных людей на
    # одном чипе, и переиспользование доноров (Задача D) не срабатывало.
    #
    # Исправление: broad_key регистрируется для каждой успешно
    # распарсенной строки (валидное число полей, валидная позиция),
    # независимо от значения генотипа.
    positions_for_signature_broad: list[tuple[str, int]] = []
    # Некоторые чипы MyHeritage измеряют одну и ту же позицию несколько раз
    # (например, разными зондами) — такие дубликаты нельзя пускать дальше
    # как есть: build_vcf() требует уникальности (chrom,pos) и падает с
    # PureCoreError. Оставляем первое встреченное измерение этой позиции,
    # остальные учитываем в result.duplicate_positions и пропускаем.
    seen_positions: set[tuple[str, int]] = set()
    seen_positions_broad: set[tuple[str, int]] = set()

    rsid_idx = columns.get("RSID")
    chrom_idx = columns.get("CHROMOSOME")
    pos_idx = columns.get("POSITION")
    result_idx = columns.get("RESULT")
    needed_idx = [i for i in (rsid_idx, chrom_idx, pos_idx, result_idx) if i is not None]
    min_fields = (max(needed_idx) + 1) if needed_idx else 0

    for row in _read_rows(csv_path, header_line, delimiter):
        if len(row) < min_fields or None in (rsid_idx, chrom_idx, pos_idx, result_idx):
            result.malformed_rows += 1
            continue

        rsid = row[rsid_idx].strip()
        chrom = row[chrom_idx].strip()
        pos_str = row[pos_idx].strip()
        genotype = row[result_idx].strip().upper()

        # --- Задача D (исправлено): парсим позицию и нормализуем
        # хромосому, регистрируем broad-позицию ДО проверки NO_CALL —
        # см. комментарий выше про причину бага.
        try:
            pos = int(pos_str)
        except ValueError:
            result.malformed_rows += 1
            continue

        chrom = _normalize_chrom(chrom)

        # Промт "HRC / TopMed": перенос координаты GRCh37 -> сборка панели
        # ДО broad-сигнатуры и до reference.base_at() — см. подробное
        # объяснение в adapters/ftdna_v3.py::parse_ftdna_v3(). Позиция,
        # которую не удалось перенести, полностью пропускается и
        # учитывается в result.lift_failed.
        if liftover is not None:
            lifted = liftover.lift(chrom, pos)
            if lifted is None:
                result.lift_failed += 1
                continue
            chrom, pos = lifted

        broad_key = (chrom, pos)
        if broad_key not in seen_positions_broad:
            seen_positions_broad.add(broad_key)
            positions_for_signature_broad.append(broad_key)

        if genotype in NO_CALL:
            result.missing += 1
            continue

        result.total_measured += 1

        ref_base = reference.base_at(chrom, pos)
        if ref_base not in COMPLEMENT:
            result.ref_non_acgt += 1
            continue

        resolved = _resolve_genotype(genotype, ref_base)

        if resolved.reject_reason == REJECT_INVALID:
            result.invalid_codes += 1
            continue
        if resolved.reject_reason == REJECT_SELF_COMPLEMENTARY:
            result.het_self_complementary += 1
            continue
        if resolved.reject_reason == REJECT_BOTH_NON_REF:
            result.both_non_ref += 1
            continue

        pos_key = (chrom, pos)
        if pos_key in seen_positions:
            # Позиция уже встречалась раньше в этом же файле — оставляем
            # первое измерение, это не отбрасывается ни в одну из
            # существующих категорий QC, поэтому отдельный счётчик.
            result.duplicate_positions += 1
            continue
        seen_positions.add(pos_key)

        result.variants.append(ParsedVariant(
            rsid=rsid, chrom=chrom, pos=pos,
            ref=resolved.ref, alt=resolved.alt, gt=resolved.gt,
        ))
        positions_for_signature.append((chrom, pos))

    if result.total_measured == 0:
        raise MyHeritageFormatError(
            "В файле MyHeritage не найдено ни одной калиброванной позиции"
        )

    if result.both_non_ref_pct > both_non_ref_threshold_pct:
        raise StrandQualityError(
            f"both_non_ref = {result.both_non_ref} "
            f"({result.both_non_ref_pct:.2f}%) — выше порога "
            f"{both_non_ref_threshold_pct}%. Вероятно, перепутана цепь."
        )

    positions_for_signature.sort()
    result.chip_signature = _chip_signature(positions_for_signature)

    # Задача D
    positions_for_signature_broad.sort()
    result.signature_positions_broad = positions_for_signature_broad
    result.chip_signature_broad = _chip_signature(positions_for_signature_broad)

    logger.info(
        "MyHeritage v5: разделитель=%r, измерено=%d, пропущено=%d, self_comp=%d, "
        "both_non_ref=%d, invalid=%d, ref_non_acgt=%d, дубликатов_позиций=%d, "
        "битых_строк=%d, lift_failed=%d, signature=%s, signature_broad=%s",
        delimiter, result.total_measured, result.missing,
        result.het_self_complementary, result.both_non_ref,
        result.invalid_codes, result.ref_non_acgt, result.duplicate_positions,
        result.malformed_rows, result.lift_failed,
        result.chip_signature, result.chip_signature_broad,
    )
    if result.duplicate_positions:
        logger.warning(
            "MyHeritage v5: обнаружено %d повторных измерений уже встречавшихся "
            "позиций — оставлено первое измерение каждой позиции, остальные "
            "пропущены (см. result.duplicate_positions).",
            result.duplicate_positions,
        )
    if liftover is not None and result.lift_failed:
        logger.info(
            "MyHeritage v5: лифтовер — не перенесено на целевую сборку: %d позиций",
            result.lift_failed,
        )
    return result
