"""
adapters/ftdna_v3.py
Адаптер FTDNA Family Finder (.csv) -> ParseResult.

=============================================================================
ИЗМЕНЕНИЯ (промт "HRC / TopMed", интеграция лифтовера в адаптеры):
=============================================================================

⚠ Найденный баг: main.py/gui/app.py вызывали parse_ftdna_v3(csv_path,
reference, liftover=liftover), а докстринг main.py утверждал, что "лифт
применяется ВНУТРИ parser_fn, ДО reference.base_at()" — но сама функция
не принимала параметр liftover вовсе. Реальный прогон падал с
    TypeError: parse_ftdna_v3() got an unexpected keyword argument 'liftover'
для ЛЮБОЙ panel (включая HRC, где liftover=None, но именованный аргумент
всё равно передаётся). core/liftover.py::ChainLiftover при этом полностью
рабочий и покрыт тестами — он просто никогда не был подключён здесь.

Исправление: добавлен параметр liftover: Optional[ChainLiftover] = None.
Перенос координаты (GRCh37 -> сборка панели) происходит СРАЗУ после
парсинга pos и нормализации chrom, ДО:
  - регистрации broad_key (widescale-сигнатура чипа, Задача D) — иначе
    chip_signature_broad считалась бы по смешанным координатам разных
    сборок и не была бы детерминированной для одного и того же чипа;
  - reference.base_at(chrom, pos) — который резолвит основание по
    референсу ЦЕЛЕВОЙ сборки (для panel="topmed" это GRCh38), поэтому ему
    обязательно нужна уже перенесённая координата, а не исходная GRCh37.

Позиции, которые лифтовер не смог перенести (ChainLiftover.lift()
вернул None — нет chain-блока для хромосомы, разрыв между блоками
выравнивания, результат вне границ целевой хромосомы), полностью
пропускаются (continue) и учитываются в result.lift_failed — они
физически не существуют в целевой сборке, включать их куда-либо дальше
(в variants, в сигнатуру) не имеет смысла.

Для liftover=None (HRC, или обратная совместимость со старыми вызовами
без этого параметра) поведение не меняется ни на йоту — блок if
liftover is not None полностью пропускается.
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
    from pyfaidx import Fasta
except ImportError as exc:
    raise ImportError("Требуется пакет pyfaidx: pip install pyfaidx") from exc

from .base import ParsedVariant, ParseResult
from core.liftover import ChainLiftover

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
EXPECTED_HEADER = ("RSID", "CHROMOSOME", "POSITION", "RESULT")
NO_CALL = "--"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
SELF_COMPLEMENTARY_PAIRS = (frozenset("AT"), frozenset("CG"))

DEFAULT_BOTH_NON_REF_THRESHOLD_PCT = 0.1

CHROM_NORMALIZE = {
    "XY": "X", "M": "MT",
    "chrX": "X", "chrY": "Y", "chrM": "MT", "chrMT": "MT",
}

REJECT_INVALID = "invalid"
REJECT_SELF_COMPLEMENTARY = "self_complementary"
REJECT_BOTH_NON_REF = "both_non_ref"


# ---------------------------------------------------------------------------
class FTDNAFormatError(ValueError):
    pass

class StrandQualityError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
def _normalize_chrom(chrom: str) -> str:
    if chrom.startswith("chr") and chrom not in CHROM_NORMALIZE:
        chrom = chrom[3:]
    return CHROM_NORMALIZE.get(chrom, chrom)


# ---------------------------------------------------------------------------
class ReferenceGenome:
    def __init__(self, fasta_path: Path):
        if not fasta_path.exists():
            raise FileNotFoundError(f"Референс не найден: {fasta_path}")
        # Задача C: сохраняем путь явно — используется bcftools norm
        # (main.py::_normalize_vcf / main() --normalize), которому нужен
        # путь к .fasta-файлу, а не сам объект pyfaidx.Fasta. Раньше этого
        # атрибута не было, и main() подстраховывался через
        # `reference.fasta_path if hasattr(...) else Path(args.reference)`
        # — при автозагрузке референса (args.reference is None) это падало
        # с TypeError: Path(None). Явный атрибут устраняет и хак, и баг.
        self.fasta_path = fasta_path
        self._fasta = Fasta(str(fasta_path), sequence_always_upper=True)
        # Промт "TopMed/HRC", п.2: GRCh37-референс (HRC) называет контиги
        # без префикса ("1", "2", ..., "X", "MT"), а GRCh38-референс
        # (TopMed) — как правило с префиксом "chr" ("chr1", "chrX"), причём
        # митохондриальный контиг в разных сборках GRCh38 может называться
        # и "chrM", и "chrMT". Строим словарь соответствия ОДИН раз при
        # инициализации (сканирование fasta.keys() — это чтение .fai-
        # индекса, а не повторный проход по самому файлу референса), а не
        # на каждый вызов base_at(), и резолвим по нему без try/except-
        # перебора вариантов на каждой позиции.
        self._contig_aliases = self._build_contig_alias_map()

    def _build_contig_alias_map(self) -> dict[str, str]:
        """
        {каноническое_имя_хромосомы: реальное_имя_контига_в_fasta}.
        Каноническое имя — то, что уже приходит на вход base_at() от
        адаптеров (_normalize_chrom() в ftdna_v3.py/myheritage_v5.py):
        "1".."22", "X", "Y", "MT" — без префикса "chr". Первое совпадение
        побеждает (setdefault) — на практике одна сборка не содержит
        одновременно контиги "1" и "chr1", коллизий не бывает.
        """
        aliases: dict[str, str] = {}
        for real_name in self._fasta.keys():
            body = real_name[3:] if real_name.lower().startswith("chr") else real_name
            canonical = "MT" if body in ("M", "MT") else body
            aliases.setdefault(canonical, real_name)
        return aliases

    def base_at(self, chrom: str, pos: int) -> str:
        contig_name = self._contig_aliases.get(chrom, chrom)
        try:
            seq = self._fasta[contig_name][pos - 1 : pos]
        except KeyError as exc:
            raise FTDNAFormatError(
                f"Хромосома '{chrom}' отсутствует в референсе"
            ) from exc
        except Exception as exc:
            # pyfaidx может бросать разные исключения (например, при позиции
            # за пределами длины контига) в зависимости от версии — ловим
            # широко и превращаем в понятную доменную ошибку.
            raise FTDNAFormatError(
                f"Не удалось прочитать референс на позиции {chrom}:{pos}: {exc}"
            ) from exc

        base = str(seq).upper()
        if len(base) != 1:
            # Некоторые версии pyfaidx не бросают исключение на выходе за
            # пределы контига, а молча возвращают пустую/укороченную
            # последовательность — ловим это явно, а не даём упасть ниже
            # по стеку с непонятной ошибкой.
            raise FTDNAFormatError(
                f"Позиция {chrom}:{pos} вне диапазона референсного контига"
            )
        return base


# ---------------------------------------------------------------------------
def _validate_header(path: Path) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        first_line = f.readline().strip().rstrip("\r")
        header = tuple(part.strip() for part in first_line.split(","))
        if header != EXPECTED_HEADER:
            raise FTDNAFormatError(
                f"Неожиданный заголовок FTDNA-файла.\n"
                f"Ожидалось: {EXPECTED_HEADER}\nПолучено: {header}"
            )


def _read_rows(path: Path) -> Iterator[list[str]]:
    """Отдаёт сырые строки CSV как есть (без валидации числа полей) —
    подсчёт и отбрасывание некорректных строк делает вызывающий код,
    чтобы это было видно в QC-статистике ParseResult, а не терялось молча."""
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        yield from reader


# ---------------------------------------------------------------------------
@dataclass
class _Resolved:
    ref: str
    alt: str
    gt: str
    reject_reason: str | None


def _resolve_genotype(result: str, ref_base: str) -> _Resolved:
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


# ---------------------------------------------------------------------------
def _chip_signature(positions: list[tuple[str, int]]) -> str:
    h = hashlib.sha256()
    for chrom, pos in positions:
        h.update(f"{chrom}:{pos}\n".encode("utf-8"))
    return h.hexdigest()[:16]


def save_position_cache(cache_dir: Path, result: ParseResult) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature}.positions.json"
    if not out_path.exists():
        payload = [(v.chrom, v.pos) for v in result.variants]
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


def save_position_cache_broad(cache_dir: Path, result: ParseResult) -> Path:
    """
    Задача D (опционально, выключено по умолчанию) — см. докстринг
    adapters/base.py. Сохраняет позиции ПО ШИРОКОЙ сигнатуре
    (chip_signature_broad, все позиции физического дизайна чипа —
    см. исправление в parse_ftdna_v3() ниже) — используется только когда
    main.py/GUI явно включили режим переиспользования доноров между
    разными людьми на одном чипе. НЕ использовать этот файл позиций
    вместе со старым флагом bcftools merge -0 — см. предупреждение в
    adapters/base.py.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature_broad}.positions.json"
    if not out_path.exists():
        payload = result.signature_positions_broad
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
def parse_ftdna_v3(
    csv_path: Path,
    reference: ReferenceGenome,
    both_non_ref_threshold_pct: float = DEFAULT_BOTH_NON_REF_THRESHOLD_PCT,
    liftover: Optional[ChainLiftover] = None,
) -> ParseResult:
    """
    liftover (промт "HRC / TopMed", интеграция лифтовера в адаптеры):
    если задан, каждая позиция переносится из GRCh37 в целевую сборку
    ЭТИМ chain-лифтовером СРАЗУ после парсинга (chrom, pos) — до
    reference.base_at() (который для panel="topmed" читает GRCh38-
    референс и должен получать уже перенесённую координату) и до
    регистрации в широкой сигнатуре чипа (Задача D). Позиции, которые
    лифтовер не смог перенести, пропускаются целиком и учитываются в
    result.lift_failed. liftover=None (по умолчанию, HRC) — поведение
    идентично версии без этого параметра.
    """
    csv_path = Path(csv_path)
    _validate_header(csv_path)
    result = ParseResult()
    positions_for_signature: list[tuple[str, int]] = []
    # Задача D: широкая сигнатура — отпечаток ФИЗИЧЕСКОГО ДИЗАЙНА чипа
    # (см. докстринг adapters/base.py).
    #
    # ⚠ ИСПРАВЛЕНИЕ БАГА (переиспользование доноров между людьми никогда
    # не срабатывало): раньше broad_key регистрировался ПОСЛЕ проверки
    # `if genotype == NO_CALL: continue` — то есть в широкую сигнатуру
    # попадали только позиции, которые у ЭТОГО КОНКРЕТНОГО человека не
    # оказались "--". Доля "--" (шум гибридизации/чтения чипа) у каждого
    # человека своя (обычно 0.05-0.3%, но всегда разная), поэтому
    # chip_signature_broad у двух разных людей на ОДНОМ И ТОМ ЖЕ чипе
    # почти никогда не совпадала, и весь смысл Задачи D (не перекачивать
    # доноров повторно для одного и того же чипа) пропадал: проверка
    # sha256(sorted(positions_broad)) просто никогда не давала совпадения.
    #
    # Исправление: broad_key теперь регистрируется для КАЖДОЙ успешно
    # распарсенной строки (валидное число полей, валидная позиция),
    # НЕЗАВИСИМО от того, "--" там или реальный генотип. Это и есть
    # настоящий отпечаток дизайна чипа: тот же набор (chrom, pos) должен
    # быть у любого человека, протестированного на этой же версии чипа,
    # вне зависимости от того, что именно у него прочиталось на каждой
    # позиции.
    positions_for_signature_broad: list[tuple[str, int]] = []
    # Защита от повторных измерений одной и той же позиции в файле (см.
    # аналогичную логику в adapters/myheritage_v5.py) — build_vcf() требует
    # уникальности (chrom,pos), иначе падает PureCoreError.
    seen_positions: set[tuple[str, int]] = set()
    seen_positions_broad: set[tuple[str, int]] = set()

    for row in _read_rows(csv_path):
        if len(row) != 4:
            result.malformed_rows += 1
            continue
        rsid, chrom, pos_str, genotype = (cell.strip() for cell in row)
        genotype = genotype.upper()

        # --- Задача D (исправлено): парсим позицию и нормализуем
        # хромосому, а затем регистрируем broad-позицию ДО проверки
        # NO_CALL — см. комментарий выше про причину бага. Строка,
        # которая не проходит даже базовый парсинг позиции, по-прежнему
        # считается malformed и не участвует ни в одной сигнатуре.
        try:
            pos = int(pos_str)
        except ValueError:
            result.malformed_rows += 1
            continue
        chrom = _normalize_chrom(chrom)

        # Промт "HRC / TopMed": перенос координаты GRCh37 -> сборка панели
        # ДО всего остального (широкой сигнатуры, reference.base_at()) —
        # ниже по коду всё уже должно работать в координатах целевой
        # сборки. Позиция, которую не удалось перенести (нет chain-блока
        # для хромосомы, разрыв между блоками выравнивания, результат вне
        # границ целевого контига), физически не существует в целевой
        # сборке — полностью пропускаем её, не регистрируя ни в одной из
        # сигнатур.
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

        if genotype == NO_CALL:
            result.missing += 1
            continue
        result.total_measured += 1

        ref_base = reference.base_at(chrom, pos)

        if ref_base not in COMPLEMENT:
            # Маскированный/неопределённый референс (например 'N') —
            # позицию нельзя надёжно сравнить с генотипом, отбрасываем.
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
            result.duplicate_positions += 1
            continue
        seen_positions.add(pos_key)

        result.variants.append(ParsedVariant(
            rsid=rsid, chrom=chrom, pos=pos,
            ref=resolved.ref, alt=resolved.alt, gt=resolved.gt,
        ))
        positions_for_signature.append((chrom, pos))

    if result.total_measured == 0:
        raise FTDNAFormatError("В файле не найдено ни одной калированной позиции")
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
        "FTDNA v3: измерено=%d, пропущено=%d, self_comp=%d, "
        "both_non_ref=%d, invalid=%d, ref_non_acgt=%d, дубликатов_позиций=%d, "
        "битых_строк=%d, lift_failed=%d, signature=%s, signature_broad=%s",
        result.total_measured, result.missing,
        result.het_self_complementary, result.both_non_ref,
        result.invalid_codes, result.ref_non_acgt, result.duplicate_positions,
        result.malformed_rows, result.lift_failed,
        result.chip_signature, result.chip_signature_broad,
    )
    if result.duplicate_positions:
        logger.warning(
            "FTDNA v3: обнаружено %d повторных измерений уже встречавшихся "
            "позиций — оставлено первое измерение каждой позиции.",
            result.duplicate_positions,
        )
    if liftover is not None and result.lift_failed:
        logger.info(
            "FTDNA v3: лифтовер — не перенесено на целевую сборку: %d позиций",
            result.lift_failed,
        )
    return result
