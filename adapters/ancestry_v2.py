"""
adapters/ancestry_v2.py
Адаптер AncestryDNA -> ParseResult.

=============================================================================
ДВА ПРИНИМАЕМЫХ ОФОРМЛЕНИЯ
=============================================================================

Адаптер читает файл в ЛЮБОМ из двух видов и даёт на них одинаковый
результат (одни и те же варианты, одну и ту же chip_signature):

  1. Сырой экспорт AncestryDNA — пять колонок:

        #AncestryDNA raw data download
        ... (~18 строк комментариев, начинающихся с '#')
        rsid<TAB>chromosome<TAB>position<TAB>allele1<TAB>allele2
        rs3131972<TAB>1<TAB>752721<TAB>G<TAB>G

  2. Результат Этапа 0 (core/ancestry_convert.py) — четыре колонки, в
     оформлении 23andMe v3:

        # ... '#'-шапка 23andMe ...
        # rsid<TAB>chromosome<TAB>position<TAB>genotype
        rs3131972<TAB>1<TAB>752721<TAB>GG

В штатном прогоне пайплайн подсовывает сюда именно (2): Этап 0 сначала
приводит файл к оформлению v3 (чтобы промежуточный файл можно было
проверить глазами и залить в Генотек как есть), а дальше идут обычные
этапы. Умение читать (1) при этом сохранено намеренно — оно нужно самому
Этапу 0 для проверок и позволяет разобрать сырой файл напрямую, минуя
конвертацию.

⚠ Почему конвертированный файл читает ЭТОТ адаптер, а не гибкий
parse_myheritage_v5(), который тоже понимает TSV с '#'-шапкой: 23andMe
пишет гаплоидные контиги ОДНОЙ буквой ('G' вместо 'GG' на Y и MT, а у
мужчин и на X), а _resolve_genotype() во всех адаптерах требует ровно
две буквы и отбрасывает всё остальное в invalid_codes. На реальном файле
это молча потеряло бы 259 измеренных позиций MT/Y. Здесь односимвольный
генотип на гаплоидном контиге разворачивается в гомозиготу.

=============================================================================
ЧТО НОРМАЛИЗУЕТСЯ ПРИ ЧТЕНИИ СЫРОГО ФОРМАТА
=============================================================================

  1. Генотип разложен на ДВЕ колонки (allele1, allele2) — склеиваем.
  2. Хромосомы закодированы числами: 23=X, 24=Y, 25=PAR, 26=MT.
     ⚠ 25 у Ancestry — это ПСЕВДОАУТОСОМНАЯ ОБЛАСТЬ X (PAR), а не
     митохондрия. В adapters/myheritage_v5.py::CHROM_NORMALIZE тот же код
     "25" означает "MT" — это РАЗНЫЕ соглашения разных производителей,
     поэтому карту оттуда переиспользовать нельзя, здесь заведена своя.
     Ошибка в эту сторону тихо перебросила бы ~36 позиций X в
     митохондриальный контиг.
  3. Пропуск обозначается аллелем "0", а не "--".
  4. Строка с названиями колонок идёт БЕЗ ведущего '#' — в отличие от
     23andMe, где она оформлена как комментарий.
  5. Переводы строк CRLF (открываем с newline="" и режем сами).

Координаты — GRCh37 ("build 37.1" в шапке Ancestry), та же сборка, что у
FTDNA/MyHeritage и у трафаретов, поэтому лифтовер для panel="hrc" не
нужен, а для panel="topmed" работает ровно так же, как в ftdna_v3.

Аллели Ancestry даны на forward-strand относительно референса — как и у
FTDNA. Это НЕ повод пропускать разрешение ориентации: тот же
reference.base_at() + _resolve_genotype() ниже выполняют ту же проверку
и дают ту же QC-статистику (both_non_ref), которая ловит перепутанную
цепь у конкретного файла, чем бы она ни была вызвана.

=============================================================================
ЧТО СОЗНАТЕЛЬНО СДЕЛАНО ТАК ЖЕ, КАК В ftdna_v3.py
=============================================================================

  - Инделы (I/D в колонках аллелей — у Ancestry V2 их порядка 17 тысяч)
    отбрасываются как invalid_codes, потому что _resolve_genotype()
    требует две буквы ACGT. Это ровно то же поведение, что у FTDNA/
    MyHeritage: VCF для импутации строится только по SNP.
  - broad_key (широкая сигнатура чипа, Задача D) регистрируется ДО
    проверки на пропуск — иначе сигнатура зависела бы от того, где
    именно у ЭТОГО человека не прочиталось, и переиспользование доноров
    между людьми на одном чипе не срабатывало бы никогда (подробности —
    в докстринге adapters/ftdna_v3.py::parse_ftdna_v3).
  - Повторные измерения одной и той же (chrom, pos) отбрасываются с
    подсчётом в duplicate_positions: build_vcf() требует уникальности,
    иначе падает PureCoreError. Для Ancestry это не теория — в реальном
    файле V2.0 встречается ~650 повторяющихся позиций.
  - liftover применяется сразу после нормализации (chrom, pos), до
    broad_key и до reference.base_at().
"""
from __future__ import annotations
import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

try:
    from pyfaidx import Fasta  # noqa: F401  (нужен для ReferenceGenome ниже)
except ImportError as exc:
    raise ImportError("Требуется пакет pyfaidx: pip install pyfaidx") from exc

from .base import ParsedVariant, ParseResult
from .ftdna_v3 import ReferenceGenome, StrandQualityError
from core.liftover import ChainLiftover
# Единая точка правды для строк заголовков: их же использует Этап 0
# (core/ancestry_convert.py), и расхождение между модулями означало бы,
# что конвертер пишет файл, который адаптер не узнаёт.
from core.ancestry_convert import (
    ANCESTRY_HEADER as EXPECTED_HEADER,
    TWENTYTHREE_HEADER as CONVERTED_HEADER,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Константы
# ---------------------------------------------------------------------------
#: Оформления, которые понимает адаптер (см. докстринг модуля).
LAYOUT_RAW = "ancestry"        # 5 колонок, сырой экспорт AncestryDNA
LAYOUT_CONVERTED = "23andme"   # 4 колонки, результат Этапа 0

_COLUMNS_BY_LAYOUT = {LAYOUT_RAW: 5, LAYOUT_CONVERTED: 4}

#: Аллель-пропуск в сыром файле AncestryDNA (аналог "--" у 23andMe).
NO_CALL_ALLELE = "0"
#: Как пропуск выглядит в конвертированном файле (23andMe пишет "--";
#: остальное — терпимость к чужим/ручным правкам).
NO_CALL_23ANDME = ("--", "", "0", "00", "NA", "N/A")

#: Контиги, на которых 23andMe пишет гаплоидный генотип одной буквой.
HAPLOID_CHROMS = ("X", "Y", "MT")

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
SELF_COMPLEMENTARY_PAIRS = (frozenset("AT"), frozenset("CG"))

DEFAULT_BOTH_NON_REF_THRESHOLD_PCT = 0.1

#: Коды хромосом AncestryDNA. 23=X, 24=Y, 25=PAR (псевдоаутосомная
#: область X — сливается с X, отдельного контига под неё нет ни в
#: референсе, ни в трафаретах), 26=MT.
#: ⚠ НЕ переиспользовать CHROM_NORMALIZE из myheritage_v5.py: там 25=MT.
CHROM_NORMALIZE = {
    "23": "X", "24": "Y", "25": "X", "26": "MT",
    "X": "X", "Y": "Y", "XY": "X", "MT": "MT", "M": "MT",
    "chrX": "X", "chrY": "Y", "chrM": "MT", "chrMT": "MT",
}

REJECT_INVALID = "invalid"
REJECT_SELF_COMPLEMENTARY = "self_complementary"
REJECT_BOTH_NON_REF = "both_non_ref"

#: Сколько первых строк просматривать в поисках строки с названиями
#: колонок. У реальных файлов Ancestry шапка — 18 строк, у 23andMe — до
#: 24; запас взят на случай, если шапку расширят в новых версиях.
MAX_HEADER_SCAN_LINES = 100


# ---------------------------------------------------------------------------
# Исключения
# ---------------------------------------------------------------------------
class AncestryFormatError(ValueError):
    pass


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------
def _normalize_chrom(chrom: str) -> str:
    """Нормализует хромосому AncestryDNA: 23->X, 24->Y, 25->X (PAR), 26->MT."""
    c = chrom.strip().strip('"').strip("'")
    if c.lower().startswith("chr") and c not in CHROM_NORMALIZE:
        c = c[3:]
    return CHROM_NORMALIZE.get(c, c)


def _find_header(path: Path) -> tuple[int, str]:
    """
    Ищет строку с названиями колонок и возвращает (индекс_строки, layout).

    layout — LAYOUT_RAW или LAYOUT_CONVERTED, см. докстринг модуля.
    Ведущий '#' и регистр при сравнении игнорируются: у Ancestry строка
    колонок идёт без '#', у 23andMe — с ним, и оба варианта допустимы.
    """
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for i, line in enumerate(f):
            if i >= MAX_HEADER_SCAN_LINES:
                break
            raw = line.rstrip("\r\n")
            if not raw.strip():
                continue
            candidate = raw.lstrip("#").strip()
            tokens = tuple(t.strip().strip('"').lower() for t in candidate.split("\t"))
            if tokens == EXPECTED_HEADER:
                logger.info("AncestryDNA: сырой формат, заголовок на строке %d", i)
                return i, LAYOUT_RAW
            if tokens == CONVERTED_HEADER:
                logger.info("AncestryDNA: формат 23andMe, заголовок на строке %d", i)
                return i, LAYOUT_CONVERTED
            if not raw.lstrip().startswith("#"):
                # Первая же не-комментарийная строка и она не заголовок —
                # дальше искать бессмысленно, это данные другого формата.
                raise AncestryFormatError(
                    f"Неожиданный заголовок на строке {i} файла {path}.\n"
                    f"Ожидалось {EXPECTED_HEADER} (сырой AncestryDNA) "
                    f"или {CONVERTED_HEADER} (формат 23andMe).\n"
                    f"Получено: {tokens}"
                )
    raise AncestryFormatError(
        f"В первых {MAX_HEADER_SCAN_LINES} строках файла {path} не найдена "
        f"строка с названиями колонок: ни {EXPECTED_HEADER}, ни {CONVERTED_HEADER}."
    )


def _read_rows(path: Path, header_line: int) -> Iterator[list[str]]:
    """
    Отдаёт строки данных как есть (без валидации числа полей) — подсчёт и
    отбрасывание некорректных строк делает вызывающий код, чтобы это было
    видно в QC-статистике ParseResult, а не терялось молча.

    csv.reader здесь не используется намеренно: файл строго TSV без
    кавычек и экранирования, а split("\\t") на 677 тысячах строк заметно
    дешевле.
    """
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for _ in range(header_line + 1):
            next(f, None)
        for line in f:
            raw = line.rstrip("\r\n")
            if not raw.strip():
                continue
            yield raw.split("\t")


def _genotype_from_row(layout: str, row: list[str], chrom: str) -> Optional[str]:
    """
    Достаёт генотип из строки в том оформлении, в котором она пришла.
    Возвращает None, если это пропуск.

    Односимвольный генотип на гаплоидном контиге ('G' на MT) разворачи-
    вается в гомозиготу ('GG') — так его и понимает _resolve_genotype(),
    и так он записан в сыром файле Ancestry, где гаплоидность не
    кодируется вовсе. На аутосоме одна буква — не гаплоидность, а
    испорченная строка: разворачивать её нельзя, она уйдёт в
    invalid_codes на общих основаниях.
    """
    if layout == LAYOUT_RAW:
        a1, a2 = row[3].strip().upper(), row[4].strip().upper()
        if a1 == NO_CALL_ALLELE or a2 == NO_CALL_ALLELE:
            return None
        return a1 + a2

    genotype = row[3].strip().upper()
    if genotype in NO_CALL_23ANDME:
        return None
    if len(genotype) == 1 and chrom in HAPLOID_CHROMS:
        return genotype + genotype
    return genotype


# ---------------------------------------------------------------------------
@dataclass
class _Resolved:
    ref: str
    alt: str
    gt: str
    reject_reason: str | None


def _resolve_genotype(result: str, ref_base: str) -> _Resolved:
    """
    Разрешает генотип относительно референса. Логика в точности та же,
    что в ftdna_v3.py/myheritage_v5.py — дублируется, а не импортируется,
    по тому же соглашению, что уже принято между этими двумя модулями.
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


# ---------------------------------------------------------------------------
def _chip_signature(positions: list[tuple[str, int]]) -> str:
    h = hashlib.sha256()
    for chrom, pos in positions:
        h.update(f"{chrom}:{pos}\n".encode("utf-8"))
    return h.hexdigest()[:16]


def save_position_cache(cache_dir: Path, result: ParseResult) -> Path:
    """
    Идентична save_position_cache() в ftdna_v3.py/myheritage_v5.py/
    vcf_source.py — единый интерфейс для
    main.py::SOURCES[...]["save_position_cache"].
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature}.positions.json"
    if not out_path.exists():
        payload = [(v.chrom, v.pos) for v in result.variants]
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


def save_position_cache_broad(cache_dir: Path, result: ParseResult) -> Path:
    """Задача D — аналог save_position_cache() по широкой сигнатуре чипа."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature_broad}.positions.json"
    if not out_path.exists():
        payload = result.signature_positions_broad
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# Главная функция парсинга
# ---------------------------------------------------------------------------
def parse_ancestry_v2(
    csv_path: Path,
    reference: ReferenceGenome,
    both_non_ref_threshold_pct: float = DEFAULT_BOTH_NON_REF_THRESHOLD_PCT,
    liftover: Optional[ChainLiftover] = None,
) -> ParseResult:
    """
    Парсит данные AncestryDNA — сырой экспорт или результат Этапа 0 в
    оформлении 23andMe (оформление определяется автоматически по строке
    заголовка) — и возвращает ParseResult.

    Сигнатура и семантика полностью совпадают с parse_ftdna_v3()/
    parse_myheritage_v5(), включая параметр liftover — см. докстринг
    parse_ftdna_v3() про место и причину переноса координаты.
    """
    csv_path = Path(csv_path)
    header_line, layout = _find_header(csv_path)
    expected_columns = _COLUMNS_BY_LAYOUT[layout]

    result = ParseResult()
    positions_for_signature: list[tuple[str, int]] = []
    positions_for_signature_broad: list[tuple[str, int]] = []
    seen_positions: set[tuple[str, int]] = set()
    seen_positions_broad: set[tuple[str, int]] = set()

    for row in _read_rows(csv_path, header_line):
        if len(row) != expected_columns:
            result.malformed_rows += 1
            continue

        rsid = row[0].strip()
        try:
            pos = int(row[2].strip())
        except ValueError:
            result.malformed_rows += 1
            continue
        chrom = _normalize_chrom(row[1])

        if liftover is not None:
            lifted = liftover.lift(chrom, pos)
            if lifted is None:
                result.lift_failed += 1
                continue
            chrom, pos = lifted

        # Широкая сигнатура (Задача D) — ДО проверки на пропуск, см.
        # докстринг модуля.
        broad_key = (chrom, pos)
        if broad_key not in seen_positions_broad:
            seen_positions_broad.add(broad_key)
            positions_for_signature_broad.append(broad_key)

        genotype = _genotype_from_row(layout, row, chrom)
        if genotype is None:
            result.missing += 1
            continue
        result.total_measured += 1

        ref_base = reference.base_at(chrom, pos)
        if ref_base not in COMPLEMENT:
            # Маскированный/неопределённый референс ('N') — позицию нельзя
            # надёжно сравнить с генотипом.
            result.ref_non_acgt += 1
            continue

        resolved = _resolve_genotype(genotype, ref_base)

        if resolved.reject_reason == REJECT_INVALID:
            # Сюда же попадают инделы I/D — см. докстринг модуля.
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
        raise AncestryFormatError("В файле не найдено ни одной калированной позиции")
    if result.both_non_ref_pct > both_non_ref_threshold_pct:
        raise StrandQualityError(
            f"both_non_ref = {result.both_non_ref} "
            f"({result.both_non_ref_pct:.2f}%) — выше порога "
            f"{both_non_ref_threshold_pct}%. Вероятно, перепутана цепь."
        )

    positions_for_signature.sort()
    result.chip_signature = _chip_signature(positions_for_signature)

    positions_for_signature_broad.sort()
    result.signature_positions_broad = positions_for_signature_broad
    result.chip_signature_broad = _chip_signature(positions_for_signature_broad)

    logger.info(
        "AncestryDNA (%s): измерено=%d, пропущено=%d, self_comp=%d, "
        "both_non_ref=%d, invalid(вкл. инделы I/D)=%d, ref_non_acgt=%d, "
        "дубликатов_позиций=%d, битых_строк=%d, lift_failed=%d, "
        "signature=%s, signature_broad=%s",
        layout,
        result.total_measured, result.missing,
        result.het_self_complementary, result.both_non_ref,
        result.invalid_codes, result.ref_non_acgt, result.duplicate_positions,
        result.malformed_rows, result.lift_failed,
        result.chip_signature, result.chip_signature_broad,
    )
    if result.duplicate_positions:
        logger.warning(
            "AncestryDNA: обнаружено %d повторных измерений уже встречавшихся "
            "позиций — оставлено первое измерение каждой позиции.",
            result.duplicate_positions,
        )
    if liftover is not None and result.lift_failed:
        logger.info(
            "AncestryDNA: лифтовер — не перенесено на целевую сборку: %d позиций",
            result.lift_failed,
        )
    return result
