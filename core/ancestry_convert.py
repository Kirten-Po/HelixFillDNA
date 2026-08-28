"""
core/ancestry_convert.py
Этап 0 для источника AncestryDNA: приведение сырого экспорта Ancestry к
ОФОРМЛЕНИЮ 23andMe v3 — тому же, что у трафарета samples/template_v3.txt.

=============================================================================
ЗАЧЕМ ЭТО ОТДЕЛЬНЫЙ ШАГ, А НЕ ЧАСТЬ ПАРСЕРА
=============================================================================

adapters/ancestry_v2.py и так умеет читать сырой пятиколоночный Ancestry
напрямую, поэтому для самого пайплайна конвертация не обязательна. Она
сделана отдельным, видимым шагом по двум причинам:

  1. Промежуточный файл — самостоятельный результат. Генотек (и другие
     сервисы, принимающие «сырые данные 23andMe») отказывается принимать
     файл Ancestry именно из-за оформления: пять колонок вместо четырёх,
     хромосомы числами 23-26 вместо X/Y/MT, пропуск как аллель '0'
     вместо '--'. Конвертированный файл можно залить туда как есть, не
     дожидаясь импутации.
  2. Его видно и можно проверить глазами. Нормализация «на лету» внутри
     парсера ничего не оставляет на диске, и если что-то пойдёт не так
     (например, чужой формат под видом Ancestry), разбираться придётся
     по логам, а не по файлу.

⚠ Конвертация НЕ меняет генотипы. Аллели Ancestry уже даны на forward-
strand в GRCh37 — той же цепи и той же сборке, что у трафаретов, поэтому
менять здесь нечего; преобразуется только оформление. Разрешение
ориентации по референсу происходит дальше, в адаптере, как и для любого
другого источника.

=============================================================================
ЧТО ИМЕННО МЕНЯЕТСЯ
=============================================================================

  шапка       -> '#'-строки трафарета 23andMe v3 (или встроенная, если
                 трафарет не передан/не читается)
  колонки     -> 4 (rsid, chromosome, position, genotype): allele1 и
                 allele2 склеиваются
  хромосомы   -> 23->X, 24->Y, 25->X (PAR!), 26->MT
                 ⚠ 25 у Ancestry — псевдоаутосомная область X, а НЕ
                 митохондрия (у MyHeritage тот же код означает MT).
  пропуск     -> аллель '0' превращается в '--'
  гаплоиды    -> Y и MT (а у мужчин и X) 23andMe пишет ОДНОЙ буквой:
                 'G', а не 'GG'. Пол определяется по гетерозиготности X
                 в самом файле, а не спрашивается у пользователя.
  порядок     -> строки пересортировываются: у Ancestry PAR (код 25) идёт
                 ПОСЛЕ Y, и после слияния с X он обязан встать на своё
                 место внутри X, иначе файл не отсортирован по (хромосома,
                 позиция), как ожидают потребители 23andMe-формата.
  переводы    -> LF (как в v3; в v5 — CRLF, но целевой формат здесь
                 строк      именно v3)

Дубликаты позиций и инделы (I/D) НЕ трогаются: они есть и в настоящих
файлах 23andMe, а отбрасывает их дальше адаптер, с подсчётом в QC.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

#: Суффикс имени конвертированного файла.
CONVERTED_SUFFIX = "_23andme_v3.txt"

#: Заголовок колонок сырого файла AncestryDNA (без ведущего '#').
ANCESTRY_HEADER = ("rsid", "chromosome", "position", "allele1", "allele2")

#: Заголовок колонок 23andMe (в самих файлах оформлен как комментарий).
TWENTYTHREE_HEADER = ("rsid", "chromosome", "position", "genotype")

#: 23=X, 24=Y, 25=PAR (часть X), 26=MT.
CHROM_MAP = {"23": "X", "24": "Y", "25": "X", "26": "MT"}

#: Порядок хромосом в выходном файле.
CHROM_ORDER = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "MT": 25}

#: Контиги, которые 23andMe пишет одной буквой, когда они гаплоидны.
HAPLOID_CHROMS = ("Y", "MT")

NO_CALL_ALLELE = "0"
NO_CALL_23ANDME = "--"

#: Порог гетерозиготности X (в процентах от прочитанных позиций X), ниже
#: которого образец считается мужским. У женщин доля гетерозигот по X
#: порядка 20-25 %, у мужчин — доли процента (только шум и PAR), так что
#: порог в 1 % разделяет их с огромным запасом и не требует спрашивать
#: пол у пользователя.
MALE_X_HET_THRESHOLD_PCT = 1.0

#: Встроенная шапка на случай, если трафарет не передан или нечитаем.
#: Намеренно НЕ выдаёт себя за настоящий экспорт 23andMe: строки про
#: сборку генома и TAB-разделение нужны потребителям формата, а вот
#: происхождение файла честно указано.
_FALLBACK_HEADER = (
    "# This data file was converted from AncestryDNA raw data by HelixFillDNA.",
    "#",
    "# Below is a text version of your data.  Fields are TAB-separated",
    "# Each line corresponds to a single SNP.  For each SNP, we provide its identifier",
    "# (an rsid or an internal id), its location on the reference human genome, and the",
    "# genotype call oriented with respect to the plus strand on the human reference sequence.",
    "# We are using reference human assembly build 37 (also known as Annotation Release 104).",
    "#",
    "# rsid\tchromosome\tposition\tgenotype",
)


class AncestryConvertError(ValueError):
    pass


@dataclass
class ConversionStats:
    """Что получилось на выходе — печатается в лог этапа и в run_info.json."""
    out_path: Optional[Path] = None
    skipped: bool = False           # вход уже был в формате 23andMe
    rows: int = 0
    malformed_rows: int = 0
    no_call: int = 0
    haploid_single_letter: int = 0
    par_merged_into_x: int = 0
    male: bool = False
    x_het_pct: float = 0.0
    header_lines: int = 0
    unknown_chroms: list[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.skipped:
            return f"файл уже в формате 23andMe, конвертация не требуется: {self.out_path}"
        return (
            f"{self.rows} строк -> {self.out_path.name if self.out_path else '?'} "
            f"(пропусков {self.no_call}, гаплоидных одной буквой "
            f"{self.haploid_single_letter}, PAR слито с X {self.par_merged_into_x}, "
            f"битых строк {self.malformed_rows}, пол по X "
            f"{'мужской' if self.male else 'женский'} при гетерозиготности "
            f"{self.x_het_pct:.1f}%)"
        )


# ---------------------------------------------------------------------------
def _header_tokens(line: str) -> tuple[str, ...]:
    return tuple(t.strip().strip('"').lower() for t in line.lstrip("#").strip().split("\t"))


def _locate_header(path: Path, max_scan_lines: int = 100) -> Optional[tuple[int, str]]:
    """
    Ищет строку с названиями колонок и возвращает (индекс_строки,
    оформление) либо None, если ни одно оформление не опознано.

    Индекс нужен не только для диагностики: строку заголовка сырого
    Ancestry ('rsid<TAB>chromosome<TAB>...') нельзя отличить от строки
    данных по одному лишь количеству полей — в ней тоже пять колонок и
    нет ведущего '#'. Без явного пропуска всего, что до заголовка
    включительно, конвертер принимал бы её за испорченную строку данных
    (int("position") -> ValueError) и накручивал malformed_rows.

    Заголовок 23andMe оформлен как комментарий ('# rsid\\tchromosome...'),
    заголовок Ancestry — без '#', поэтому ведущий '#' при сравнении
    игнорируется, а отличают форматы сами названия колонок.
    """
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as f:
            for i, line in enumerate(f):
                if i >= max_scan_lines:
                    break
                raw = line.rstrip("\r\n")
                if not raw.strip():
                    continue
                tokens = _header_tokens(raw)
                if tokens == ANCESTRY_HEADER:
                    return i, "ancestry"
                if tokens == TWENTYTHREE_HEADER:
                    return i, "23andme"
                if not raw.lstrip().startswith("#"):
                    # Первая не-комментарийная строка и она не заголовок —
                    # это уже данные, дальше искать нечего.
                    return None
    except OSError as e:
        logger.warning("Не удалось прочитать %s для определения формата: %s", path, e)
        return None
    return None


def detect_layout(path: Path, max_scan_lines: int = 100) -> Optional[str]:
    """
    Оформление файла: "ancestry" (сырой пятиколоночный экспорт),
    "23andme" (четырёхколоночный формат 23andMe, в т.ч. результат
    convert_ancestry_to_23andme_v3()) или None.
    """
    located = _locate_header(path, max_scan_lines=max_scan_lines)
    return located[1] if located else None


def _read_template_header(template_path: Optional[Path]) -> list[str]:
    """'#'-строки трафарета. При любой проблеме — встроенная шапка."""
    if template_path is None:
        return list(_FALLBACK_HEADER)
    try:
        header: list[str] = []
        with Path(template_path).open("r", encoding="utf-8-sig", newline="") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                header.append(line.rstrip("\r\n"))
        if not header:
            raise ValueError("в трафарете нет '#'-строк шапки")
        return header
    except (OSError, ValueError) as e:
        logger.warning(
            "Шапку из трафарета %s взять не удалось (%s) — использую встроенную.",
            template_path, e,
        )
        return list(_FALLBACK_HEADER)


# ---------------------------------------------------------------------------
def convert_ancestry_to_23andme_v3(
    src: Path,
    dst: Path,
    template_path: Optional[Path] = None,
) -> ConversionStats:
    """
    Переоформляет сырой файл AncestryDNA в вид 23andMe v3.

    src            — сырой экспорт AncestryDNA (.txt).
    dst            — куда писать результат.
    template_path  — откуда взять '#'-шапку; обычно samples/template_v3.txt.
                     None или нечитаемый файл — используется встроенная.

    Генотипы не изменяются, см. докстринг модуля.
    """
    src, dst = Path(src), Path(dst)
    located = _locate_header(src)
    layout = located[1] if located else None
    if layout != "ancestry":
        raise AncestryConvertError(
            f"Файл {src} не похож на сырой экспорт AncestryDNA: не найдена "
            f"строка с названиями колонок {ANCESTRY_HEADER}."
        )

    header = _read_template_header(template_path)
    stats = ConversionStats(out_path=dst, header_lines=len(header))

    header_line = located[0]
    rows: list[tuple[str, str, int, str, str]] = []
    unknown: set[str] = set()
    with src.open("r", encoding="utf-8-sig", newline="") as f:
        # Пропускаем шапку вместе со строкой названий колонок — см.
        # докстринг _locate_header() про то, почему её нельзя отличить
        # от строки данных по содержимому.
        for _ in range(header_line + 1):
            next(f, None)
        for line in f:
            s = line.rstrip("\r\n")
            if not s.strip() or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) != 5:
                stats.malformed_rows += 1
                continue
            rsid, chrom_raw, pos_str, a1, a2 = (p.strip() for p in parts)
            try:
                pos = int(pos_str)
            except ValueError:
                stats.malformed_rows += 1
                continue
            chrom = CHROM_MAP.get(chrom_raw, chrom_raw)
            if chrom_raw == "25":
                stats.par_merged_into_x += 1
            if chrom not in CHROM_ORDER:
                unknown.add(chrom_raw)
            rows.append((rsid, chrom, pos, a1.upper(), a2.upper()))

    stats.rows = len(rows)
    stats.unknown_chroms = sorted(unknown)
    if not rows:
        raise AncestryConvertError(f"В файле {src} не найдено ни одной строки данных.")
    if unknown:
        logger.warning(
            "AncestryDNA: неизвестные коды хромосом %s — строки перенесены как есть "
            "и окажутся в конце файла.", stats.unknown_chroms,
        )

    # --- Пол по гетерозиготности X (нужен только чтобы решить, писать X
    # одной буквой или двумя) ---------------------------------------------
    x_called = x_het = 0
    for _rsid, chrom, _pos, a1, a2 in rows:
        if chrom != "X" or NO_CALL_ALLELE in (a1, a2):
            continue
        x_called += 1
        if a1 != a2:
            x_het += 1
    stats.x_het_pct = (100.0 * x_het / x_called) if x_called else 0.0
    stats.male = x_called > 0 and stats.x_het_pct < MALE_X_HET_THRESHOLD_PCT
    haploid = set(HAPLOID_CHROMS) | ({"X"} if stats.male else set())

    # --- Сборка строк ------------------------------------------------------
    out: list[tuple[int, int, str, str, str]] = []
    for rsid, chrom, pos, a1, a2 in rows:
        if a1 == NO_CALL_ALLELE or a2 == NO_CALL_ALLELE:
            genotype = NO_CALL_23ANDME
            stats.no_call += 1
        elif chrom in haploid and a1 == a2:
            genotype = a1
            stats.haploid_single_letter += 1
        else:
            genotype = a1 + a2
        out.append((CHROM_ORDER.get(chrom, 99), pos, rsid, chrom, genotype))

    out.sort(key=lambda r: (r[0], r[1]))

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="utf-8", newline="\n") as f:
        for h in header:
            f.write(h + "\n")
        for _order, pos, rsid, chrom, genotype in out:
            f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")

    logger.info("AncestryDNA -> 23andMe v3: %s", stats.summary())
    return stats


# ---------------------------------------------------------------------------
def prepare_ancestry_file(
    src: Path,
    out_dir: Path,
    template_path: Optional[Path] = None,
) -> ConversionStats:
    """
    Обёртка для пайплайна (Этап 0): конвертирует src в
    out_dir/<имя>_23andme_v3.txt и возвращает статистику, в которой
    out_path — файл, который надо отдать парсеру дальше.

    Идемпотентна: если на вход дали УЖЕ конвертированный файл (или любой
    другой файл в формате 23andMe), конвертация пропускается и out_path
    указывает на исходный файл. Это не редкий случай — пользователь может
    выбрать источник «AncestryDNA» и подсунуть результат прошлого запуска,
    и переконвертировать его было бы нечем: пяти колонок там уже нет.
    """
    src, out_dir = Path(src), Path(out_dir)
    layout = detect_layout(src)

    if layout == "23andme":
        logger.info(
            "AncestryDNA: файл %s уже в формате 23andMe — Этап 0 пропущен.", src.name,
        )
        return ConversionStats(out_path=src, skipped=True)

    if layout is None:
        raise AncestryConvertError(
            f"Файл {src} не похож ни на сырой экспорт AncestryDNA "
            f"(колонки {ANCESTRY_HEADER}), ни на формат 23andMe "
            f"(колонки {TWENTYTHREE_HEADER}). Проверьте, что выбран верный "
            f"источник данных."
        )

    dst = out_dir / (src.stem + CONVERTED_SUFFIX)
    return convert_ancestry_to_23andme_v3(src, dst, template_path=template_path)
