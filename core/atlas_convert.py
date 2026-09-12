"""
core/atlas_convert.py
Этап 0 для источника «Атлас»: перенос координат файла Атласа из GRCh38 в
GRCh37 и приведение к оформлению 23andMe v3 (то же, что у трафаретов).

=============================================================================
ЧЕМ ФАЙЛ АТЛАСА ОТЛИЧАЕТСЯ ОТ ОСТАЛЬНЫХ ИСТОЧНИКОВ
=============================================================================

Оформлением — почти ничем: одна '#'-строка с названиями колонок, дальше
четыре колонки через табуляцию (rsid, chromosome, position, genotype),
хромосомы уже как X/Y/MT, гаплоидные вызовы уже одной буквой, пропусков
('--') в файле вообще нет. То есть это тот же формат, что пишет 23andMe,
и adapters/ancestry_v2.py читает его напрямую (LAYOUT_CONVERTED).

Отличие ровно одно, и оно принципиальное: **координаты в GRCh38**, а не в
GRCh37. Проверено на реальном файле (582 151 маркер): из 524 650 общих с
трафаретом v5 rsID позиция совпала лишь у 0,95 %, а после переноса
hg19 -> hg38 совпадение — 99,7 %; сдвиги систематические и разные по
хромосомам (медиана от -2,6 Мб на chr9 до +2,3 Мб на chr18) — характерная
подпись смены сборки. Отдельно подтверждает вывод то, что 14 маркеров в
файле названы не rsID, а HGVS по accession'ам GRCh38 (NC_000006.12,
NC_000010.11 — в GRCh37 это .11 и .10).

Если такой файл отдать пайплайну как есть (а автодетект до этого промта
определял его как MyHeritage), позиции будут сверяться с GRCh37-референсом
и совпадать со случайными местами генома: both_non_ref взлетит, парсер
упадёт на StrandQualityError, а в лучшем случае — молча соберётся мусор.

=============================================================================
ПОЧЕМУ ЭТО ОТДЕЛЬНЫЙ ЭТАП, А НЕ ПАРАМЕТР liftover У ПАРСЕРА
=============================================================================

adapters/ancestry_v2.py умеет лифтовать координаты на лету (параметр
liftover), и технически можно было бы просто передать ему chain-файл
hg38ToHg19. Но параметр liftover уже занят другим смыслом: main.py/
gui/app.py передают туда ПРЯМОЙ лифтовер GRCh37 -> сборка панели (для
panel="topmed"). Смешивать в одном параметре два разных преобразования,
выбирая направление по источнику, — прямая дорога к тому, что для
Атласа + TopMed однажды применится не то направление.

Здесь же файл один раз нормализуется к GRCh37 — и дальше Атлас для всего
пайплайна неотличим от любого другого чипа: трафареты, доноры, HRC — всё
в GRCh37, а для panel="topmed" обычный прямой лифтовер отработает как для
FTDNA/MyHeritage. Двойной перенос (38 -> 37 -> 38) для TopMed обходится
десятыми долями процента позиций, которые chain-файл не переносит.

Плюс тот же бонус, что у AncestryDNA: промежуточный файл — самостоятельный
результат. Это валидный «сырой файл 23andMe v3» в GRCh37, который можно
залить в Генотек как есть, не дожидаясь импутации, и можно проверить
глазами.

=============================================================================
ЧТО ИМЕННО МЕНЯЕТСЯ
=============================================================================

  координаты  -> GRCh38 -> GRCh37 по UCSC chain-файлу hg38ToHg19
                 (core/liftover.py, тот же, которым Этап 7 возвращает
                 результат TopMed в GRCh37)
  шапка       -> '#'-строки трафарета 23andMe v3
  порядок     -> строки пересортировываются по (хромосома, позиция):
                 лифтовер меняет порядок внутри хромосомы (а на
                 перестроенных участках и сильно), а потребители
                 23andMe-формата ждут сортированный файл
  не-rsID     -> строки с HGVS-именами (NC_000006.12:g...) отбрасываются:
                 трафареты адресуются по rsID, такой строке в них всё
                 равно нет соответствия
  дубликаты   -> если после переноса две строки попали в одну (хромосома,
                 позиция), остаётся первая

⚠ Генотипы НЕ меняются. Аллели Атласа даны на plus-strand, как и у
23andMe; смена сборки не меняет цепь. Инделы (I/D) и гаплоидные вызовы
одной буквой переносятся как есть — их разбирает дальше адаптер, с
подсчётом в QC.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from core.ancestry_convert import read_template_header
from core.liftover import ChainLiftover

logger = logging.getLogger(__name__)

#: Суффикс имени сконвертированного файла. Отличается от
#: ancestry_convert.CONVERTED_SUFFIX намеренно: там меняется только
#: оформление, здесь ещё и сборка генома, и это должно быть видно по имени.
CONVERTED_SUFFIX = "_grch37_23andme_v3.txt"

#: Заголовок колонок файла Атласа (он же — заголовок 23andMe).
ATLAS_HEADER = ("rsid", "chromosome", "position", "genotype")

#: Порядок хромосом в выходном файле (как в 23andMe/трафаретах).
CHROM_ORDER = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "MT": 25}

#: Сколько первых строк просматривать в поисках строки с названиями колонок.
MAX_HEADER_SCAN_LINES = 100

#: Максимум ведущих '#'-строк, при котором файл ещё считается «атласовым»
#: по оформлению. У Атласа шапка — ровно одна строка (сама строка
#: колонок), у 23andMe/MyHeritage — десятки; это не запрет, а часть
#: сигнала для автодетекта (см. detect_layout()).
MAX_ATLAS_COMMENT_LINES = 3


class AtlasConvertError(ValueError):
    pass


@dataclass
class AtlasConversionStats:
    """Что получилось на выходе — печатается в лог этапа и в run_info.json."""
    out_path: Optional[Path] = None
    skipped: bool = False           # вход уже был в GRCh37
    rows_in: int = 0
    rows_out: int = 0
    malformed_rows: int = 0
    lift_failed: int = 0
    non_rsid_dropped: int = 0
    duplicate_positions: int = 0
    header_lines: int = 0
    unknown_chroms: list[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.skipped:
            return (
                f"файл уже в GRCh37 и в оформлении 23andMe, "
                f"конвертация не требуется: {self.out_path}"
            )
        return (
            f"{self.rows_in} строк -> {self.rows_out} в "
            f"{self.out_path.name if self.out_path else '?'} "
            f"(GRCh38 -> GRCh37; не перенеслось {self.lift_failed}, "
            f"не-rsID отброшено {self.non_rsid_dropped}, дубликатов позиций "
            f"{self.duplicate_positions}, битых строк {self.malformed_rows})"
        )


# ---------------------------------------------------------------------------
def _header_tokens(line: str) -> tuple[str, ...]:
    return tuple(
        t.strip().strip('"').lower()
        for t in line.lstrip("#").strip().split("\t")
    )


def locate_header(path: Path) -> Optional[tuple[int, int]]:
    """
    (индекс строки с названиями колонок, число ведущих '#'-строк) или None,
    если строку колонок 23andMe/Атласа найти не удалось.
    """
    comments = 0
    with Path(path).open("r", encoding="utf-8-sig", newline="") as f:
        for i, line in enumerate(f):
            if i >= MAX_HEADER_SCAN_LINES:
                break
            raw = line.rstrip("\r\n")
            if not raw.strip():
                continue
            if _header_tokens(raw) == ATLAS_HEADER:
                return i, comments
            if raw.lstrip().startswith("#"):
                comments += 1
                continue
            # Первая же не-комментарийная строка, и она не заголовок —
            # дальше искать бессмысленно.
            break
    return None


def detect_layout(path: Path) -> Optional[str]:
    """
    "atlas" — файл оформлен как 23andMe (4 колонки через табуляцию) и при
    этом шапка короткая (не больше MAX_ATLAS_COMMENT_LINES '#'-строк), как
    у Атласа. None — не похоже.

    Сборку генома здесь НЕ проверяем: это дело автодетекта источника
    (main.py::detect_source_from_file(), где есть трафарет для сверки) —
    задача этой функции только оформление.
    """
    found = locate_header(Path(path))
    if found is None:
        return None
    _, comments = found
    return "atlas" if comments <= MAX_ATLAS_COMMENT_LINES else None


def _iter_rows(path: Path, header_line: int):
    with Path(path).open("r", encoding="utf-8-sig", newline="") as f:
        for i, line in enumerate(f):
            if i <= header_line:
                continue
            raw = line.rstrip("\r\n")
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            yield raw.split("\t")


def convert_atlas_to_grch37(
    src: Path,
    dst: Path,
    liftover: ChainLiftover,
    template_path: Optional[Path] = None,
) -> AtlasConversionStats:
    """
    Переносит файл Атласа в GRCh37 и пишет его в оформлении 23andMe v3.

    src           — файл Атласа (.txt).
    dst           — куда писать результат.
    liftover      — ChainLiftover на hg38ToHg19.over.chain.gz (ОБЯЗАТЕЛЕН:
                    без него смысла в этом этапе нет, и молча пропускать
                    перенос нельзя — получился бы файл с координатами
                    GRCh38 под видом GRCh37).
    template_path — откуда взять '#'-шапку; обычно samples/template_v3.txt.
    """
    src, dst = Path(src), Path(dst)
    found = locate_header(src)
    if found is None:
        raise AtlasConvertError(
            f"В файле {src} не найдена строка колонок "
            f"{'/'.join(ATLAS_HEADER)} — проверьте, что выбран верный "
            f"источник данных."
        )
    header_line, _ = found

    stats = AtlasConversionStats(out_path=dst)
    unknown: set[str] = set()
    rows: list[tuple[int, int, str, str, str]] = []

    for row in _iter_rows(src, header_line):
        stats.rows_in += 1
        if len(row) != 4:
            stats.malformed_rows += 1
            continue
        rsid = row[0].strip()
        chrom = row[1].strip()
        try:
            pos = int(row[2].strip())
        except ValueError:
            stats.malformed_rows += 1
            continue
        genotype = row[3].strip()

        if not rsid.lower().startswith("rs"):
            # HGVS-имена (NC_000006.12:g.18143597T>G) — в трафаретах таким
            # строкам соответствия нет, адресация там по rsID.
            stats.non_rsid_dropped += 1
            continue

        lifted = liftover.lift(chrom, pos)
        if lifted is None:
            stats.lift_failed += 1
            continue
        new_chrom, new_pos = lifted
        new_chrom = new_chrom[3:] if new_chrom.lower().startswith("chr") else new_chrom
        if new_chrom not in CHROM_ORDER:
            unknown.add(new_chrom)
            continue

        rows.append((CHROM_ORDER[new_chrom], new_pos, rsid, new_chrom, genotype))

    if not rows:
        raise AtlasConvertError(
            f"После переноса координат из {src.name} не осталось ни одной "
            f"строки — похоже, chain-файл не подходит к этому файлу."
        )

    rows.sort(key=lambda r: (r[0], r[1]))

    header = read_template_header(template_path)
    stats.header_lines = len(header)
    stats.unknown_chroms = sorted(unknown)

    seen: set[tuple[str, int]] = set()
    dst.parent.mkdir(parents=True, exist_ok=True)
    # newline="" + '\n' в строках: перевод строки у v3 — LF, и Windows не
    # должен превращать его в CRLF (у v5/genotek наоборот, но целевое
    # оформление здесь именно v3).
    with dst.open("w", encoding="utf-8", newline="") as f:
        for line in header:
            f.write(line + "\n")
        for _, pos, rsid, chrom, genotype in rows:
            key = (chrom, pos)
            if key in seen:
                stats.duplicate_positions += 1
                continue
            seen.add(key)
            f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            stats.rows_out += 1

    logger.info("Атлас: %s", stats.summary())
    if stats.unknown_chroms:
        logger.warning(
            "Атлас: после переноса встретились неизвестные контиги (строки "
            "отброшены): %s", ", ".join(stats.unknown_chroms),
        )
    return stats


def prepare_atlas_file(
    src: Path,
    out_dir: Path,
    liftover: ChainLiftover,
    template_path: Optional[Path] = None,
) -> AtlasConversionStats:
    """
    Обёртка для пайплайна (Этап 0): конвертирует src в
    out_dir/<имя>_grch37_23andme_v3.txt и возвращает статистику, в которой
    out_path — файл, который надо отдать парсеру дальше.

    Идемпотентна по имени: если на вход дали файл, который уже прошёл этот
    этап (суффикс CONVERTED_SUFFIX), конвертация пропускается — иначе
    координаты уехали бы в GRCh36. По содержимому «уже ли это GRCh37»
    надёжно судить здесь нельзя (нужен трафарет для сверки, и это работа
    автодетекта источника), поэтому проверка именно по имени файла, а не
    по догадке.
    """
    src, out_dir = Path(src), Path(out_dir)
    if src.name.endswith(CONVERTED_SUFFIX):
        logger.info(
            "Атлас: файл %s уже прошёл перенос в GRCh37 — Этап 0 пропущен.",
            src.name,
        )
        return AtlasConversionStats(out_path=src, skipped=True)

    if detect_layout(src) is None:
        raise AtlasConvertError(
            f"Файл {src} не похож на файл Атласа: ожидались 4 колонки "
            f"через табуляцию ({', '.join(ATLAS_HEADER)}) и короткая шапка. "
            f"Проверьте, что выбран верный источник данных."
        )

    dst = out_dir / (src.stem + CONVERTED_SUFFIX)
    return convert_atlas_to_grch37(src, dst, liftover, template_path=template_path)
