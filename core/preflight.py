"""
core/preflight.py
Отчёт "состояние исходника" — дешёвый анализ сырого файла ДО запуска
пайплайна: один проход по файлу, ни одной закачки, ни референса, ни
доноров.

Зачем. Половина разборов постфактум ("почему Генотек не принял", "почему
на X дыры", "почему гаплогруппы не определились") была бы видна сразу,
если бы кто-то посмотрел на исходник до того, как потратить час на
скачивание доноров и сутки на очередь Michigan Imputation Server.

Чего этот модуль сознательно НЕ делает:

  * не выбирает панель. Качество файла и выбор панели — разные оси.
    Плохой файл (низкий call rate чипа, контаминация, перепутанный пол)
    большая панель не спасёт — она аккуратно доимпутирует мусор. Такой
    файл надо ПОМЕТИТЬ и предупредить, а не молча переключить на TopMed.
    Выбор панели живёт в core/panel_advisor.py и опирается на состав
    чипа по частотам, а не на качество;
  * не отбраковывает и не правит файл. Все проверки — информационные:
    решение продолжать остаётся за человеком. Единственное исключение из
    этого правила по замыслу отсутствует: даже "красный" отчёт не мешает
    нажать "Продолжить";
  * не требует референсного генома. Доля позиций, где ОБА аллеля не
    референсные (прямой признак проблем с ориентацией цепи), считается
    позже, самим адаптером (ParseResult.both_non_ref_pct) — здесь для
    неё нет данных, и врать о ней нечем. Зато доля палиндромных (A/T и
    C/G) гетерозигот считается уже тут: она видна по одному генотипу.

Определение сборки генома. Никаких зашитых по памяти координат реперных
rsID: сравниваем позиции файла с ТРАФАРЕТОМ (samples/template_v3.txt) —
это настоящий экспорт 23andMe в GRCh37, который и так лежит в комплекте
и по которому потом собирается итоговый файл. Совпадают позиции у общих
rsID — GRCh37; не совпадают — файл в какой-то другой сборке (почти
наверняка GRCh38), и его нельзя подавать в пайплайн, рассчитанный на
GRCh37-трафарет, без лифтовера.
"""
from __future__ import annotations

import contextlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Пороги. Все — ориентиры для ПРЕДУПРЕЖДЕНИЙ, а не критерии отбраковки.
# Их придётся калибровать по накопленным запускам (см. core/metrics_log.py):
# на двух файлах границу не вывести, поэтому здесь стоят осознанно широкие
# значения, отсекающие только явную патологию.
# ---------------------------------------------------------------------------
MIN_MARKERS = 100_000          # меньше — это не автосомный чип целиком
LOW_CALL_RATE_PCT = 97.0       # ниже — чип отработал плохо
BAD_CALL_RATE_PCT = 95.0       # ниже — файл сомнителен целиком
# Полоса аутосомной гетерозиготности. Ориентир "25-35 %" верен для чипа,
# набранного из ОБЩИХ SNP; у чипа, насыщенного редкими маркерами, доля
# гетерозигот естественно ниже — на файле FTDNA (13 % маркеров с EUR_AF
# ниже 0,5 %) это 21-22 %, и это норма, а не находка. Поэтому нижняя
# граница здесь заметно шире учебной: цель — поймать явную патологию
# (контаминация задирает гетерозиготность вверх), а не ругаться на
# нормальный файл. Точную границу считать по накопленному
# runs_metrics.csv, а не по двум файлам.
HET_BAND_PCT = (20.0, 38.0)
MALE_X_HET_MAX_PCT = 1.0       # как в core/pure_python_core.py
FEMALE_X_HET_MIN_PCT = 10.0
MALE_Y_CALL_RATE_MIN_PCT = 15.0   # у мужчины Y реально прочитан
FEMALE_Y_CALL_RATE_MAX_PCT = 5.0  # у женщины по Y почти сплошь пропуски
MIN_X_CALLS_FOR_SEX = 200
MIN_Y_MARKERS_FOR_SEX = 100
MIN_MT_MARKERS = 100           # меньше — на mtDNA-гаплогруппу рассчитывать нечего
BUILD_MATCH_PCT = 90.0         # доля совпавших с трафаретом позиций -> GRCh37
BUILD_MISMATCH_PCT = 20.0      # ниже — сборка точно другая

AUTOSOMES = {str(i) for i in range(1, 23)}

#: Сколько rsID выбирается из файла для сверки сборки с трафаретом.
#: Больше не нужно: доля совпадений на 5000 маркерах отличает 100 % от 0 %
#: с колоссальным запасом, а память и время прохода по трафарету экономит.
BUILD_PROBE_SIZE = 5000

#: Минимум сверенных с трафаретом маркеров, при котором вывод о сборке
#: вообще имеет смысл делать.
MIN_BUILD_COMPARISONS = 50

_NO_CALL_CHARS = set("-0N.?")
_BASES = set("ACGT")
_PALINDROMIC = {frozenset("AT"), frozenset("CG")}


# ---------------------------------------------------------------------------
# Чтение сырого файла
# ---------------------------------------------------------------------------
# Карты кодов хромосом РАЗНЫЕ у разных производителей: у AncestryDNA код
# "25" — это псевдоаутосомный регион X, а у MyHeritage тот же "25"
# означает MT. Перепутать их — тихо перебросить десятки позиций X в
# митохондриальный контиг, поэтому карта выбирается по источнику, а не
# берётся "общая".
_CHROM_MAPS: dict[str, dict[str, str]] = {
    "ancestry": {"23": "X", "24": "Y", "25": "X", "26": "MT"},
    "myheritage": {"23": "X", "24": "Y", "25": "MT"},
    "ftdna": {"23": "X", "24": "Y", "25": "MT"},
    "vcf": {},
}
_CHROM_COMMON = {
    "X": "X", "Y": "Y", "XY": "X", "MT": "MT", "M": "MT",
}

_COLUMN_SYNONYMS: dict[str, set[str]] = {
    "RSID": {"rsid", "rs_id", "rs id", "snp_id", "snp id", "marker",
             "markername", "id", "snp"},
    "CHROMOSOME": {"chromosome", "chrom", "chr"},
    "POSITION": {"position", "pos", "bp", "basepair", "base_pair"},
    "RESULT": {"result", "genotype", "gt", "allele", "alleles", "call"},
    "ALLELE1": {"allele1", "allele 1", "allele_1"},
    "ALLELE2": {"allele2", "allele 2", "allele_2"},
}


def normalise_chrom(raw: str, source: str = "ftdna") -> str:
    """Код хромосомы -> каноническое имя ("1".."22", "X", "Y", "MT")."""
    c = str(raw).strip().strip('"').strip("'")
    if c.lower().startswith("chr") and len(c) > 3:
        c = c[3:]
    if c in _CHROM_COMMON:
        return _CHROM_COMMON[c]
    upper = c.upper()
    if upper in _CHROM_COMMON:
        return _CHROM_COMMON[upper]
    return _CHROM_MAPS.get(source, {}).get(c, c)


def _split_row(line: str, delimiter: str) -> list[str]:
    return [cell.strip().strip('"').strip("'") for cell in line.split(delimiter)]


def _find_header(path: Path) -> tuple[int, str, dict[str, int]]:
    """
    (номер строки заголовка, разделитель, {каноническое имя: индекс}).

    Тот же приём, что в adapters/myheritage_v5.py::_find_header(): формат
    23andMe-подобных экспортов у всех производителей "почти одинаковый",
    но заголовок бывает и с ведущим '#', и без него, и через табуляцию, и
    через запятую, и с разным набором колонок (RESULT одной колонкой либо
    ALLELE1/ALLELE2 двумя).
    """
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        for lineno, line in enumerate(f):
            if lineno > 200:
                break
            stripped = line.strip().lstrip("#").strip()
            if not stripped:
                continue
            for delimiter in ("\t", ",", ";"):
                if delimiter not in stripped:
                    continue
                tokens = [t.lower() for t in _split_row(stripped, delimiter)]
                found: dict[str, int] = {}
                for canonical, synonyms in _COLUMN_SYNONYMS.items():
                    for idx, token in enumerate(tokens):
                        if token in synonyms:
                            found[canonical] = idx
                            break
                has_result = "RESULT" in found or (
                    "ALLELE1" in found and "ALLELE2" in found
                )
                if {"RSID", "CHROMOSOME", "POSITION"} <= set(found) and has_result:
                    return lineno, delimiter, found
    raise PreflightError(
        f"В файле {path.name} не нашлась строка заголовка с колонками "
        f"rsid/chromosome/position/genotype — возможно, выбран не тот "
        f"источник данных или файл повреждён."
    )


class PreflightError(RuntimeError):
    """Файл не удалось прочитать настолько, чтобы что-то о нём сказать."""


def iter_raw_calls(path: Path, source: str) -> Iterator[tuple[str, str, int, str]]:
    """
    (rsid, chrom, pos, genotype) для каждой строки данных. Генотип — как в
    файле, в верхнем регистре, слитой строкой ("AG", "--", "0", "II").

    Битые строки просто пропускаются: цель прохода — статистика по файлу,
    а не его валидация; отдельный счётчик битых строк ведёт analyse_file().
    """
    if source == "vcf":
        yield from _iter_vcf_calls(path)
        return

    header_line, delimiter, cols = _find_header(path)
    i_rsid, i_chrom, i_pos = cols["RSID"], cols["CHROMOSOME"], cols["POSITION"]
    i_result = cols.get("RESULT")
    i_a1, i_a2 = cols.get("ALLELE1"), cols.get("ALLELE2")
    need = max(x for x in (i_rsid, i_chrom, i_pos, i_result, i_a1, i_a2)
               if x is not None) + 1

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        for lineno, line in enumerate(f):
            if lineno <= header_line:
                continue
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            row = _split_row(line, delimiter)
            if len(row) < need:
                yield ("", "", -1, "")   # маркер битой строки для счётчика
                continue
            try:
                pos = int(row[i_pos])
            except ValueError:
                yield ("", "", -1, "")
                continue
            if i_result is not None:
                genotype = row[i_result].upper()
            else:
                genotype = (row[i_a1] + row[i_a2]).upper()
            yield (row[i_rsid], normalise_chrom(row[i_chrom], source), pos, genotype)


_GT_SPLIT_RE = re.compile(r"[/|]")


def _iter_vcf_calls(path: Path) -> Iterator[tuple[str, str, int, str]]:
    """
    То же самое для готового VCF. Генотип разворачивается обратно в буквы
    (0/1 + REF/ALT -> "AG"), чтобы вся статистика ниже считалась одним и
    тем же кодом для всех источников. Мультиаллельные и не-SNP строки
    (инделы, символические ALT) пропускаются — чип их не измеряет, и в
    статистике call rate им не место.

    errors="replace" — как и везде, где читается прошедший через bcftools
    файл: служебная строка заголовка на Windows может быть в кодировке
    консоли, и одна такая строка не должна ронять анализ.
    """
    import gzip

    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:  # type: ignore[operator]
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                yield ("", "", -1, "")
                continue
            chrom, pos_s, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                pos = int(pos_s)
            except ValueError:
                yield ("", "", -1, "")
                continue
            alts = alt.split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts if a != "."):
                continue
            alleles = [ref.upper()] + [a.upper() for a in alts]
            gt_field = parts[9].split(":")[0]
            calls = _GT_SPLIT_RE.split(gt_field)
            genotype = ""
            for c in calls:
                if c in (".", ""):
                    genotype += "-"
                    continue
                try:
                    genotype += alleles[int(c)]
                except (ValueError, IndexError):
                    genotype += "-"
            yield (rsid, normalise_chrom(chrom, "vcf"), pos, genotype)


# ---------------------------------------------------------------------------
# Отчёт
# ---------------------------------------------------------------------------
@dataclass
class ChromStats:
    markers: int = 0
    called: int = 0

    @property
    def call_rate(self) -> float:
        return 100.0 * self.called / self.markers if self.markers else 0.0


@dataclass
class Finding:
    """
    Одна строка отчёта. level: "ok" | "warn" | "bad" — цвет и порядок
    показа; detail — то, что видно всегда, hint — что с этим делать.
    """
    level: str
    title: str
    detail: str = ""
    hint: str = ""


@dataclass
class PreflightReport:
    path: str = ""
    source: str = ""
    total_rows: int = 0
    malformed_rows: int = 0
    markers: int = 0
    called: int = 0
    indel_calls: int = 0
    per_chrom: dict[str, ChromStats] = field(default_factory=dict)

    autosomal_called: int = 0
    autosomal_het: int = 0
    palindromic_het: int = 0

    x_nonpar_calls: int = 0
    x_nonpar_het: int = 0
    y_markers: int = 0
    y_called: int = 0
    mt_markers: int = 0
    mt_called: int = 0

    duplicate_rsids: int = 0
    duplicate_positions: int = 0
    unsorted_positions: int = 0

    build: str = "неизвестно"
    build_match_pct: Optional[float] = None
    build_probe_used: int = 0

    findings: list[Finding] = field(default_factory=list)

    # -- производные величины ------------------------------------------------
    @property
    def call_rate(self) -> float:
        return 100.0 * self.called / self.markers if self.markers else 0.0

    @property
    def autosomal_het_pct(self) -> float:
        if not self.autosomal_called:
            return 0.0
        return 100.0 * self.autosomal_het / self.autosomal_called

    @property
    def palindromic_het_pct(self) -> float:
        """
        Доля палиндромных (A/T и C/G) среди ГЕТЕРОЗИГОТНЫХ вызовов.

        Именно среди гетерозиготных, а не среди всех: по гомозиготному
        вызову ("AA") невозможно понять, палиндромный это маркер или нет —
        второго аллеля мы не видим. Гетерозиготный "AT" — видно сразу.
        Отсюда и ожидаемый порядок величины: на типовом чипе это единицы
        процентов от гетерозигот.
        """
        if not self.autosomal_het:
            return 0.0
        return 100.0 * self.palindromic_het / self.autosomal_het

    @property
    def x_het_pct(self) -> float:
        if not self.x_nonpar_calls:
            return 0.0
        return 100.0 * self.x_nonpar_het / self.x_nonpar_calls

    @property
    def y_call_rate(self) -> float:
        return 100.0 * self.y_called / self.y_markers if self.y_markers else 0.0

    @property
    def sex_by_x(self) -> str:
        if self.x_nonpar_calls < MIN_X_CALLS_FOR_SEX:
            return "неизвестно"
        if self.x_het_pct < MALE_X_HET_MAX_PCT:
            return "мужской"
        if self.x_het_pct >= FEMALE_X_HET_MIN_PCT:
            return "женский"
        return "неоднозначно"

    @property
    def sex_by_y(self) -> str:
        if self.y_markers < MIN_Y_MARKERS_FOR_SEX:
            return "неизвестно"
        if self.y_call_rate >= MALE_Y_CALL_RATE_MIN_PCT:
            return "мужской"
        if self.y_call_rate <= FEMALE_Y_CALL_RATE_MAX_PCT:
            return "женский"
        return "неоднозначно"

    @property
    def worst_level(self) -> str:
        levels = {f.level for f in self.findings}
        if "bad" in levels:
            return "bad"
        if "warn" in levels:
            return "warn"
        return "ok"


# ---------------------------------------------------------------------------
# Сборка отчёта
# ---------------------------------------------------------------------------
def _is_no_call(genotype: str) -> bool:
    return not genotype or all(ch in _NO_CALL_CHARS for ch in genotype)


def _is_indel(genotype: str) -> bool:
    return any(ch in ("I", "D") for ch in genotype)


def analyse_file(
    path: Path,
    source: str,
    *,
    template_path: Optional[Path] = None,
    par_regions: Optional[tuple[tuple[int, int], tuple[int, int]]] = None,
) -> PreflightReport:
    """
    Один проход по файлу + (при наличии трафарета) один проход по
    трафарету для определения сборки. Ничего не скачивает.

    par_regions — границы псевдоаутосомных регионов X той сборки, в
    которой файл. По умолчанию берутся GRCh37 из
    core.pure_python_core.PAR_REGIONS_BY_BUILD: гетерозиготность в PAR
    нормальна и у мужчин, и её нельзя мешать в оценку пола.
    """
    path = Path(path)
    report = PreflightReport(path=str(path), source=source)

    if par_regions is None:
        from core.pure_python_core import PAR_REGIONS_BY_BUILD
        par_regions = PAR_REGIONS_BY_BUILD["grch37"]
    par1, par2 = par_regions

    seen_rsids: set[str] = set()
    seen_positions: set[tuple[str, int]] = set()
    last_key: Optional[tuple[str, int]] = None
    # Проба для сверки сборки. Прореживание АДАПТИВНОЕ: пока проба не
    # переполнилась, берём подряд; при переполнении выбрасываем каждую
    # вторую запись и удваиваем шаг. Так проба равномерно размазана по
    # всему файлу (а не утыкается в начало chr1), память ограничена, и
    # это одинаково работает и на файле в 200 строк из теста, и на
    # реальных 750 тысячах — фиксированный шаг "каждый сотый" на первом
    # дал бы всего две записи и молчаливое "сборка неизвестна".
    probe: dict[str, int] = {}
    probe_step = 1
    probe_counter = 0

    for rsid, chrom, pos, genotype in iter_raw_calls(path, source):
        if pos < 0:
            report.malformed_rows += 1
            continue
        report.total_rows += 1
        report.markers += 1

        stats = report.per_chrom.get(chrom)
        if stats is None:
            stats = report.per_chrom[chrom] = ChromStats()
        stats.markers += 1

        # дубли и отсортированность
        if rsid and rsid.startswith("rs"):
            if rsid in seen_rsids:
                report.duplicate_rsids += 1
            else:
                seen_rsids.add(rsid)
                probe_counter += 1
                if probe_counter % probe_step == 0:
                    probe[rsid] = pos
                    if len(probe) >= BUILD_PROBE_SIZE * 2:
                        probe = {k: v for i, (k, v) in enumerate(probe.items())
                                 if i % 2 == 0}
                        probe_step *= 2
        key = (chrom, pos)
        if key in seen_positions:
            report.duplicate_positions += 1
        else:
            seen_positions.add(key)
        if last_key is not None and last_key[0] == chrom and pos < last_key[1]:
            report.unsorted_positions += 1
        last_key = key

        if _is_indel(genotype):
            report.indel_calls += 1
            continue
        if _is_no_call(genotype):
            continue

        report.called += 1
        stats.called += 1
        bases = [b for b in genotype if b in _BASES]

        if chrom in AUTOSOMES:
            report.autosomal_called += 1
            if len(bases) == 2 and bases[0] != bases[1]:
                report.autosomal_het += 1
                if frozenset(bases) in _PALINDROMIC:
                    report.palindromic_het += 1
        elif chrom == "X":
            in_par = (par1[0] <= pos <= par1[1]) or (par2[0] <= pos <= par2[1])
            if not in_par:
                report.x_nonpar_calls += 1
                if len(bases) == 2 and bases[0] != bases[1]:
                    report.x_nonpar_het += 1

    report.y_markers = report.per_chrom.get("Y", ChromStats()).markers
    report.y_called = report.per_chrom.get("Y", ChromStats()).called
    report.mt_markers = report.per_chrom.get("MT", ChromStats()).markers
    report.mt_called = report.per_chrom.get("MT", ChromStats()).called

    if template_path is not None:
        _detect_build(report, probe, Path(template_path))

    _add_findings(report)
    return report


def build_match_pct(probe: dict[str, int],
                    template_path: Path) -> tuple[Optional[float], int]:
    """
    Доля пробных rsID, чья позиция совпала с трафаретом (реальный экспорт
    23andMe в GRCh37), и число сверенных маркеров.

    Отдельная публичная функция, а не часть _detect_build(), потому что тем
    же способом определяет сборку автодетект источника в
    main.py::detect_source_from_file() (файлы Атласа приходят в GRCh38, и
    отличить их от файлов в GRCh37 можно только по координатам —
    оформление у них одинаковое). Возвращает (None, 0), если сверить не с
    чем: нет пробы, нет трафарета, сверено меньше MIN_BUILD_COMPARISONS
    маркеров.
    """
    if not probe or not Path(template_path).is_file():
        return None, 0
    matched = compared = 0
    try:
        with open(template_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t") if "\t" in line else line.split(",")
                if len(parts) < 3:
                    continue
                rsid = parts[0].strip()
                want = probe.get(rsid)
                if want is None:
                    continue
                try:
                    tpos = int(parts[2].strip())
                except ValueError:
                    continue
                compared += 1
                if tpos == want:
                    matched += 1
    except OSError as e:
        logger.info("Определение сборки пропущено: %s", e)
        return None, 0

    if compared < MIN_BUILD_COMPARISONS:
        return None, compared
    return 100.0 * matched / compared, compared


def _detect_build(report: PreflightReport, probe: dict[str, int],
                  template_path: Path) -> None:
    """
    Сверяет позиции пробных rsID с трафаретом (реальный экспорт 23andMe в
    GRCh37). Совпало почти всё — GRCh37; почти ничего — другая сборка.
    Промежуточный результат честно называется неоднозначным, а не
    округляется в удобную сторону.
    """
    pct, compared = build_match_pct(probe, Path(template_path))
    if pct is None:
        return
    report.build_match_pct = pct
    report.build_probe_used = compared
    if pct >= BUILD_MATCH_PCT:
        report.build = "GRCh37"
    elif pct <= BUILD_MISMATCH_PCT:
        report.build = "не GRCh37 (вероятно GRCh38)"
    else:
        report.build = "неоднозначно"


def _add_findings(r: PreflightReport) -> None:
    """Превращает числа в понятные строки отчёта."""
    add = r.findings.append

    # --- объём и общий call rate ---
    if r.markers < MIN_MARKERS:
        add(Finding("bad", "Слишком мало маркеров",
                    f"{r.markers:,} строк".replace(",", " "),
                    "Похоже, это не полный автосомный экспорт чипа — "
                    "проверьте, что выбран нужный файл и он не обрезан."))
    else:
        add(Finding("ok", "Объём файла",
                    f"{r.markers:,} маркеров".replace(",", " ")))

    cr = r.call_rate
    if cr < BAD_CALL_RATE_PCT:
        add(Finding("bad", "Низкий call rate чипа", f"{cr:.2f}%",
                    "Импутация не чинит плохое чтение чипа: большая "
                    "панель просто аккуратно доимпутирует пропуски. "
                    "Если файл можно перезаказать у лаборатории — это "
                    "даст больше, чем любая настройка здесь."))
    elif cr < LOW_CALL_RATE_PCT:
        add(Finding("warn", "Call rate чипа ниже обычного", f"{cr:.2f}%",
                    "Обычно у исправного чипа 98-99,5%. Запуск имеет "
                    "смысл, но итоговая заполняемость будет ниже."))
    else:
        add(Finding("ok", "Call rate чипа", f"{cr:.2f}%"))

    # --- гетерозиготность ---
    het = r.autosomal_het_pct
    lo, hi = HET_BAND_PCT
    if r.autosomal_called < MIN_MARKERS // 2:
        pass
    elif het > hi:
        add(Finding("bad", "Повышенная аутосомная гетерозиготность",
                    f"{het:.1f}% (ожидается {lo:.0f}-{hi:.0f}%)",
                    "Так выглядит контаминация — смесь ДНК двух людей. "
                    "Импутация её не разделит, а результат будет "
                    "правдоподобным на вид и неверным по сути."))
    elif het < lo:
        add(Finding("warn", "Пониженная аутосомная гетерозиготность",
                    f"{het:.1f}% (обычно {lo:.0f}-{hi:.0f}%)",
                    "Либо файл прошёл нестандартную обработку, либо часть "
                    "гетерозигот потеряна при экспорте. Учтите, что у "
                    "чипов, насыщенных редкими маркерами, значение и в "
                    "норме ниже — это ориентир, а не диагноз."))
    else:
        add(Finding("ok", "Аутосомная гетерозиготность", f"{het:.1f}%"))

    # --- пол ---
    by_x, by_y = r.sex_by_x, r.sex_by_y
    detail = (f"по X: {by_x} (гетерозиготность nonPAR {r.x_het_pct:.2f}% "
              f"на {r.x_nonpar_calls:,} позициях); "
              f"по Y: {by_y} (call rate {r.y_call_rate:.1f}% "
              f"на {r.y_markers:,} маркерах)").replace(",", " ")
    if by_x in ("мужской", "женский") and by_y in ("мужской", "женский"):
        if by_x != by_y:
            add(Finding("bad", "Пол по X и по Y не сходится", detail,
                        "Это красный флаг: перепутанный образец, "
                        "контаминация или сбой экспорта. X будет "
                        "импутирован по одной из версий, и она может "
                        "оказаться не той."))
        else:
            add(Finding("ok", f"Пол образца: {by_x}", detail))
    elif by_x in ("мужской", "женский"):
        # Y в файле нет или он почти пустой — это не противоречие, а
        # просто отсутствие перекрёстной проверки. Пайплайн всё равно
        # определяет пол по X (infer_male_from_variants), поэтому и
        # отчёт называет ответ, а не разводит руками.
        add(Finding("ok", f"Пол образца: {by_x} (по X)", detail,
                    "Проверить по Y не по чему — в файле нет "
                    "Y-маркеров. Пайплайн ориентируется на X, как и "
                    "обычно."))
    elif by_y in ("мужской", "женский"):
        add(Finding("warn", f"Пол образца: {by_y} (только по Y)", detail,
                    "По X данных не хватает, поэтому пайплайн запишет X "
                    "диплоидно (как женский) — это безопасный вариант по "
                    "умолчанию, но для мужского образца он хуже "
                    "гаплоидного."))
    else:
        add(Finding("warn", "Пол определить не по чему", detail,
                    "Данных и по X, и по Y в файле мало. X будет "
                    "импутирован как женский (диплоидно) — это "
                    "безопасный вариант по умолчанию."))

    # --- Y и MT: сразу сказать, что гаплогрупп не будет ---
    if r.y_markers < MIN_Y_MARKERS_FOR_SEX:
        add(Finding("warn", "Y-хромосомы в файле практически нет",
                    f"{r.y_markers} маркеров",
                    "Y-гаплогруппу по этому файлу определить не выйдет."))
    if r.mt_markers < MIN_MT_MARKERS:
        add(Finding("warn", "Митохондриальной ДНК в файле практически нет",
                    f"{r.mt_markers} маркеров",
                    "mtDNA-гаплогруппу по этому файлу определить не выйдет."))
    else:
        add(Finding("ok", "Митохондриальная ДНК",
                    f"{r.mt_markers} маркеров, call rate {100.0 * r.mt_called / r.mt_markers:.1f}%"))

    # --- палиндромные маркеры ---
    add(Finding(
        "ok" if r.palindromic_het_pct < 15 else "warn",
        "Палиндромные (A/T и C/G) гетерозиготы",
        f"{r.palindromic_het_pct:.1f}% от всех гетерозигот "
        f"({r.palindromic_het:,} шт.)".replace(",", " "),
        "Это те маркеры, для которых ориентацию цепи нельзя разрешить по "
        "самому генотипу — адаптер отбрасывает их осознанно. Доля выше "
        "обычной означает, что в файле их непропорционально много.",
    ))

    # --- структура файла ---
    if r.duplicate_positions or r.duplicate_rsids:
        add(Finding("warn", "Дубли в файле",
                    f"позиций: {r.duplicate_positions}, "
                    f"rsID: {r.duplicate_rsids}",
                    "Повторы отбрасываются при парсинге (остаётся первое "
                    "вхождение) — на результат влияет мало, но говорит о "
                    "нестандартной обработке файла."))
    if r.unsorted_positions:
        add(Finding("warn", "Файл не отсортирован по позициям",
                    f"{r.unsorted_positions} нарушений порядка",
                    "Пайплайн сортирует сам, это не ошибка — но признак, "
                    "что файл собран не производителем."))
    if r.malformed_rows:
        add(Finding("warn", "Битые строки", f"{r.malformed_rows} шт.",
                    "Пропускаются при чтении."))

    # --- сборка генома ---
    if r.build == "GRCh37":
        add(Finding("ok", "Сборка генома: GRCh37",
                    f"совпало {r.build_match_pct:.1f}% позиций из "
                    f"{r.build_probe_used} сверенных с трафаретом"))
    elif r.build.startswith("не GRCh37"):
        if r.source == "atlas":
            # Для Атласа GRCh38 — норма, а не поломка: координаты
            # переносятся в GRCh37 на Этапе 0 (core/atlas_convert.py).
            # Показывать здесь "✗" значило бы пугать пользователя тем,
            # что программа сама и исправит через несколько секунд.
            add(Finding("ok", "Сборка генома: GRCh38 (ожидаемо для Атласа)",
                        f"совпало всего {r.build_match_pct:.1f}% позиций из "
                        f"{r.build_probe_used} сверенных с трафаретом GRCh37",
                        "Координаты будут перенесены в GRCh37 на Этапе 0 — "
                        "отдельным видимым файлом, который можно проверить "
                        "глазами и залить в Генотек как есть."))
        else:
            add(Finding("bad", f"Сборка генома: {r.build}",
                        f"совпало всего {r.build_match_pct:.1f}% позиций из "
                        f"{r.build_probe_used} сверенных с трафаретом",
                        "Весь пайплайн и трафарет рассчитаны на GRCh37. "
                        "Файл в другой сборке даст почти пустой результат — "
                        "его нужно сначала перевести в GRCh37."))
    elif r.build == "неоднозначно":
        add(Finding("warn", "Сборку генома определить не удалось",
                    f"совпало {r.build_match_pct:.1f}% позиций из "
                    f"{r.build_probe_used} сверенных с трафаретом",
                    "Проверьте в шапке файла, в какой сборке даны "
                    "координаты (у GRCh37 обычно написано build 37)."))


# ---------------------------------------------------------------------------
# Текстовое представление (лог, run.log, письмо об ошибке)
# ---------------------------------------------------------------------------
_LEVEL_MARK = {"ok": "✓", "warn": "⚠", "bad": "✗"}


def format_report(r: PreflightReport) -> str:
    lines = [
        "=" * 70,
        "СОСТОЯНИЕ ИСХОДНИКА (пре-флайт, до запуска пайплайна)",
        "=" * 70,
        f"Файл: {Path(r.path).name}",
        f"Источник: {r.source}",
        "",
    ]
    order = {"bad": 0, "warn": 1, "ok": 2}
    for f in sorted(r.findings, key=lambda x: order.get(x.level, 3)):
        mark = _LEVEL_MARK.get(f.level, "·")
        lines.append(f"{mark} {f.title}" + (f": {f.detail}" if f.detail else ""))
        if f.hint:
            lines.append(f"    {f.hint}")
    lines.append("")
    lines.append("Call rate по хромосомам:")
    from core.pure_python_core import _chrom_sort_key
    for chrom in sorted(r.per_chrom, key=_chrom_sort_key):
        st = r.per_chrom[chrom]
        lines.append(f"  chr{chrom:<3} {st.markers:>8} маркеров   "
                     f"{st.call_rate:6.2f}%")
    lines.append("=" * 70)
    return "\n".join(lines)


def metrics_row(r: PreflightReport) -> dict:
    """
    Плоский набор чисел для core/metrics_log.py — то, из чего потом, после
    десятка запусков, можно будет откалибровать пороги вместо того, чтобы
    угадывать их на двух файлах.
    """
    return {
        "markers": r.markers,
        "call_rate": round(r.call_rate, 4),
        "autosomal_het_pct": round(r.autosomal_het_pct, 4),
        "palindromic_het_pct": round(r.palindromic_het_pct, 4),
        "x_nonpar_calls": r.x_nonpar_calls,
        "x_het_pct": round(r.x_het_pct, 4),
        "y_markers": r.y_markers,
        "y_call_rate": round(r.y_call_rate, 4),
        "mt_markers": r.mt_markers,
        "sex_by_x": r.sex_by_x,
        "sex_by_y": r.sex_by_y,
        "duplicate_positions": r.duplicate_positions,
        "duplicate_rsids": r.duplicate_rsids,
        "unsorted_positions": r.unsorted_positions,
        "malformed_rows": r.malformed_rows,
        "indel_calls": r.indel_calls,
        "build": r.build,
        "build_match_pct": (round(r.build_match_pct, 2)
                            if r.build_match_pct is not None else ""),
        "preflight_level": r.worst_level,
    }
