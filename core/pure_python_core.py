"""
core/pure_python_core.py
Чистая Python-реализация общего ядра конвертации.
Базовая функциональность работает без внешних бинарников.
Опционально: если переданы пути к bgzip/tabix/bcftools — используются они.

Заменяет:
  build_vcf()          — генерация VCF из ParseResult
  make_common_positions() — список позиций доноров
  intersect_with_donors() — фильтрация по пересечению
  merge_with_donors()  — объединение sample + donors (через bcftools merge)
  split_autosomes()    — разбивка по хромосомам (в BGZF)
  qc_imputed_vcf()     — QC по Rsq

=============================================================================
ИЗМЕНЕНИЯ (промт "UnicodeDecodeError на нелатинских путях/именах запуска"):
=============================================================================

Предыстория: при имени запуска на кириллице (например, папка
output/runs/мама/) split_autosomes() падал с
    UnicodeDecodeError: 'utf-8' codec can't decode byte 0xec ...
Причина не в самих данных генотипов (они всегда чистый ASCII), а в том, что
bcftools при вызовах norm/view дописывает в СЛУЖЕБНЫЙ заголовок VCF команду,
которой он был вызван (##bcftools_normCommand=..., ##bcftools_viewCommand=...),
включая полный путь к файлу. На Windows эта строка формируется в кодировке
консоли процесса (обычно cp1251 для кириллицы, а не UTF-8) и попадает в
файл как есть — набор байт, невалидный как UTF-8. Чтение файла строгим
`encoding="utf-8"` (без errors=) роняет всю сборку из-за одной
диагностической строки заголовка, которая никак не влияет на сами данные.

Исправление: все места, где ЧИТАЕТСЯ VCF/info-файл, потенциально прошедший
через внешний bcftools (или пришедший с сервера Michigan Imputation
Server), теперь используют errors="replace" — невалидные байты в
служебных строках заголовка (или где угодно ещё) заменяются на U+FFFD
вместо падения программы; сами координаты/генотипы (чистый ASCII) не
затрагиваются. Файлы, которые пишет и сразу же сам читает ЭТОТ модуль (наш
собственный common_pos.txt и т.п., где мы гарантированно контролируем
кодировку на записи), errors= не получили — там строгая проверка полезнее
как признак реальной порчи данных, а не побочный эффект стороннего
инструмента.
"""
from __future__ import annotations
import gzip
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from adapters.base import ParsedVariant, ParseResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Структуры данных
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class QCResult:
    total: int
    retained: int
    rejected: int
    threshold: float


class PureCoreError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------
def _normalise_chrom(chrom: str) -> str:
    """Нормализует хромосому: убирает 'chr', приводит к единому формату."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def _vcf_gt(gt: str) -> str:
    """Приводит GT к VCF-формату (0/0, 0/1, 1/1)."""
    gt = gt.strip()
    allowed = {"0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1"}
    if gt not in allowed:
        raise PureCoreError(f"Неподдерживаемый GT: {gt!r}")
    return gt.replace("|", "/")


# ---------------------------------------------------------------------------
# X-хромосома: псевдоаутосомные регионы и определение пола
#
# Промт "Покрытие X-хромосомы" (жалоба Генотека: ~30% пропусков на X).
# Michigan Imputation Server делит присланный chrX на PAR1/nonPAR/PAR2 сам,
# но требует, чтобы В ПРЕДЕЛАХ nonPAR у каждого образца была ОДНА плоидность
# ("Ploidy Check: verifies if all variants in the nonPAR region are either
# haploid or diploid"). Для мужчины биологически верный и ожидаемый сервером
# вариант — гаплоидный nonPAR (GT "0"/"1"), в PAR — диплоидный, как у
# аутосом. Именно так устроены и донорские мужские образцы 1000 Genomes,
# с которыми наш sample.vcf.gz потом объединяется.
#
# Границы PAR берутся по сборке: чип FTDNA/23andMe помечает часть PAR как
# отдельную "хромосому" XY, адаптеры сводят её к X (см. adapters/*.py
# CHROM_MAP), поэтому опираться на исходную метку нельзя — только на
# координату.
PAR_REGIONS_BY_BUILD: dict[str, tuple[tuple[int, int], tuple[int, int]]] = {
    # (PAR1_start, PAR1_end), (PAR2_start, PAR2_end) — включительно
    "grch37": ((60001, 2699520), (154931044, 155260560)),
    "grch38": ((10001, 2781479), (155701383, 156030895)),
}

# Доля гетерозиготных вызовов в nonPAR X, ниже которой образец считается
# мужским. То же пороговое значение, что и в core/ancestry_convert.py
# (MALE_X_HET_THRESHOLD_PCT) — у мужчин это чистый шум чипа (доли процента),
# у женщин счёт идёт на десятки процентов, промежуточных значений на
# практике не бывает.
MALE_X_HET_THRESHOLD_PCT = 1.0
# Меньше этого числа калиброванных позиций в nonPAR X — определять пол не по
# чему (например, файл вообще без X): считаем образец женским, то есть
# пишем X диплоидно. Это безопасный вариант по умолчанию: диплоидный nonPAR
# сервер тоже принимает, просто импутирует его как женский.
MIN_X_CALLS_FOR_SEX = 200


def is_par_position(pos: int, genome_build: str = "grch37") -> bool:
    """Попадает ли координата X в псевдоаутосомный регион этой сборки."""
    par1, par2 = PAR_REGIONS_BY_BUILD.get(
        genome_build, PAR_REGIONS_BY_BUILD["grch37"],
    )
    return par1[0] <= pos <= par1[1] or par2[0] <= pos <= par2[1]


def infer_male_from_variants(variants, genome_build: str = "grch37") -> tuple[bool, float, int]:
    """
    Определяет пол по гетерозиготности nonPAR X.

    Возвращает (male, het_pct, x_nonpar_calls). male=False при недостатке
    данных (см. MIN_X_CALLS_FOR_SEX) — это осознанно безопасный вариант
    по умолчанию, а не утверждение, что образец женский.
    """
    called = het = 0
    for v in variants:
        if _normalise_chrom(v.chrom) != "X":
            continue
        if is_par_position(int(v.pos), genome_build):
            continue
        called += 1
        if _vcf_gt(v.gt) in ("0/1", "1/0"):
            het += 1
    het_pct = (100.0 * het / called) if called else 0.0
    male = called >= MIN_X_CALLS_FOR_SEX and het_pct < MALE_X_HET_THRESHOLD_PCT
    return male, het_pct, called


def _chrom_sort_key(chrom: str) -> tuple[int, int]:
    """Ключ сортировки для хромосом: 1-22, X, Y, MT."""
    c = _normalise_chrom(chrom)
    if c.isdigit():
        return (int(c), 0)
    if c == "X":
        return (23, 0)
    if c == "Y":
        return (24, 0)
    if c in ("MT", "M"):
        return (25, 0)
    return (99, 0)


# ---------------------------------------------------------------------------
# Валидация
# ---------------------------------------------------------------------------
def validate_variants(result: ParseResult) -> None:
    """Проверяет контракт адаптера перед построением VCF."""
    seen: set[tuple[str, int]] = set()
    allowed_chroms = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    allowed_bases = {"A", "C", "G", "T"}

    for v in result.variants:
        chrom = _normalise_chrom(v.chrom)
        if chrom not in allowed_chroms:
            raise PureCoreError(f"Неподдерживаемая хромосома: {v.chrom}")
        if int(v.pos) <= 0:
            raise PureCoreError(f"Некорректная позиция: {v.chrom}:{v.pos}")
        if v.ref not in allowed_bases:
            raise PureCoreError(f"Некорректный REF: {v.ref!r}")
        if v.alt != "." and v.alt not in allowed_bases:
            raise PureCoreError(f"Некорректный ALT: {v.alt!r}")
        _vcf_gt(v.gt)
        key = (chrom, int(v.pos))
        if key in seen:
            raise PureCoreError(f"Дублирующая позиция в ParseResult: {chrom}:{v.pos}")
        seen.add(key)


# ---------------------------------------------------------------------------
# ParseResult -> VCF
# ---------------------------------------------------------------------------
def build_vcf(
    result: ParseResult,
    output_vcf: Path,
    sample_name: str = "sample",
    compress: bool = True,
    bgzip_path: Optional[str] = None,
    chrom_prefix: str = "",
    haploid_x: bool = False,
    genome_build: str = "grch37",
) -> Path:
    """
    Записывает VCF из ParseResult.

    Режимы сжатия:
      - compress=False → обычный .vcf без сжатия.
      - bgzip_path задан → настоящий BGZF через внешний бинарник (пригоден
        для tabix/bcftools).
      - bgzip_path не задан, compress=True → встроенный gzip (валидный gzip,
        но НЕ BGZF — tabix/bcftools такой файл не примут).

    chrom_prefix (промт "HRC / TopMed"): подставляется ТОЛЬКО в момент
    записи строки в выходной VCF (см. _write_vcf_line()) — "" для HRC/GRCh37
    (поведение не меняется), "chr" для TopMed/GRCh38 (REFERENCE_PANELS[panel]
    ["chrom_prefix"] в main.py), чтобы CHROM-колонка совпадала с тем, как
    называет контиги референс/доноры этой сборки. Сортировка (_chrom_sort_key)
    и валидация (validate_variants) по-прежнему работают с каноническим
    именем хромосомы (без префикса) через _normalise_chrom() — эта функция
    не трогается и не должна получать chrom_prefix.

    haploid_x (промт "Покрытие X-хромосомы"): если True (мужской образец,
    см. infer_male_from_variants()), позиции nonPAR X пишутся ГАПЛОИДНО
    (GT "0"/"1") — этого ждёт Michigan Imputation Server от мужского X
    (Ploidy Check), так же устроены донорские мужские образцы 1000
    Genomes, и только так сервер импутирует X как мужской, а не выдаёт
    биологически невозможные гетерозиготы. PAR остаётся диплоидным.

    Редкие гетерозиготные вызовы в nonPAR X у мужчины (шум гибридизации
    чипа — обычно единицы позиций из десятков тысяч) в гаплоидный GT не
    переводятся однозначно и ОТБРАСЫВАЮТСЯ: оставить их диплоидными
    нельзя (смешанная плоидность в nonPAR — прямой провал Ploidy Check и
    отказ всего задания), а выбирать за чип один из двух аллелей —
    выдумывать данные. Их число логируется.

    genome_build определяет границы PAR (см. PAR_REGIONS_BY_BUILD).
    """
    validate_variants(result)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    variants = sorted(
        result.variants,
        key=lambda v: (_chrom_sort_key(v.chrom), int(v.pos)),
    )

    haploid_positions: set[int] = set()
    if haploid_x:
        dropped_het = 0
        kept: list = []
        for v in variants:
            if (_normalise_chrom(v.chrom) == "X"
                    and not is_par_position(int(v.pos), genome_build)):
                if _vcf_gt(v.gt) in ("0/1", "1/0"):
                    dropped_het += 1
                    continue
                haploid_positions.add(int(v.pos))
            kept.append(v)
        variants = kept
        logger.info(
            "X-хромосома: мужской образец — %d позиций nonPAR записаны "
            "гаплоидно, %d гетерозиготных вызовов nonPAR отброшено "
            "(шум чипа, см. докстринг build_vcf)",
            len(haploid_positions), dropped_het,
        )

    if not compress:
        with output_vcf.open("w", encoding="utf-8", newline="\n") as f:
            _write_vcf_header(f, sample_name)
            for v in variants:
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix,
                                haploid_positions=haploid_positions)
        final_path = output_vcf
    elif bgzip_path:
        # Настоящий BGZF через внешний бинарник
        tmp_vcf = output_vcf.with_suffix(".vcf")
        with tmp_vcf.open("w", encoding="utf-8", newline="\n") as f:
            _write_vcf_header(f, sample_name)
            for v in variants:
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix,
                                haploid_positions=haploid_positions)
        result_proc = subprocess.run(
            [bgzip_path, "-f", str(tmp_vcf)],
            capture_output=True, text=True,
        )
        if result_proc.returncode != 0:
            tmp_vcf.unlink(missing_ok=True)
            raise PureCoreError(f"bgzip завершился с ошибкой:\n{result_proc.stderr}")
        tmp_vcf.unlink(missing_ok=True)
        compressed = Path(str(tmp_vcf) + ".gz")
        compressed.replace(output_vcf)
        final_path = output_vcf
    else:
        # Встроенный gzip (не BGZF!)
        with gzip.open(output_vcf, "wt", encoding="utf-8", newline="\n") as f:
            _write_vcf_header(f, sample_name)
            for v in variants:
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix,
                                haploid_positions=haploid_positions)
        final_path = output_vcf

    logger.info("VCF собран: %d вариантов в %s (chrom_prefix=%r)",
                len(variants), final_path, chrom_prefix)
    return final_path


def _write_vcf_header(f, sample_name: str) -> None:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=PurePythonCore\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")


def _write_vcf_line(f, v: ParsedVariant, chrom_prefix: str = "",
                    haploid_positions: Optional[set] = None) -> None:
    """
    chrom_prefix подставляется ЗДЕСЬ, в момент записи, и только здесь —
    v.chrom/_normalise_chrom() остаются каноническими (без префикса) везде
    в остальном модуле (сортировка, dict-ключи, валидация).

    haploid_positions — координаты X, которые нужно записать гаплоидно
    (мужской nonPAR, см. build_vcf(haploid_x=...)). Набор уже отфильтрован
    вызывающим кодом по хромосоме, поэтому здесь достаточно сверки
    "X + позиция в наборе".
    """
    canonical = _normalise_chrom(v.chrom)
    chrom = f"{chrom_prefix}{canonical}"
    gt = _vcf_gt(v.gt)
    if haploid_positions and canonical == "X" and int(v.pos) in haploid_positions:
        # "0/0" -> "0", "1/1" -> "1"; гетерозиготы сюда не попадают —
        # build_vcf() их отбрасывает до записи.
        gt = gt.split("/")[0]
    f.write(
        f"{chrom}\t{int(v.pos)}\t{v.rsid}\t{v.ref}\t{v.alt}"
        f"\t.\tPASS\t.\tGT\t{gt}\n"
    )


# ---------------------------------------------------------------------------
# Донорские позиции (чистый Python)
# ---------------------------------------------------------------------------
def make_common_positions(
    donor_vcfs: Iterable[Path],
    output_positions: Path,
) -> Path:
    """Извлекает все позиции из донорских VCF."""
    positions: set[tuple[str, int]] = set()
    for donor in donor_vcfs:
        donor = Path(donor)
        if not donor.exists():
            logger.warning("Файл донора не найден: %s", donor)
            continue
        opener = gzip.open if str(donor).endswith(".gz") else open
        # errors="replace" — донорский VCF мог пройти через bcftools,
        # который дописывает в заголовок команду вызова с полным путём
        # (см. докстринг модуля); нас интересуют только строки с
        # координатами, они всегда чистый ASCII.
        with opener(donor, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    chrom = _normalise_chrom(parts[0])
                    try:
                        pos = int(parts[1])
                        positions.add((chrom, pos))
                    except ValueError:
                        continue

    output_positions = Path(output_positions)
    output_positions.parent.mkdir(parents=True, exist_ok=True)
    with output_positions.open("w", encoding="utf-8", newline="\n") as f:
        for chrom, pos in sorted(positions, key=lambda x: (_chrom_sort_key(x[0]), x[1])):
            f.write(f"{chrom}\t{pos}\n")

    logger.info("Извлечено %d уникальных позиций из доноров", len(positions))
    return output_positions


# ---------------------------------------------------------------------------
# Пересечение с донорами (чистый Python)
# ---------------------------------------------------------------------------
def intersect_with_donors(
    sample_vcf: Path,
    common_positions: Path,
    output_vcf: Path,
) -> Path:
    """Оставляет в sample VCF только позиции из common_positions."""
    allowed: set[tuple[str, int]] = set()
    with common_positions.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                chrom, pos = parts
                allowed.add((_normalise_chrom(chrom), int(pos)))
    logger.info("Загружено %d разрешённых позиций", len(allowed))

    sample_vcf = Path(sample_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(sample_vcf).endswith(".gz") else open
    out_opener = gzip.open if str(output_vcf).endswith(".gz") else open

    retained = 0
    # errors="replace" на чтении — sample_vcf мог пройти через bcftools
    # (merge/norm), который дописывает в заголовок команду вызова с полным
    # путём (см. докстринг модуля про кириллицу/не-ASCII в именах запусков).
    # Запись (fout) — наш собственный вывод, кодировку контролируем сами,
    # errors= там не нужен.
    with opener(sample_vcf, "rt", encoding="utf-8", errors="replace") as fin, \
         out_opener(output_vcf, "wt", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                chrom = _normalise_chrom(parts[0])
                try:
                    pos = int(parts[1])
                    if (chrom, pos) in allowed:
                        fout.write(line)
                        retained += 1
                except ValueError:
                    continue

    logger.info("После фильтрации: %d позиций (было %d)", retained, len(allowed))
    return output_vcf


# ---------------------------------------------------------------------------
# Merge с донорами (через bcftools merge)
# ---------------------------------------------------------------------------
def merge_with_donors(
    sample_vcf: Path,
    donor_vcf: Path,
    output_vcf: Path,
    sample_name: str = "sample",
    bcftools_path: Optional[str] = None,
    tabix_path: Optional[str] = None,
) -> Path:
    """
    Объединяет sample VCF с донорским VCF через bcftools merge.
    Это правильно нормализует аллели и обрабатывает несовпадения REF/ALT,
    что критично для Michigan Imputation Server (иначе получаем 400k+
    "Invalid alleles").

    Оба входных файла должны быть BGZF-сжатыми и проиндексированными
    (с .tbi). Это обеспечивается на стороне main.py: build_vcf вызывается
    с bgzip_path, после чего вызывается _index_vcf().
    """
    bcftools = bcftools_path or shutil.which("bcftools")
    tabix = tabix_path or shutil.which("tabix")

    if not bcftools:
        raise PureCoreError(
            "bcftools не найден. Укажите --bin-dir или добавьте в PATH."
        )

    sample_vcf = Path(sample_vcf)
    donor_vcf = Path(donor_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Проверяем наличие индексов
    sample_tbi = Path(str(sample_vcf) + ".tbi")
    donor_tbi = Path(str(donor_vcf) + ".tbi")
    if not sample_tbi.exists():
        raise PureCoreError(
            f"Не найден индекс для {sample_vcf.name}. "
            f"Запустите: tabix -p vcf {sample_vcf}"
        )
    if not donor_tbi.exists():
        raise PureCoreError(
            f"Не найден индекс для {donor_vcf.name}. "
            f"Запустите: tabix -p vcf {donor_vcf}"
        )

    cmd = [
        bcftools, "merge",
        "--force-samples",
        "-0",  # missing генотипы как "./."
        str(sample_vcf),
        str(donor_vcf),
        "-Oz", "-o", str(output_vcf),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise PureCoreError(f"bcftools merge failed:\n{result.stderr}")

    # Индексация выходного файла
    if tabix:
        subprocess.run(
            [tabix, "-p", "vcf", "-f", str(output_vcf)],
            check=True, capture_output=True,
        )

    logger.info("✓ Merged VCF создан: %s", output_vcf)
    return output_vcf


# ---------------------------------------------------------------------------
# Разбивка по хромосомам (BGZF)
# ---------------------------------------------------------------------------
def _write_bgzf(lines: Iterable[str], out_path: Path, bgzip_path: Optional[str] = None) -> None:
    """
    Пишет строки в НАСТОЯЩИЙ BGZF-файл.
    Приоритет:
      1. bgzip_path — внешний бинарник.
      2. Bio.bgzf (biopython) — чистый Python fallback.
    """
    if bgzip_path:
        tmp = out_path.with_suffix(out_path.suffix + ".tmp")
        proc = subprocess.Popen(
            [bgzip_path, "-c"],
            stdin=subprocess.PIPE,
            stdout=open(tmp, "wb"),
        )
        assert proc.stdin is not None
        try:
            for line in lines:
                proc.stdin.write(line.encode("utf-8"))
        finally:
            proc.stdin.close()
            ret = proc.wait()
        if ret != 0:
            tmp.unlink(missing_ok=True)
            raise PureCoreError(f"bgzip завершился с ошибкой (код {ret}) при записи {out_path.name}")
        tmp.replace(out_path)
        return

    try:
        from Bio import bgzf  # type: ignore
    except ImportError as exc:
        raise PureCoreError(
            f"Не удалось записать BGZF для {out_path.name}: не передан bgzip_path "
            f"и не установлен biopython. Установите: pip install biopython"
        ) from exc

    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with bgzf.BgzfWriter(str(tmp), "wb") as f:
        for line in lines:
            f.write(line.encode("utf-8"))
    tmp.replace(out_path)


def _reprefix_chrom_field(line: str, chrom_prefix: str) -> str:
    """
    Заменяет значение поля CHROM (первая колонка строки данных VCF) на
    "{chrom_prefix}{каноническое_имя}" — используется split_autosomes()
    при chrom_prefix != "" (TopMed/GRCh38), чтобы CHROM в выходных
    chr*.vcf.gz совпадал с тем, как называет контиги референс/доноры
    этой сборки. Входная строка уже прошла через merge/bcftools и может
    содержать любой префикс (или отсутствие такового) во входном merged_vcf
    — здесь он нормализуется через _normalise_chrom(), затем подставляется
    заново.
    """
    parts = line.split("\t", 1)
    if len(parts) != 2:
        return line
    canonical = _normalise_chrom(parts[0])
    return f"{chrom_prefix}{canonical}\t{parts[1]}"


# Хромосомы, которые уходят на сервер импутации. X добавлена вместе с
# поддержкой импутации X-хромосомы (жалоба Генотека на ~30% пропусков на
# X): и HRC r1.1, и 1000G Phase 3 v5 на Michigan Imputation Server
# поддерживают X, сервер сам делит присланный chrX.vcf.gz на PAR1/nonPAR/
# PAR2, импутирует их независимо и возвращает одним файлом — от нас
# требуется только НЕ выбрасывать X на этом шаге (раньше выбрасывалась:
# by_chrom содержал только "1".."22", и все ~31 тыс. позиций X молча
# исчезали между batch_merged.vcf.gz и upload/).
#
# Y и MT сюда не входят: этих хромосом нет ни в одной поддерживаемой
# панели импутации, они попадают в финальный файл только как прямые
# измерения чипа (если они есть в исходном файле).
UPLOAD_CHROMS: list[str] = [str(i) for i in range(1, 23)] + ["X"]


_MISSING_GT = {"./.", ".|.", ".", "./.|.", ""}


def _sample_gt(field: str) -> str:
    """GT — всегда первая подколонка поля образца (FORMAT=GT:... )."""
    return field.split(":", 1)[0]


def _replace_gt(field: str, new_gt: str) -> str:
    parts = field.split(":", 1)
    return new_gt if len(parts) == 1 else f"{new_gt}:{parts[1]}"


def normalise_x_ploidy(lines: list[str], genome_build: str = "grch37") -> tuple[list[str], int]:
    """
    Приводит плоидность каждого образца в nonPAR X к его собственной
    преобладающей плоидности. Возвращает (строки, число_исправленных_полей).

    Зачем (найдено живым прогоном на Michigan Imputation Server, задание
    провалилось на QC):

        Error: ChrX nonPAR region includes ambiguous samples (haploid and
        diploid positions). Imputation cannot be started!

    Сервер считает ПРОПУСК "./." диплоидной записью. В нашем мужском
    образце nonPAR писался гаплоидно (см. build_vcf(haploid_x=...)), но
    `bcftools merge` подставлял "./." на позициях, которые есть у доноров
    и отсутствуют на чипе — и этих 13 записей из 28 900 хватило, чтобы
    весь образец стал "ambiguous" и задание было отвергнуто целиком.
    Донорские образцы 1000 Genomes при этом безупречны: каждый либо
    целиком гаплоиден (мужчина), либо целиком диплоиден (женщина).

    Поэтому плоидность нормализуется ПОСЛЕ merge, на готовых строках, и
    для КАЖДОГО образца отдельно — по факту его собственных вызовов, а не
    по нашему предположению о его поле:
      * преобладающая плоидность считается только по непропущенным GT;
      * пропуск переписывается в ту же плоидность ("./." -> "." у
        гаплоидного образца);
      * гомозиготный диплоидный вызов у гаплоидного образца сжимается
        ("0/0" -> "0");
      * гетерозиготный вызов у гаплоидного образца гаплоидным быть не
        может — становится пропуском "." (выбирать за прибор один из двух
        аллелей значило бы выдумывать данные);
      * гаплоидный вызов у диплоидного образца удваивается ("1" -> "1/1").
    PAR не трогается вовсе: там две копии у всех, и проверка сервера на
    него не распространяется.
    """
    data_idx = [i for i, line in enumerate(lines) if not line.startswith("#")]
    if not data_idx:
        return lines, 0

    # --- проход 1: преобладающая плоидность каждого образца в nonPAR ---
    hap_counts: dict[int, int] = {}
    dip_counts: dict[int, int] = {}
    for i in data_idx:
        parts = lines[i].rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        if is_par_position(int(parts[1]), genome_build):
            continue
        for col in range(9, len(parts)):
            gt = _sample_gt(parts[col])
            if gt in _MISSING_GT:
                continue
            if "/" in gt or "|" in gt:
                dip_counts[col] = dip_counts.get(col, 0) + 1
            else:
                hap_counts[col] = hap_counts.get(col, 0) + 1

    haploid_cols = {
        col for col in set(hap_counts) | set(dip_counts)
        if hap_counts.get(col, 0) > dip_counts.get(col, 0)
    }
    if not haploid_cols and not dip_counts:
        return lines, 0

    # --- проход 2: правка ---
    fixed = 0
    out = list(lines)
    for i in data_idx:
        raw = out[i]
        newline_suffix = "\n" if raw.endswith("\n") else ""
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        if is_par_position(int(parts[1]), genome_build):
            continue
        changed = False
        for col in range(9, len(parts)):
            gt = _sample_gt(parts[col])
            want_haploid = col in haploid_cols
            new_gt = None
            if want_haploid:
                if gt in _MISSING_GT:
                    if gt != ".":
                        new_gt = "."
                elif "/" in gt or "|" in gt:
                    a, b = gt.replace("|", "/").split("/", 1)
                    new_gt = a if a == b else "."
            else:
                if gt == ".":
                    new_gt = "./."
                elif gt not in _MISSING_GT and "/" not in gt and "|" not in gt:
                    new_gt = f"{gt}/{gt}"
            if new_gt is not None and new_gt != gt:
                parts[col] = _replace_gt(parts[col], new_gt)
                changed = True
                fixed += 1
        if changed:
            out[i] = "\t".join(parts) + newline_suffix
    return out, fixed


def split_autosomes(
    merged_vcf: Path,
    output_dir: Path,
    bgzip_path: Optional[str] = None,
    chrom_prefix: str = "",
    chroms: Optional[list[str]] = None,
    genome_build: str = "grch37",
) -> list[Path]:
    """
    Делит merged VCF на chr1..chr22 + chrX (см. UPLOAD_CHROMS).
    Историческое имя функции (split_autosomes) сохранено ради обратной
    совместимости с вызывающим кодом и тестами; chroms= позволяет сузить
    набор явно (например, вернуть поведение "только аутосомы").
    Выходные файлы пишутся как настоящий BGZF — они загружаются на
    Michigan Imputation Server, который проверяет формат через tabix.

    chrom_prefix (промт "HRC / TopMed"): "" для HRC/GRCh37 (поведение не
    меняется — CHROM пишется как есть после _normalise_chrom()), "chr" для
    TopMed/GRCh38 (REFERENCE_PANELS[panel]["chrom_prefix"] в main.py).
    Группировка by_chrom по-прежнему ведётся по КАНОНИЧЕСКОМУ имени
    ("1".."22", без префикса) — это только ключи словаря/имена выходных
    файлов (chr{chrom}.vcf.gz уже содержит "chr" в самом имени файла
    независимо от chrom_prefix, это разные вещи). chrom_prefix влияет
    только на содержимое поля CHROM внутри самих строк VCF.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chrom_names = list(chroms) if chroms else list(UPLOAD_CHROMS)
    by_chrom: dict[str, list[str]] = {c: [] for c in chrom_names}
    header_lines: list[str] = []
    opener = gzip.open if str(merged_vcf).endswith(".gz") else open

    # errors="replace" — merged_vcf прошёл через несколько вызовов
    # bcftools (merge/norm/view), которые дописывают в заголовок VCF
    # команду вызова с полным путём к файлу; на Windows при не-ASCII
    # символах в пути (например, кириллица в имени запуска, см.
    # докстринг модуля) эта строка может оказаться не в UTF-8. Строки с
    # координатами/генотипами всегда чистый ASCII и не затрагиваются.
    with opener(merged_vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            parts = line.split("\t")
            if len(parts) >= 1:
                chrom = _normalise_chrom(parts[0])
                if chrom in by_chrom:
                    if chrom_prefix:
                        line = _reprefix_chrom_field(line, chrom_prefix)
                    by_chrom[chrom].append(line)

    # Плоидность мужского nonPAR X: сервер отвергает задание целиком,
    # если у образца в nonPAR встречаются и гаплоидные, и диплоидные
    # записи — а "./." от bcftools merge считается диплоидной. См.
    # normalise_x_ploidy().
    if "X" in by_chrom and by_chrom["X"]:
        by_chrom["X"], fixed = normalise_x_ploidy(by_chrom["X"], genome_build)
        if fixed:
            logger.info(
                "chrX: плоидность приведена к единой в пределах nonPAR — "
                "исправлено %d полей образцов", fixed,
            )

    outputs: list[Path] = []
    for chrom_str in chrom_names:
        chrom = chrom_str
        out_path = output_dir / f"chr{chrom}.vcf.gz"

        def _lines_for_chrom(h=header_lines, d=by_chrom[chrom_str]):
            yield from h
            yield from d

        _write_bgzf(_lines_for_chrom(), out_path, bgzip_path=bgzip_path)
        outputs.append(out_path)
        logger.info("Хромосома %s: %d позиций (chrom_prefix=%r)",
                    chrom, len(by_chrom[chrom_str]), chrom_prefix)

    return outputs


# ---------------------------------------------------------------------------
# QC импутированных данных
# ---------------------------------------------------------------------------
def read_rsq(info_gz: Path) -> dict[tuple[str, int], float]:
    """Читает Rsq из .info.gz."""
    result: dict[tuple[str, int], float] = {}
    # errors="replace" для единообразия с остальными местами модуля,
    # читающими файлы, потенциально прошедшие через внешние инструменты
    # (см. докстринг модуля) — сами Rsq-значения всегда чистый ASCII.
    with gzip.open(info_gz, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n\r").split("\t")
        lower = [x.lower() for x in header]
        pos_idx = next((i for i, x in enumerate(lower) if x in {"position", "pos"}), None)
        chr_idx = next((i for i, x in enumerate(lower) if x in {"chromosome", "chrom", "chr"}), None)
        rsq_idx = next((i for i, x in enumerate(lower) if x == "rsq"), None)
        if pos_idx is None or chr_idx is None or rsq_idx is None:
            raise PureCoreError(f"Не найдены CHROM/POS/Rsq в {info_gz}. Заголовок: {header}")
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            try:
                key = (_normalise_chrom(fields[chr_idx]), int(fields[pos_idx]))
                result[key] = float(fields[rsq_idx])
            except (ValueError, IndexError):
                continue
    return result


def qc_imputed_vcf(
    imputed_vcf: Path,
    info_gz: Path,
    output_vcf: Path,
    rsq_threshold: float = 0.30,
) -> QCResult:
    """QC: если Rsq < threshold, GT заменяется на './.'."""
    rsq = read_rsq(info_gz)
    total = retained = rejected = 0
    missing_in_info = 0

    imputed_vcf = Path(imputed_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(imputed_vcf).endswith(".gz") else open
    out_opener = gzip.open if str(output_vcf).endswith(".gz") else open

    # errors="replace" на чтении — imputed_vcf мог пройти через bcftools
    # локально (см. докстринг модуля про кириллицу/не-ASCII в путях).
    with opener(imputed_vcf, "rt", encoding="utf-8", errors="replace") as fin, \
         out_opener(output_vcf, "wt", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 10:
                continue
            key = (_normalise_chrom(fields[0]), int(fields[1]))
            total += 1
            value = rsq.get(key)
            if value is None:
                missing_in_info += 1
                sample = fields[9].split(":")
                sample[0] = "./."
                fields[9] = ":".join(sample)
                rejected += 1
            elif value < rsq_threshold:
                sample = fields[9].split(":")
                sample[0] = "./."
                fields[9] = ":".join(sample)
                rejected += 1
            else:
                retained += 1
            fout.write("\t".join(fields) + "\n")

    if missing_in_info:
        logger.warning("QC: %d позиций из VCF не найдены в info.gz и обнулены", missing_in_info)
    logger.info("QC: total=%d, retained=%d, rejected=%d", total, retained, rejected)
    return QCResult(total, retained, rejected, rsq_threshold)