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
    """
    validate_variants(result)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    variants = sorted(
        result.variants,
        key=lambda v: (_chrom_sort_key(v.chrom), int(v.pos)),
    )

    if not compress:
        with output_vcf.open("w", encoding="utf-8", newline="\n") as f:
            _write_vcf_header(f, sample_name)
            for v in variants:
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix)
        final_path = output_vcf
    elif bgzip_path:
        # Настоящий BGZF через внешний бинарник
        tmp_vcf = output_vcf.with_suffix(".vcf")
        with tmp_vcf.open("w", encoding="utf-8", newline="\n") as f:
            _write_vcf_header(f, sample_name)
            for v in variants:
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix)
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
                _write_vcf_line(f, v, chrom_prefix=chrom_prefix)
        final_path = output_vcf

    logger.info("VCF собран: %d вариантов в %s (chrom_prefix=%r)",
                len(variants), final_path, chrom_prefix)
    return final_path


def _write_vcf_header(f, sample_name: str) -> None:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=PurePythonCore\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")


def _write_vcf_line(f, v: ParsedVariant, chrom_prefix: str = "") -> None:
    """
    chrom_prefix подставляется ЗДЕСЬ, в момент записи, и только здесь —
    v.chrom/_normalise_chrom() остаются каноническими (без префикса) везде
    в остальном модуле (сортировка, dict-ключи, валидация).
    """
    chrom = f"{chrom_prefix}{_normalise_chrom(v.chrom)}"
    f.write(
        f"{chrom}\t{int(v.pos)}\t{v.rsid}\t{v.ref}\t{v.alt}"
        f"\t.\tPASS\t.\tGT\t{_vcf_gt(v.gt)}\n"
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


def split_autosomes(
    merged_vcf: Path,
    output_dir: Path,
    bgzip_path: Optional[str] = None,
    chrom_prefix: str = "",
) -> list[Path]:
    """
    Делит merged VCF на chr1..chr22.
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

    by_chrom: dict[str, list[str]] = {str(i): [] for i in range(1, 23)}
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

    outputs: list[Path] = []
    for chrom in range(1, 23):
        chrom_str = str(chrom)
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