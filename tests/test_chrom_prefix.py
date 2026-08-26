"""
tests/test_chrom_prefix.py

Регрессионный тест на промт "HRC / TopMed" (v4-фикс): main.py вызывает
core/pure_python_core.py::build_vcf()/split_autosomes() с параметром
chrom_prefix, которого раньше не было в сигнатуре этих функций — реальный
прогон падал с TypeError на Этапе 2 для ЛЮБОЙ панели, включая HRC. Этот
тест ловит именно такую регрессию: он не про корректность самого лифтовера
координат, а про то, что "" (HRC) не меняет поведение, а "chr" (TopMed)
реально попадает в CHROM выходного VCF, и только туда — не в места
сравнения/сортировки/индексации (_normalise_chrom()/_chrom_sort_key()
остаются каноническими).
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.base import ParsedVariant, ParseResult  # noqa: E402
from core.pure_python_core import build_vcf, split_autosomes  # noqa: E402


def _sample_result() -> ParseResult:
    result = ParseResult()
    result.variants = [
        ParsedVariant(rsid="rs1", chrom="1", pos=100, ref="A", alt="G", gt="0/1"),
        ParsedVariant(rsid="rs2", chrom="2", pos=200, ref="C", alt=".", gt="0/0"),
        ParsedVariant(rsid="rsX", chrom="X", pos=300, ref="T", alt="C", gt="1/1"),
    ]
    result.total_measured = len(result.variants)
    return result


def _data_lines(vcf_path: Path) -> list[str]:
    with vcf_path.open("r", encoding="utf-8") as f:
        return [line.rstrip("\n") for line in f if not line.startswith("#")]


# ---------------------------------------------------------------------------
# build_vcf()
# ---------------------------------------------------------------------------
def test_build_vcf_default_chrom_prefix_is_empty(tmp_path):
    """Без явного chrom_prefix (или chrom_prefix="") поведение идентично
    тому, что было ДО добавления параметра — HRC не должен измениться
    ни на байт."""
    out = tmp_path / "sample.vcf"
    build_vcf(_sample_result(), out, compress=False)

    lines = _data_lines(out)
    chroms = [line.split("\t")[0] for line in lines]
    assert chroms == ["1", "2", "X"], f"CHROM не должен содержать префикс по умолчанию: {chroms}"


def test_build_vcf_explicit_empty_prefix_matches_default(tmp_path):
    out_default = tmp_path / "default.vcf"
    out_explicit = tmp_path / "explicit_empty.vcf"
    build_vcf(_sample_result(), out_default, compress=False)
    build_vcf(_sample_result(), out_explicit, compress=False, chrom_prefix="")

    assert _data_lines(out_default) == _data_lines(out_explicit)


def test_build_vcf_with_chr_prefix(tmp_path):
    out = tmp_path / "sample.vcf"
    build_vcf(_sample_result(), out, compress=False, chrom_prefix="chr")

    lines = _data_lines(out)
    chroms = [line.split("\t")[0] for line in lines]
    assert chroms == ["chr1", "chr2", "chrX"], (
        f"CHROM должен содержать префикс 'chr' для каждой строки: {chroms}"
    )

    # Остальные поля строки (POS/ID/REF/ALT/GT) не должны меняться от
    # наличия префикса — сравниваем построчно с "беспрефиксной" версией.
    out_no_prefix = tmp_path / "sample_no_prefix.vcf"
    build_vcf(_sample_result(), out_no_prefix, compress=False)
    lines_no_prefix = _data_lines(out_no_prefix)

    for with_prefix, without_prefix in zip(lines, lines_no_prefix):
        parts_with = with_prefix.split("\t")
        parts_without = without_prefix.split("\t")
        assert parts_with[1:] == parts_without[1:], (
            "Только CHROM должен отличаться из-за chrom_prefix, остальные "
            "поля обязаны совпадать"
        )
        assert parts_with[0] == f"chr{parts_without[0]}"


def test_build_vcf_chrom_prefix_does_not_affect_sort_order(tmp_path):
    """Сортировка вариантов (_chrom_sort_key) работает по каноническому
    имени — порядок строк в файле не должен зависеть от chrom_prefix."""
    out_prefixed = tmp_path / "prefixed.vcf"
    out_plain = tmp_path / "plain.vcf"
    build_vcf(_sample_result(), out_prefixed, compress=False, chrom_prefix="chr")
    build_vcf(_sample_result(), out_plain, compress=False)

    rsids_prefixed = [line.split("\t")[2] for line in _data_lines(out_prefixed)]
    rsids_plain = [line.split("\t")[2] for line in _data_lines(out_plain)]
    assert rsids_prefixed == rsids_plain == ["rs1", "rs2", "rsX"]


# ---------------------------------------------------------------------------
# split_autosomes()
# ---------------------------------------------------------------------------
def _write_plain_merged_vcf(path: Path, chrom_lines: dict[str, list[str]]) -> None:
    """
    Пишет НЕсжатый (без .gz) merged_vcf-подобный файл — split_autosomes()
    сама выбирает open()/gzip.open() по расширению имени файла, так что
    для теста достаточно обычного текстового файла с расширением .vcf,
    без реальной BGZF-компрессии на входе.
    """
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
        for chrom, lines in chrom_lines.items():
            for line in lines:
                f.write(line + "\n")


def _read_bgzf_or_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return [l.rstrip("\n") for l in f if not l.startswith("#")]


@pytest.fixture
def _bgzip_writer_available():
    """
    split_autosomes()/_write_bgzf() без bgzip_path откатывается на
    Bio.bgzf (biopython) — если он не установлен в тестовом окружении,
    пропускаем тест записи (сам чтение/расчёт группировки по хромосомам
    в build_vcf-тестах выше это не затрагивает, они не зависят от
    _write_bgzf вовсе).
    """
    pytest.importorskip("Bio.bgzf", reason="biopython не установлен — split_autosomes() недоступна без него/bgzip_path")


def test_split_autosomes_default_no_prefix(tmp_path, _bgzip_writer_available):
    merged = tmp_path / "merged.vcf"  # без .gz — читается как обычный текст
    _write_plain_merged_vcf(merged, {
        "1": ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"],
        "2": ["2\t200\trs2\tC\t.\t.\tPASS\t.\tGT\t0/0"],
    })
    output_dir = tmp_path / "upload"
    outputs = split_autosomes(merged, output_dir)

    chr1_lines = _read_bgzf_or_gzip_lines(output_dir / "chr1.vcf.gz")
    assert chr1_lines == ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"]
    chr2_lines = _read_bgzf_or_gzip_lines(output_dir / "chr2.vcf.gz")
    assert chr2_lines == ["2\t200\trs2\tC\t.\t.\tPASS\t.\tGT\t0/0"]
    assert len(outputs) == 22


def test_split_autosomes_with_chr_prefix(tmp_path, _bgzip_writer_available):
    merged = tmp_path / "merged.vcf"
    _write_plain_merged_vcf(merged, {
        "1": ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"],
    })
    output_dir = tmp_path / "upload"
    split_autosomes(merged, output_dir, chrom_prefix="chr")

    chr1_lines = _read_bgzf_or_gzip_lines(output_dir / "chr1.vcf.gz")
    assert chr1_lines == ["chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"], (
        "CHROM в выходной строке должен получить префикс 'chr', "
        "остальные поля — без изменений"
    )


def test_split_autosomes_input_already_has_chr_prefix_is_normalised(tmp_path, _bgzip_writer_available):
    """
    Входной merged_vcf может уже содержать 'chr1' в CHROM (например, если
    он собирался под TopMed) — split_autosomes() должен привести его к
    каноническому виду через _normalise_chrom() перед перегруппировкой и
    подстановкой chrom_prefix заново, а не задваивать префикс ('chrchr1').
    """
    merged = tmp_path / "merged.vcf"
    _write_plain_merged_vcf(merged, {
        "1": ["chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"],
    })
    output_dir = tmp_path / "upload"
    split_autosomes(merged, output_dir, chrom_prefix="chr")

    chr1_lines = _read_bgzf_or_gzip_lines(output_dir / "chr1.vcf.gz")
    assert chr1_lines == ["chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"], (
        f"Ожидался ровно один префикс 'chr', без задвоения: {chr1_lines}"
    )


if __name__ == "__main__":
    import tempfile

    print("Запускаю тесты build_vcf() без pytest...")
    with tempfile.TemporaryDirectory() as td:
        test_build_vcf_default_chrom_prefix_is_empty(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_build_vcf_explicit_empty_prefix_matches_default(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_build_vcf_with_chr_prefix(Path(td))
    with tempfile.TemporaryDirectory() as td:
        test_build_vcf_chrom_prefix_does_not_affect_sort_order(Path(td))
    print("OK: build_vcf() тесты прошли. Тесты split_autosomes() запускайте через pytest "
          "(требуют biopython или pytest.importorskip корректно их пропустит).")
