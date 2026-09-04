"""
tests/test_rsq_info_formats.py

Регрессионные тесты на разбор chr*.info.gz — файла, из которого берётся
качество импутации (Rsq / R2).

Michigan Imputation Server отдаёт его в двух форматах: классическом TSV
с колонкой Rsq и sites-only VCF, где качество лежит в поле INFO как
"R2=". Прежний разбор умел только первый. Встретив второй, он читал
первой строкой "##fileformat=VCFv4.2", не находил ни одной нужной
колонки и МОЛЧА возвращал пустую карту — после чего вызывающий код
подставлял Rsq=1.0 каждой позиции, и порог качества переставал отсекать
что-либо вообще, без единого сообщения в лог.

Найдено разбором реального результата: в итоговый файл прошло 10-12 %
импутированных позиций с R2 < 0.3 при выставленном пороге 0.30.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from core.pure_python_core import (  # noqa: E402
    PureCoreError, read_rsq, read_rsq_map,
)


def _write_gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as f:
        f.write(text)
    return path


# --- формат 1: классический TSV -------------------------------------------
TSV_INFO = """SNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped
22:16050435\tT\tC\t0.0001\t0.0001\t0.99\t0.812\tImputed
22:16050783\tA\tG\t0.0001\t0.0001\t0.99\t0.145\tImputed
"""


def test_reads_classic_tsv(tmp_path):
    info = _write_gz(tmp_path / "chr22.info.gz",
                     TSV_INFO.replace("SNP\t", "chromosome\tposition\t")
                             .replace("22:16050435\t", "22\t16050435\t")
                             .replace("22:16050783\t", "22\t16050783\t"))
    m = read_rsq_map(info)
    assert m == {("22", 16050435): 0.812, ("22", 16050783): 0.145}


# --- формат 2: sites-only VCF ---------------------------------------------
VCF_INFO = """##fileformat=VCFv4.2
##filedate=20260828
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
22\t16050435\t22:16050435\tT\tC\t.\t.\tAF=0.0001;MAF=0.0001;AVG_CS=0.99;R2=0.812;IMPUTED
22\t16050783\t22:16050783\tA\tG\t.\t.\tAF=0.0001;MAF=0.0001;AVG_CS=0.99;R2=0.145;IMPUTED
"""


def test_reads_sites_only_vcf(tmp_path):
    """Тот самый формат, на котором прежний разбор молча сдавался."""
    info = _write_gz(tmp_path / "chr22.info.gz", VCF_INFO)
    assert read_rsq_map(info) == {
        ("22", 16050435): 0.812,
        ("22", 16050783): 0.145,
    }


def test_chr_prefix_is_normalised(tmp_path):
    """Для TopMed/GRCh38 CHROM записан как "chr22" — ключи всё равно
    должны быть каноническими, иначе поиск по ним промахнётся."""
    info = _write_gz(tmp_path / "chrX.info.gz",
                     VCF_INFO.replace("\n22\t", "\nchrX\t").replace("22:", "X:"))
    assert set(read_rsq_map(info)) == {("X", 16050435), ("X", 16050783)}


def test_r2_is_not_confused_with_other_info_fields(tmp_path):
    """В INFO есть поля, содержащие "R2" как подстроку (например
    AVG_CS_R2 или ER2) — regexp обязан цепляться только за само поле."""
    info = _write_gz(tmp_path / "chr1.info.gz",
                     "##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     "1\t100\t.\tA\tG\t.\t.\tER2=0.999;R2=0.42;IMPUTED\n")
    assert read_rsq_map(info) == {("1", 100): 0.42}


def test_typed_only_sites_without_r2_are_skipped(tmp_path):
    """У реальных измерений чипа качества импутации нет — их отсутствие
    в карте означает "фильтровать нечего", и это верно."""
    info = _write_gz(tmp_path / "chr1.info.gz",
                     "##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     "1\t100\t.\tA\tG\t.\t.\tAF=0.5;TYPED_ONLY\n"
                     "1\t200\t.\tC\tT\t.\t.\tAF=0.5;R2=0.9;IMPUTED\n")
    assert read_rsq_map(info) == {("1", 200): 0.9}


# --- поведение при нераспознанном файле -----------------------------------
def test_unrecognised_format_warns_loudly(tmp_path, caplog):
    """Молчаливое "качество у всех идеальное" — худший вариант; пустая
    карта обязана быть заметна в логе."""
    info = _write_gz(tmp_path / "chr1.info.gz", "какая-то ерунда\nи ещё строка\n")
    with caplog.at_level("WARNING"):
        assert read_rsq_map(info) == {}
    assert any("НЕ РАБОТАЕТ" in r.message for r in caplog.records), caplog.text


def test_missing_file_is_not_an_error(tmp_path):
    """Хромосомы может не быть в результатах (Y, MT, иногда X) — это
    штатная ситуация, а не сбой."""
    assert read_rsq_map(tmp_path / "chr99.info.gz") == {}


def test_strict_wrapper_raises_on_unrecognised(tmp_path):
    """qc_imputed_vcf() при пустой карте забраковал бы КАЖДУЮ позицию,
    поэтому для него нераспознанный файл — явная ошибка."""
    info = _write_gz(tmp_path / "chr1.info.gz", "мусор\n")
    with pytest.raises(PureCoreError, match="формат не распознан"):
        read_rsq(info)


def test_strict_wrapper_passes_through_valid_file(tmp_path):
    info = _write_gz(tmp_path / "chr22.info.gz", VCF_INFO)
    assert read_rsq(info)[("22", 16050435)] == 0.812
