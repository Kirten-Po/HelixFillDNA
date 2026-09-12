"""
tests/test_donor_sample_count.py

Регрессия на "галочка «использовать всех доступных EUR-доноров» ничего не
делает: программа просто пропускает скачивание и оставляет обычные 20".

Причина была не в галочке, а в том, что актуальность кэша доноров
проверялась ТОЛЬКО по chip_signature.txt — то есть по чипу. Размер
донорской подвыборки в сигнатуру не входит, поэтому:
  * main.check_donor_cache() признавал кэш на 20 образцах актуальным, и
    этап скачивания доноров пропускался целиком;
  * а если скачивание всё-таки запускалось, process_chromosome() видел
    готовые kgp_sub_*.vcf.gz и печатал "chr1 уже готов" — файлы на 20
    образцах оставались на месте.

Проверяем оба рубежа плюс сам подсчёт образцов.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import download_donors
import main
from core import donor_cache

PANEL_LINES = (
    ["HG0000%d\tGBR\tEUR\tfemale" % i for i in range(1, 6)]
    + ["NA1900%d\tYRI\tAFR\tmale" % i for i in range(1, 4)]
)


def _write_donor_vcf(path: Path, samples: list[str], chrom: str = "1") -> None:
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
    )
    row = (
        f"{chrom}\t10000\trs1\tA\tG\t.\tPASS\t.\tGT\t"
        + "\t".join(["0|1"] * len(samples))
        + "\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(header + row)


def _donors_dir(tmp_path: Path, n_samples: int) -> Path:
    d = tmp_path / "donors" / "ftdna" / "hrc"
    d.mkdir(parents=True, exist_ok=True)
    (d / donor_cache.SAMPLES_PANEL_FILENAME).write_text(
        "\n".join(PANEL_LINES) + "\n", encoding="utf-8")
    samples = [f"S{i:03d}" for i in range(n_samples)]
    for chrom in list(range(1, 23)) + ["X"]:
        _write_donor_vcf(d / f"kgp_sub_{chrom}.vcf.gz", samples, chrom=str(chrom))
    (d / "chip_signature.txt").write_text("SIG123", encoding="utf-8")
    return d


# --- подсчёт образцов ------------------------------------------------------

def test_vcf_sample_count_reads_chrom_line(tmp_path):
    vcf = tmp_path / "kgp_sub_1.vcf.gz"
    _write_donor_vcf(vcf, ["A", "B", "C"])
    assert donor_cache.vcf_sample_count(vcf) == 3


def test_vcf_sample_count_none_for_missing_or_broken(tmp_path):
    assert donor_cache.vcf_sample_count(tmp_path / "нет.vcf.gz") is None
    broken = tmp_path / "broken.vcf.gz"
    broken.write_bytes(b"not a gzip at all")
    assert donor_cache.vcf_sample_count(broken) is None


def test_available_eur_count_counts_only_eur(tmp_path):
    d = tmp_path / "donors"
    d.mkdir()
    (d / donor_cache.SAMPLES_PANEL_FILENAME).write_text(
        "\n".join(PANEL_LINES) + "\n", encoding="utf-8")
    assert donor_cache.available_eur_count(d) == 5


def test_verdict_unchecked_by_default(tmp_path):
    d = _donors_dir(tmp_path, 2)
    verdict = donor_cache.eur_count_verdict(d, d / "kgp_sub_1.vcf.gz")
    assert verdict.ok is True


def test_verdict_all_available_rejects_small_cache(tmp_path):
    d = _donors_dir(tmp_path, 2)
    verdict = donor_cache.eur_count_verdict(d, d / "kgp_sub_1.vcf.gz", None)
    assert verdict == (False, 2, 5)


def test_verdict_all_available_accepts_full_cache(tmp_path):
    d = _donors_dir(tmp_path, 5)
    assert donor_cache.eur_count_verdict(d, d / "kgp_sub_1.vcf.gz", None).ok


def test_verdict_explicit_count_requires_exact_match(tmp_path):
    d = _donors_dir(tmp_path, 20)
    assert donor_cache.eur_count_verdict(d, d / "kgp_sub_1.vcf.gz", 20).ok
    assert not donor_cache.eur_count_verdict(d, d / "kgp_sub_1.vcf.gz", 30).ok


# --- рубеж 1: проверка кэша перед запуском ---------------------------------

def test_check_donor_cache_accepts_cache_when_count_not_checked(tmp_path):
    _donors_dir(tmp_path, 2)
    donors = main.check_donor_cache("SIG123", "ftdna", tmp_path / "donors", panel="hrc")
    assert len(donors) == 23


def test_check_donor_cache_rejects_small_subsample_when_all_requested(tmp_path):
    _donors_dir(tmp_path, 2)
    with pytest.raises(RuntimeError) as exc:
        main.check_donor_cache(
            "SIG123", "ftdna", tmp_path / "donors", panel="hrc",
            eur_sample_count=None,
        )
    text = str(exc.value)
    assert "2 донорских образцах" in text
    assert "все доступные (5)" in text


def test_check_donor_cache_accepts_matching_explicit_count(tmp_path):
    _donors_dir(tmp_path, 20)
    donors = main.check_donor_cache(
        "SIG123", "ftdna", tmp_path / "donors", panel="hrc",
        eur_sample_count=20,
    )
    assert len(donors) == 23


# --- рубеж 2: инвалидация кэша перед перекачкой ----------------------------

def test_invalidate_by_sample_count_purges_stale_cache(tmp_path):
    d = _donors_dir(tmp_path, 2)
    assert download_donors._invalidate_donor_cache_by_sample_count(d, None) is True
    assert list(d.glob("kgp_sub_*.vcf.gz")) == []
    assert not (d / "chip_signature.txt").exists()
    # Файл панели образцов не привязан к числу доноров — его не трогаем.
    assert (d / donor_cache.SAMPLES_PANEL_FILENAME).exists()


def test_invalidate_by_sample_count_keeps_matching_cache(tmp_path):
    d = _donors_dir(tmp_path, 5)
    assert download_donors._invalidate_donor_cache_by_sample_count(d, None) is False
    assert len(list(d.glob("kgp_sub_*.vcf.gz"))) == 23


def test_invalidate_by_sample_count_does_nothing_when_unchecked(tmp_path):
    d = _donors_dir(tmp_path, 2)
    assert download_donors._invalidate_donor_cache_by_sample_count(d) is False
    assert len(list(d.glob("kgp_sub_*.vcf.gz"))) == 23


def test_invalidate_by_sample_count_on_empty_dir(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    assert download_donors._invalidate_donor_cache_by_sample_count(d, None) is False
