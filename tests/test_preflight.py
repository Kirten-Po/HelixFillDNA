"""
tests/test_preflight.py

Отчёт "состояние исходника" — единственное место, где приложение говорит
про файл до того, как потрачен час на доноров и сутки на очередь MIS.
Ошибка здесь тихая: отчёт будет показан, просто с неверными числами, и
человек примет по ним неверное решение. Поэтому проверяем на файлах с
ЗАРАНЕЕ ИЗВЕСТНЫМ ответом, а не на "смотри, не падает".
"""
from __future__ import annotations

import gzip

import pytest

from core import preflight


def _write_ftdna(path, rows):
    lines = ["RSID,CHROMOSOME,POSITION,RESULT"]
    lines += [f"{r},{c},{p},{g}" for r, c, p, g in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_call_rate_counts_no_calls(tmp_path):
    rows = [(f"rs{i}", "1", 1000 + i, "AG") for i in range(8)]
    rows += [(f"rs9{i}", "1", 5000 + i, "--") for i in range(2)]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    assert r.markers == 10
    assert r.called == 8
    assert r.call_rate == pytest.approx(80.0)


def test_heterozygosity_and_palindromes(tmp_path):
    rows = [
        ("rs1", "1", 100, "AG"),   # гетерозигота, не палиндром
        ("rs2", "1", 200, "AT"),   # гетерозигота, палиндром A/T
        ("rs3", "1", 300, "CG"),   # гетерозигота, палиндром C/G
        ("rs4", "1", 400, "AA"),   # гомозигота
        ("rs5", "1", 500, "GG"),
    ]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    assert r.autosomal_called == 5
    assert r.autosomal_het == 3
    assert r.autosomal_het_pct == pytest.approx(60.0)
    # Палиндромность считается ОТ ГЕТЕРОЗИГОТ: по "AA" о ней судить нельзя.
    assert r.palindromic_het == 2
    assert r.palindromic_het_pct == pytest.approx(200.0 / 3)


def test_par_positions_excluded_from_sex_call(tmp_path):
    """
    Гетерозиготность в псевдоаутосомных регионах нормальна и у мужчин —
    смешать её в оценку пола значит объявить каждого мужчину женщиной.
    """
    rows = [("rsp", "X", 1_000_000, "AG")] * 1          # PAR1 (GRCh37)
    rows = [(f"rsp{i}", "X", 1_000_000 + i, "AG") for i in range(300)]
    rows += [(f"rsn{i}", "X", 50_000_000 + i, "AA") for i in range(300)]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    assert r.x_nonpar_calls == 300, "в оценку пола должен идти только nonPAR"
    assert r.x_het_pct == pytest.approx(0.0)
    assert r.sex_by_x == "мужской"


def test_female_by_x(tmp_path):
    rows = [(f"rs{i}", "X", 50_000_000 + i,
             "AG" if i % 2 else "AA") for i in range(400)]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    assert r.sex_by_x == "женский"


def test_chrom_code_map_is_vendor_specific(tmp_path):
    """
    У AncestryDNA код "25" — это псевдоаутосомный регион X, у MyHeritage
    тот же "25" означает MT. Перепутать карты — тихо перебросить десятки
    позиций X в митохондриальный контиг.
    """
    assert preflight.normalise_chrom("25", "ancestry") == "X"
    assert preflight.normalise_chrom("25", "myheritage") == "MT"
    assert preflight.normalise_chrom("chr7", "vcf") == "7"
    assert preflight.normalise_chrom("XY", "ftdna") == "X"


def test_duplicates_and_sort_order(tmp_path):
    rows = [
        ("rs1", "1", 100, "AA"),
        ("rs1", "1", 200, "AA"),   # дубль rsID
        ("rs2", "1", 200, "AA"),   # дубль позиции
        ("rs3", "1", 150, "AA"),   # нарушение порядка
    ]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    assert r.duplicate_rsids == 1
    assert r.duplicate_positions == 1
    assert r.unsorted_positions == 1


def test_missing_y_and_mt_are_reported(tmp_path):
    rows = [(f"rs{i}", "1", 100 + i, "AA") for i in range(50)]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    titles = " ".join(f.title for f in r.findings)
    assert "Y-хромосом" in titles
    assert "Митохондриальной" in titles


def test_build_detected_against_template(tmp_path):
    """
    Сборка определяется сверкой с трафаретом — настоящим экспортом 23andMe
    в GRCh37, который и так лежит в комплекте. Никаких зашитых по памяти
    координат реперных rsID.
    """
    template = tmp_path / "template_v3.txt"
    template.write_text(
        "# header\n" + "\n".join(f"rs{i}\t1\t{1000 + i}\tAA" for i in range(200)),
        encoding="utf-8",
    )
    same = [(f"rs{i}", "1", 1000 + i, "AA") for i in range(200)]
    other = [(f"rs{i}", "1", 9_000_000 + i, "AA") for i in range(200)]

    r_ok = preflight.analyse_file(
        _write_ftdna(tmp_path / "a.csv", same), "ftdna",
        template_path=template,
    )
    assert r_ok.build == "GRCh37"
    assert r_ok.build_match_pct == pytest.approx(100.0)

    r_bad = preflight.analyse_file(
        _write_ftdna(tmp_path / "b.csv", other), "ftdna",
        template_path=template,
    )
    assert r_bad.build.startswith("не GRCh37")
    assert any(f.level == "bad" and "Сборка" in f.title for f in r_bad.findings)


def test_vcf_source_reads_genotypes(tmp_path):
    vcf = tmp_path / "s.vcf.gz"
    body = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
        "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1",
        "1\t200\trs2\tA\tG\t.\t.\t.\tGT\t./.",
        "1\t300\trs3\tAT\tA\t.\t.\t.\tGT\t0/1",   # индель — не SNP, пропускается
    ]
    with gzip.open(vcf, "wt", encoding="utf-8") as f:
        f.write("\n".join(body) + "\n")
    r = preflight.analyse_file(vcf, "vcf")
    assert r.markers == 2, "инделя в статистике чипа быть не должно"
    assert r.called == 1
    assert r.autosomal_het == 1


def test_unreadable_file_raises_clear_error(tmp_path):
    junk = tmp_path / "j.csv"
    junk.write_text("это не экспорт чипа\nи не станет им\n", encoding="utf-8")
    with pytest.raises(preflight.PreflightError):
        preflight.analyse_file(junk, "ftdna")


def test_metrics_row_is_flat(tmp_path):
    rows = [(f"rs{i}", "1", 100 + i, "AG") for i in range(20)]
    r = preflight.analyse_file(_write_ftdna(tmp_path / "f.csv", rows), "ftdna")
    row = preflight.metrics_row(r)
    assert all(not isinstance(v, (dict, list)) for v in row.values()), (
        "строка для CSV должна быть плоской"
    )
    assert row["markers"] == 20
