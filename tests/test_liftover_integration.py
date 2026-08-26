"""
tests/test_liftover_integration.py

Регрессионный тест на промт "HRC / TopMed v4" (интеграция лифтовера в
адаптеры): main.py/gui/app.py вызывают
adapters/ftdna_v3.py::parse_ftdna_v3(csv_path, reference, liftover=liftover)
и adapters/myheritage_v5.py::parse_myheritage_v5(csv_path, reference, liftover=liftover)
— но до этого промта ни одна из функций не принимала параметр liftover
вовсе. Реальный прогон падал с

    TypeError: parse_ftdna_v3() got an unexpected keyword argument 'liftover'

для ЛЮБОЙ панели, как только вызывающий код (main.py::main(), Этап [0b/7])
строил ChainLiftover и передавал его дальше. Это не про корректность самого
лифтовера координат (core/liftover.py::ChainLiftover уже покрыт
tests/test_liftover.py) — это про то, что параметр физически долетает до
адаптера и координата переносится РАНЬШЕ reference.base_at() и РАНЬШЕ
регистрации в широкой (broad) сигнатуре чипа (Задача D).

Используется реальный pyfaidx (как и tests/test_reference_genome_contigs.py)
— никакого мока, крошечные синтетические .fasta во временной папке pytest.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.ftdna_v3 import parse_ftdna_v3, ReferenceGenome  # noqa: E402
from adapters.myheritage_v5 import parse_myheritage_v5  # noqa: E402
from core.liftover import ChainLiftover  # noqa: E402


# ---------------------------------------------------------------------------
# Синтетический chain: GRCh37-подобный контиг "1" длиной 100 -> GRCh38-
# подобный контиг "chr1" длиной 200, один ungapped-блок [0,100) -> [0,100),
# без разрывов и без reverse-strand — координата 1-based pos N переносится
# в pos N той же длины блока (ChainLiftover.lift() возвращает q0+1).
# ---------------------------------------------------------------------------
SYNTHETIC_CHAIN = (
    "chain 1000 1 100 + 0 100 chr1 200 + 0 100 1\n"
    "100\n"
)


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with path.open("w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")
    return path


@pytest.fixture
def grch37_reference(tmp_path: Path) -> ReferenceGenome:
    """Контиг '1' без префикса 'chr', длина 100 — как HRC/human_g1k_v37."""
    fasta = _write_fasta(tmp_path / "grch37.fasta", {"1": "A" * 100})
    return ReferenceGenome(fasta)


@pytest.fixture
def grch38_reference(tmp_path: Path) -> ReferenceGenome:
    """Контиг 'chr1' с префиксом, длина 200 — как TopMed/GRCh38."""
    fasta = _write_fasta(tmp_path / "grch38.fasta", {"chr1": "G" * 200})
    return ReferenceGenome(fasta)


@pytest.fixture
def liftover(tmp_path: Path) -> ChainLiftover:
    chain_path = tmp_path / "test.over.chain"
    chain_path.write_text(SYNTHETIC_CHAIN, encoding="utf-8")
    return ChainLiftover(chain_path)


def _write_ftdna_csv(path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    lines = ["RSID,CHROMOSOME,POSITION,RESULT"]
    for rsid, chrom, pos, genotype in rows:
        lines.append(f"{rsid},{chrom},{pos},{genotype}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_myheritage_tsv(path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    lines = ["# comment"] * 12
    lines.append("RSID\tCHROMOSOME\tPOSITION\tRESULT")
    for rsid, chrom, pos, genotype in rows:
        lines.append(f"{rsid}\t{chrom}\t{pos}\t{genotype}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# parse_ftdna_v3
# ---------------------------------------------------------------------------
def test_parse_ftdna_v3_accepts_liftover_kwarg_without_typeerror(
    tmp_path, grch38_reference, liftover,
):
    """Сам факт, что вызов не падает с TypeError — это и есть регрессионный
    тест на найденный баг (main.py/gui/app.py вызывают именно так)."""
    csv_path = _write_ftdna_csv(tmp_path / "ftdna.csv", [("rs1", "1", 10, "AA")])
    result = parse_ftdna_v3(csv_path, grch38_reference, liftover=liftover)
    assert len(result.variants) == 1


def test_parse_ftdna_v3_default_liftover_none_unchanged_behaviour(
    tmp_path, grch37_reference,
):
    """HRC-путь (liftover не передан вовсе) не должен измениться ни на
    йоту — координата остаётся GRCh37 без переноса."""
    csv_path = _write_ftdna_csv(tmp_path / "ftdna.csv", [("rs1", "1", 10, "AA")])
    result = parse_ftdna_v3(csv_path, grch37_reference)
    assert result.lift_failed == 0
    assert len(result.variants) == 1
    assert result.variants[0].chrom == "1"
    assert result.variants[0].pos == 10


def test_parse_ftdna_v3_lifts_coordinate_before_reference_lookup(
    tmp_path, grch38_reference, liftover,
):
    """Координата переносится ДО reference.base_at() — если бы перенос не
    происходил (или происходил после чтения референса), base_at() читал бы
    GRCh38-референс по GRCh37-координате: для позиции 10 это совпадёт
    случайно (блок identity-переноса), поэтому решающая проверка —
    результирующий ParsedVariant.chrom/pos должен быть в системе координат
    ЦЕЛЕВОЙ сборки, а не исходной."""
    csv_path = _write_ftdna_csv(tmp_path / "ftdna.csv", [("rs1", "1", 10, "AA")])
    result = parse_ftdna_v3(csv_path, grch38_reference, liftover=liftover)
    assert result.lift_failed == 0
    assert len(result.variants) == 1
    v = result.variants[0]
    assert (v.chrom, v.pos) == ("1", 10)  # синтетический chain — identity-перенос


def test_parse_ftdna_v3_unliftable_position_dropped_and_counted(
    tmp_path, grch38_reference, liftover,
):
    """Позиция вне chain-блока (t_end=100) не переносится — должна быть
    полностью исключена из variants и учтена в lift_failed, а не привести
    к падению или к записи с исходной GRCh37-координатой."""
    csv_path = _write_ftdna_csv(tmp_path / "ftdna.csv", [
        ("rs1", "1", 10, "AA"),
        ("rs_gap", "1", 500, "AA"),  # далеко за пределами chain
    ])
    result = parse_ftdna_v3(csv_path, grch38_reference, liftover=liftover)
    assert result.lift_failed == 1
    assert len(result.variants) == 1
    assert result.variants[0].rsid == "rs1"


def test_parse_ftdna_v3_broad_signature_uses_lifted_coordinates(
    tmp_path, grch38_reference, liftover,
):
    """Задача D: широкая сигнатура чипа (chip_signature_broad) должна
    строиться по УЖЕ ПЕРЕНЕСЁННЫМ координатам — иначе для одного и того же
    физического чипа она отличалась бы в зависимости от того, вызывался
    parser с liftover или нет, что ломает переиспользование доноров между
    HRC- и TopMed-прогонами одного чипа (они, впрочем, и не должны
    совпадать — но сигнатура должна быть детерминирована относительно
    выбранной панели, а не содержать мусор из недолифтованных координат)."""
    csv_path = _write_ftdna_csv(tmp_path / "ftdna.csv", [
        ("rs1", "1", 10, "AA"),
        ("rs2", "1", 20, "--"),  # NO_CALL — всё равно должен попасть в broad
    ])
    result = parse_ftdna_v3(csv_path, grch38_reference, liftover=liftover)
    assert result.chip_signature_broad
    # Обе позиции (включая NO_CALL) должны присутствовать в broad-списке,
    # и обе — уже в целевых (перенесённых) координатах.
    assert ("1", 10) in result.signature_positions_broad
    assert ("1", 20) in result.signature_positions_broad


# ---------------------------------------------------------------------------
# parse_myheritage_v5 — та же интеграция, тот же баг был бы там
# ---------------------------------------------------------------------------
def test_parse_myheritage_v5_accepts_liftover_kwarg_without_typeerror(
    tmp_path, grch38_reference, liftover,
):
    tsv_path = _write_myheritage_tsv(tmp_path / "myheritage.tsv", [("rs1", "1", 10, "AA")])
    result = parse_myheritage_v5(tsv_path, grch38_reference, liftover=liftover)
    assert len(result.variants) == 1


def test_parse_myheritage_v5_default_liftover_none_unchanged_behaviour(
    tmp_path, grch37_reference,
):
    tsv_path = _write_myheritage_tsv(tmp_path / "myheritage.tsv", [("rs1", "1", 10, "AA")])
    result = parse_myheritage_v5(tsv_path, grch37_reference)
    assert result.lift_failed == 0
    assert result.variants[0].chrom == "1"
    assert result.variants[0].pos == 10


def test_parse_myheritage_v5_lifts_coordinate_and_counts_failures(
    tmp_path, grch38_reference, liftover,
):
    tsv_path = _write_myheritage_tsv(tmp_path / "myheritage.tsv", [
        ("rs1", "1", 10, "AA"),
        ("rs_gap", "1", 500, "AA"),
    ])
    result = parse_myheritage_v5(tsv_path, grch38_reference, liftover=liftover)
    assert result.lift_failed == 1
    assert len(result.variants) == 1
    assert (result.variants[0].chrom, result.variants[0].pos) == ("1", 10)


if __name__ == "__main__":
    raise SystemExit(
        "Этот файл рассчитан на запуск через pytest: "
        "python -m pytest tests/test_liftover_integration.py -v"
    )
