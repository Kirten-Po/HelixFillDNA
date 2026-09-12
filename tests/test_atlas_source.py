"""
tests/test_atlas_source.py

Поддержка файлов компании «Атлас» (промт "вот файл компании Атласа").

Файл Атласа оформлен ровно как сырые данные 23andMe — одна '#'-строка с
названиями колонок и четыре колонки через табуляцию, — а координаты в нём
GRCh38. Из этого два требования, которые тут и проверяются:

  1. Автодетект источника обязан отличать его от MyHeritage. До этого
     промта правило "табов больше, чем запятых" опознавало файл Атласа как
     MyHeritage, и GRCh38-координаты уходили на сверку с GRCh37-референсом:
     both_non_ref взлетал, парсер падал на StrandQualityError. Отличить по
     шапке невозможно (она совпадает с 23andMe), поэтому решает проба
     координат по трафарету.
  2. Этап 0 обязан перенести координаты в GRCh37 и НЕ трогать генотипы —
     включая гаплоидные вызовы одной буквой (X/Y/MT) и инделы I/D, которые
     отбраковывает дальше уже адаптер, с подсчётом в QC.

Лифтовер везде синтетический (как в tests/test_liftover.py) — без сети и
без многомегабайтного chain-файла UCSC.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import main
from core import atlas_convert
from core.liftover import ChainLiftover

# chr1: hg38 1..100 -> hg19 501..600; chr2: hg38 1..50 -> hg19 1001..1050.
# Позиции chr1 > 100 и chr2 > 50 не переносятся (за пределами chain).
SYNTHETIC_CHAIN = (
    "chain 1000 chr1 1000 + 0 100 chr1 2000 + 500 600 1\n"
    "100\n"
    "\n"
    "chain 900 chr2 1000 + 0 50 chr2 2000 + 1000 1050 2\n"
    "50\n"
)

ATLAS_HEADER_LINE = "# rsid\tchromosome\tposition\tgenotype\n"


@pytest.fixture
def liftover(tmp_path: Path) -> ChainLiftover:
    chain = tmp_path / "hg38ToHg19.test.chain"
    chain.write_text(SYNTHETIC_CHAIN, encoding="utf-8")
    return ChainLiftover(chain)


def _write_atlas(path: Path, rows: list[tuple[str, str, int, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as f:
        f.write(ATLAS_HEADER_LINE)
        for rsid, chrom, pos, gt in rows:
            f.write(f"{rsid}\t{chrom}\t{pos}\t{gt}\n")
    return path


def _data_rows(path: Path) -> list[list[str]]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        out.append(line.split("\t"))
    return out


# --- Этап 0: перенос координат ---------------------------------------------

def test_convert_lifts_positions_and_keeps_genotypes(tmp_path, liftover):
    src = _write_atlas(tmp_path / "atlas.txt", [
        ("rs1", "1", 10, "AG"),
        ("rs2", "1", 20, "CC"),
        ("rs3", "X", 30, "G"),      # гаплоид одной буквой
        ("rs4", "2", 5, "II"),      # индел
    ])
    dst = tmp_path / "out.txt"
    stats = atlas_convert.convert_atlas_to_grch37(src, dst, liftover)

    rows = _data_rows(dst)
    assert stats.rows_in == 4
    # chrX в синтетическом chain нет — эта строка не переносится.
    assert stats.lift_failed == 1
    assert [r[:4] for r in rows] == [
        ["rs1", "1", "510", "AG"],
        ["rs2", "1", "520", "CC"],
        ["rs4", "2", "1005", "II"],
    ]
    assert stats.rows_out == 3


def test_convert_writes_lf_and_template_header(tmp_path, liftover):
    src = _write_atlas(tmp_path / "atlas.txt", [("rs1", "1", 10, "AG")])
    template = tmp_path / "template_v3.txt"
    template.write_text(
        "# шапка трафарета\n# rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t510\t--\n",
        encoding="utf-8",
    )
    dst = tmp_path / "out.txt"
    atlas_convert.convert_atlas_to_grch37(src, dst, liftover, template_path=template)

    raw = dst.read_bytes()
    assert b"\r\n" not in raw, "у v3 перевод строки — LF"
    assert raw.startswith("# шапка трафарета\n".encode("utf-8"))


def test_convert_sorts_rows_after_lift(tmp_path, liftover):
    """Лифтовер меняет порядок — файл должен остаться сортированным."""
    src = _write_atlas(tmp_path / "atlas.txt", [
        ("rs2", "2", 5, "AA"),
        ("rs1", "1", 20, "GG"),
    ])
    dst = tmp_path / "out.txt"
    atlas_convert.convert_atlas_to_grch37(src, dst, liftover)
    rows = _data_rows(dst)
    assert [(r[1], int(r[2])) for r in rows] == [("1", 520), ("2", 1005)]


def test_convert_drops_non_rsid_markers(tmp_path, liftover):
    src = _write_atlas(tmp_path / "atlas.txt", [
        ("rs1", "1", 10, "AG"),
        ("NC_000006.12:g.18143597T>G", "1", 20, "TG"),
    ])
    dst = tmp_path / "out.txt"
    stats = atlas_convert.convert_atlas_to_grch37(src, dst, liftover)
    assert stats.non_rsid_dropped == 1
    assert [r[0] for r in _data_rows(dst)] == ["rs1"]


def test_convert_drops_duplicate_positions_after_lift(tmp_path, liftover):
    src = _write_atlas(tmp_path / "atlas.txt", [
        ("rs1", "1", 10, "AG"),
        ("rs2", "1", 10, "GG"),
    ])
    dst = tmp_path / "out.txt"
    stats = atlas_convert.convert_atlas_to_grch37(src, dst, liftover)
    assert stats.duplicate_positions == 1
    assert [r[0] for r in _data_rows(dst)] == ["rs1"]


def test_convert_counts_malformed_rows(tmp_path, liftover):
    src = tmp_path / "atlas.txt"
    with src.open("w", encoding="utf-8", newline="") as f:
        f.write(ATLAS_HEADER_LINE)
        f.write("rs1\t1\t10\tAG\n")
        f.write("сломанная строка\n")
        f.write("rs2\t1\tнеположение\tCC\n")
    stats = atlas_convert.convert_atlas_to_grch37(src, tmp_path / "out.txt", liftover)
    assert stats.malformed_rows == 2
    assert stats.rows_out == 1


def test_convert_raises_when_nothing_lifted(tmp_path, liftover):
    src = _write_atlas(tmp_path / "atlas.txt", [("rs1", "9", 10, "AG")])
    with pytest.raises(atlas_convert.AtlasConvertError):
        atlas_convert.convert_atlas_to_grch37(src, tmp_path / "out.txt", liftover)


# --- оформление / идемпотентность ------------------------------------------

def test_detect_layout_atlas_vs_23andme_with_long_header(tmp_path):
    atlas = _write_atlas(tmp_path / "atlas.txt", [("rs1", "1", 10, "AG")])
    assert atlas_convert.detect_layout(atlas) == "atlas"

    long_header = tmp_path / "23andme.txt"
    with long_header.open("w", encoding="utf-8", newline="") as f:
        for i in range(12):
            f.write(f"# строка шапки {i}\n")
        f.write(ATLAS_HEADER_LINE)
        f.write("rs1\t1\t10\tAG\n")
    assert atlas_convert.detect_layout(long_header) is None


def test_prepare_atlas_file_skips_already_converted(tmp_path, liftover):
    src = _write_atlas(
        tmp_path / ("sample" + atlas_convert.CONVERTED_SUFFIX),
        [("rs1", "1", 10, "AG")],
    )
    stats = atlas_convert.prepare_atlas_file(src, tmp_path, liftover)
    assert stats.skipped is True
    assert Path(stats.out_path) == src


def test_prepare_atlas_file_rejects_foreign_format(tmp_path, liftover):
    foreign = tmp_path / "ftdna.csv"
    foreign.write_text("RSID,CHROMOSOME,POSITION,RESULT\nrs1,1,10,AG\n", encoding="utf-8")
    with pytest.raises(atlas_convert.AtlasConvertError):
        atlas_convert.prepare_atlas_file(foreign, tmp_path, liftover)


# --- автодетект источника --------------------------------------------------

def _fake_project_with_template(tmp_path, monkeypatch, positions: dict[str, int]):
    """Подменяет PROJECT_ROOT на временный проект с трафаретом GRCh37."""
    samples = tmp_path / "samples"
    samples.mkdir(exist_ok=True)
    with (samples / "template_genotek.txt").open("w", encoding="utf-8", newline="") as f:
        f.write(ATLAS_HEADER_LINE)
        for rsid, pos in positions.items():
            f.write(f"{rsid}\t1\t{pos}\t--\n")
    monkeypatch.setattr(main, "PROJECT_ROOT", tmp_path)


def test_detect_source_recognises_atlas_by_build(tmp_path, monkeypatch):
    # Трафарет (GRCh37) и файл (GRCh38) расходятся по позициям на тех же rsID.
    grch37 = {f"rs{i}": 1_000_000 + i for i in range(200)}
    _fake_project_with_template(tmp_path, monkeypatch, grch37)

    atlas = _write_atlas(
        tmp_path / "atlas.txt",
        [(f"rs{i}", "1", 2_000_000 + i, "AG") for i in range(200)],
    )
    assert main.detect_source_from_file(atlas) == ("atlas", 0.95)


def test_detect_source_keeps_grch37_file_as_myheritage(tmp_path, monkeypatch):
    """Тот же формат, но координаты GRCh37 — это не Атлас."""
    grch37 = {f"rs{i}": 1_000_000 + i for i in range(200)}
    _fake_project_with_template(tmp_path, monkeypatch, grch37)

    same_build = _write_atlas(
        tmp_path / "other.txt",
        [(f"rs{i}", "1", 1_000_000 + i, "AG") for i in range(200)],
    )
    source, _ = main.detect_source_from_file(same_build)
    assert source != "atlas"


# --- подключение источника к пайплайну -------------------------------------

def test_atlas_is_wired_into_pipeline():
    assert "atlas" in main.SOURCES
    assert main.SOURCES["atlas"]["parser"] is main.parse_atlas
    assert "atlas" in main._SOURCES_NEEDING_CONVERSION
    assert main._needs_reference("atlas") is True
    assert main._supports_liftover("atlas") is True


def test_prepare_source_file_is_noop_for_other_sources(tmp_path):
    src = tmp_path / "ftdna.csv"
    src.write_text("RSID,CHROMOSOME,POSITION,RESULT\n", encoding="utf-8")
    path, stats = main.prepare_source_file("ftdna", src, tmp_path)
    assert path == src and stats is None
