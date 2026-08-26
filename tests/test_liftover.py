"""
tests/test_liftover.py
Юнит-тесты на core.liftover.ChainLiftover — на синтетическом (не UCSC)
chain-файле, без сети и без реального многомегабайтного chain.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from core.liftover import ChainLiftover, LiftoverError  # noqa: E402


SYNTHETIC_CHAIN = (
    "chain 1000 chr1 500 + 0 80 chr1 2000 + 1000 1075 1\n"
    "50\t10\t5\n"
    "20\n"
    "\n"
    "chain 500 chr2 300 + 0 30 chr2 200 - 100 130 2\n"
    "30\n"
)


@pytest.fixture
def chain_path(tmp_path: Path) -> Path:
    p = tmp_path / "test.over.chain"
    p.write_text(SYNTHETIC_CHAIN, encoding="utf-8")
    return p


def test_forward_strand_first_block(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("1", 1) == ("1", 1001)
    assert lo.lift("1", 50) == ("1", 1050)


def test_forward_strand_gap_between_blocks(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("1", 51) is None
    assert lo.lift("1", 60) is None


def test_forward_strand_second_block(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("1", 61) == ("1", 1056)
    assert lo.lift("1", 80) == ("1", 1075)


def test_forward_strand_beyond_chain_is_none(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("1", 81) is None


def test_reverse_strand_math(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("2", 1) == ("2", 100)
    assert lo.lift("2", 30) == ("2", 71)


def test_unknown_chromosome_returns_none(chain_path):
    lo = ChainLiftover(chain_path)
    assert lo.lift("99", 1) is None


def test_stats_are_tracked(chain_path):
    lo = ChainLiftover(chain_path)
    lo.lift("1", 1)
    lo.lift("1", 51)
    lo.lift("99", 1)
    assert lo.stats.total_calls == 3
    assert lo.stats.lifted == 1
    assert lo.stats.in_gap == 1
    assert lo.stats.no_chain_for_chrom == 1
    assert lo.stats.failed == 2


def test_missing_file_raises(tmp_path):
    with pytest.raises(LiftoverError):
        ChainLiftover(tmp_path / "does_not_exist.chain")


def test_bad_tstrand_raises(tmp_path):
    p = tmp_path / "bad.chain"
    p.write_text(
        "chain 100 chr1 500 - 0 10 chr1 500 + 0 10 1\n10\n",
        encoding="utf-8",
    )
    with pytest.raises(LiftoverError):
        ChainLiftover(p)


def test_malformed_header_raises(tmp_path):
    p = tmp_path / "bad2.chain"
    p.write_text("chain 100 chr1 500 + 0 10\n10\n", encoding="utf-8")
    with pytest.raises(LiftoverError):
        ChainLiftover(p)


def test_malformed_block_line_raises(tmp_path):
    p = tmp_path / "bad3.chain"
    p.write_text(
        "chain 100 chr1 500 + 0 10 chr1 500 + 0 10 1\n10\tX\tY\n",
        encoding="utf-8",
    )
    with pytest.raises(LiftoverError):
        ChainLiftover(p)


def test_gzipped_chain_file_reads(tmp_path):
    import gzip
    p = tmp_path / "test.over.chain.gz"
    with gzip.open(p, "wt", encoding="utf-8") as f:
        f.write(SYNTHETIC_CHAIN)
    lo = ChainLiftover(p)
    assert lo.lift("1", 1) == ("1", 1001)
