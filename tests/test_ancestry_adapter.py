"""
tests/test_ancestry_adapter.py

Тесты адаптера AncestryDNA (adapters/ancestry_v2.py) и его подключения к
пайплайну (main.SOURCES / detect_source_from_file).

Референс здесь подменён заглушкой (_StubReference) — настоящий
ReferenceGenome тянет за собой .fasta на несколько ГБ, а проверяем мы
разбор формата и нормализацию, а не pyfaidx.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

import main as pipeline
from adapters.ancestry_v2 import (
    AncestryFormatError,
    CHROM_NORMALIZE,
    LAYOUT_RAW,
    _find_header,
    _normalize_chrom,
    parse_ancestry_v2,
)


# ---------------------------------------------------------------------------
# Заглушка референса
# ---------------------------------------------------------------------------
class _StubReference:
    """Отдаёт заранее заданное основание для каждой (chrom, pos)."""

    def __init__(self, bases: dict[tuple[str, int], str], default: str = "A"):
        self._bases = bases
        self._default = default

    def base_at(self, chrom: str, pos: int) -> str:
        return self._bases.get((chrom, pos), self._default)


ANCESTRY_HEADER = "\n".join(
    ["#AncestryDNA raw data download"]
    + [f"#comment line {i}" for i in range(17)]
    + ["rsid\tchromosome\tposition\tallele1\tallele2"]
)


def _write(tmp_path: Path, body: str, name: str = "AncestryDNA.txt") -> Path:
    path = tmp_path / name
    # CRLF — как в настоящих файлах Ancestry.
    #
    # newline="" здесь обязателен: без него текстовый режим на Windows
    # переводит каждый '\n' в '\r\n', и записанные нами '\r\n'
    # превращаются в '\r\r\n'. Файл после этого читается как удвоенное
    # число строк (пустая строка после каждой настоящей), и индекс
    # строки заголовка уезжает: 18 -> 36. На Linux этого не видно —
    # перевода нет, поэтому тест зелёный там и красный на Windows.
    path.write_text((ANCESTRY_HEADER + "\n" + body).replace("\n", "\r\n"),
                    encoding="utf-8", newline="")
    return path


# ---------------------------------------------------------------------------
# Нормализация хромосом
# ---------------------------------------------------------------------------
def test_chrom_codes_normalized():
    assert _normalize_chrom("1") == "1"
    assert _normalize_chrom("23") == "X"
    assert _normalize_chrom("24") == "Y"
    assert _normalize_chrom("26") == "MT"


def test_code_25_is_par_of_x_not_mt():
    """
    Регрессия: у AncestryDNA код 25 — это псевдоаутосомная область X
    (PAR), а у MyHeritage тот же код 25 означает MT. Если кто-то
    переиспользует CHROM_NORMALIZE из myheritage_v5.py, ~36 позиций X
    молча уедут в митохондриальный контиг.
    """
    assert _normalize_chrom("25") == "X"
    assert CHROM_NORMALIZE["25"] == "X"

    from adapters.myheritage_v5 import CHROM_NORMALIZE as MH_MAP
    assert MH_MAP["25"] == "MT", "предпосылка теста: у MyHeritage 25 = MT"
    assert CHROM_NORMALIZE["25"] != MH_MAP["25"]


# ---------------------------------------------------------------------------
# Поиск заголовка
# ---------------------------------------------------------------------------
def test_header_found_after_comment_block(tmp_path):
    """_find_header() возвращает (индекс_строки, оформление) — оформлений
    два, см. adapters/ancestry_v2.py и tests/test_ancestry_convert.py."""
    path = _write(tmp_path, "rs1\t1\t100\tA\tA\n")
    assert _find_header(path) == (18, LAYOUT_RAW)


def test_foreign_format_rejected(tmp_path):
    path = tmp_path / "ftdna.csv"
    path.write_text("RSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,AA\n", encoding="utf-8")
    with pytest.raises(AncestryFormatError):
        _find_header(path)


# ---------------------------------------------------------------------------
# Разбор генотипов
# ---------------------------------------------------------------------------
def test_two_allele_columns_are_joined_and_resolved(tmp_path):
    path = _write(tmp_path, textwrap.dedent("""\
        rs_hom_ref\t1\t100\tA\tA
        rs_het\t1\t200\tA\tG
        rs_hom_alt\t1\t300\tG\tG
        """))
    reference = _StubReference({}, default="A")
    result = parse_ancestry_v2(path, reference)

    by_rsid = {v.rsid: v for v in result.variants}
    assert by_rsid["rs_hom_ref"].gt == "0/0"
    assert by_rsid["rs_het"].gt == "0/1"
    assert by_rsid["rs_het"].alt == "G"
    assert by_rsid["rs_hom_alt"].gt == "1/1"
    assert result.total_measured == 3


def test_zero_allele_is_no_call(tmp_path):
    """У Ancestry пропуск — аллель '0', а не '--'."""
    path = _write(tmp_path, textwrap.dedent("""\
        rs_ok\t1\t100\tA\tA
        rs_missing\t1\t200\t0\t0
        rs_half\t1\t300\tA\t0
        """))
    result = parse_ancestry_v2(path, _StubReference({}))
    assert result.missing == 2
    assert result.total_measured == 1
    assert [v.rsid for v in result.variants] == ["rs_ok"]


def test_indels_counted_as_invalid_not_crash(tmp_path):
    """I/D в колонках аллелей — не ACGT, уходят в invalid_codes (как в FTDNA)."""
    path = _write(tmp_path, textwrap.dedent("""\
        rs_ok\t1\t100\tA\tA
        rs_ins\t1\t200\tI\tI
        rs_del\t1\t300\tD\tD
        """))
    result = parse_ancestry_v2(path, _StubReference({}))
    assert result.invalid_codes == 2
    assert len(result.variants) == 1


def test_duplicate_positions_dropped(tmp_path):
    """В реальном V2.0 таких повторов ~650; build_vcf() падает на дублях."""
    path = _write(tmp_path, textwrap.dedent("""\
        rs_first\t1\t100\tA\tA
        rs_second\t1\t100\tA\tA
        """))
    result = parse_ancestry_v2(path, _StubReference({}))
    assert result.duplicate_positions == 1
    assert [v.rsid for v in result.variants] == ["rs_first"]


def test_malformed_row_counted(tmp_path):
    path = _write(tmp_path, "rs_ok\t1\t100\tA\tA\nrs_short\t1\t200\tA\n")
    result = parse_ancestry_v2(path, _StubReference({}))
    assert result.malformed_rows == 1


def test_broad_signature_includes_no_calls(tmp_path):
    """
    Широкая сигнатура (Задача D) — отпечаток ДИЗАЙНА чипа, поэтому в неё
    должны попадать и позиции, где у этого человека не прочиталось.
    Иначе переиспользование доноров между людьми не сработает никогда.
    """
    path = _write(tmp_path, "rs_ok\t1\t100\tA\tA\nrs_missing\t1\t200\t0\t0\n")
    result = parse_ancestry_v2(path, _StubReference({}))
    assert result.signature_positions_broad == [("1", 100), ("1", 200)]
    assert result.chip_signature != result.chip_signature_broad


# ---------------------------------------------------------------------------
# Подключение к пайплайну
# ---------------------------------------------------------------------------
def test_source_registered():
    assert "ancestry" in pipeline.SOURCES
    assert pipeline.SOURCES["ancestry"]["parser"] is parse_ancestry_v2
    assert pipeline._needs_reference("ancestry")
    assert pipeline._supports_liftover("ancestry")


def test_autodetect_prefers_ancestry_over_myheritage(tmp_path):
    """
    Регрессия: у Ancestry 18 ведущих '#'-строк, а правило MyHeritage
    срабатывает на '10+ комментариев подряд'. Без отдельного правила,
    проверяемого РАНЬШЕ, любой Ancestry-файл определялся бы как
    myheritage с уверенностью 0.9 — и уехал бы в чужой парсер.
    """
    path = _write(tmp_path, "rs1\t1\t100\tA\tA\n")
    source, confidence = pipeline.detect_source_from_file(path)
    assert source == "ancestry"
    assert confidence >= 0.9


def test_autodetect_by_column_header_without_ancestry_banner(tmp_path):
    """Первая строка-баннер могла быть срезана — ловим по строке колонок."""
    path = tmp_path / "raw.txt"
    body = "\n".join(
        [f"#comment {i}" for i in range(12)]
        + ["rsid\tchromosome\tposition\tallele1\tallele2", "rs1\t1\t100\tA\tA"]
    )
    path.write_text(body + "\n", encoding="utf-8")
    assert pipeline.detect_source_from_file(path)[0] == "ancestry"


def test_autodetect_still_recognises_myheritage(tmp_path):
    """Новое правило не должно перехватывать чужие файлы."""
    path = tmp_path / "mh.csv"
    body = "\n".join(
        [f"# comment {i}" for i in range(12)]
        + ['"RSID","CHROMOSOME","POSITION","RESULT"', '"rs1","1","100","AA"']
    )
    path.write_text(body + "\n", encoding="utf-8")
    assert pipeline.detect_source_from_file(path)[0] == "myheritage"


def test_autodetect_still_recognises_ftdna(tmp_path):
    path = tmp_path / "ftdna.csv"
    path.write_text("RSID,CHROMOSOME,POSITION,RESULT\nrs1,1,100,AA\n", encoding="utf-8")
    assert pipeline.detect_source_from_file(path)[0] == "ftdna"
