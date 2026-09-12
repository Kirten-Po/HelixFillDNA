"""
tests/test_genotek_template.py

Трафарет `template_genotek.txt` собран из пересечения позиций трёх
реальных VCF от Генотека (брат, сестра и её муж). Смысл его в том, что
сравнение с родными данными Генотека идёт по ПЕРЕСЕЧЕНИЮ наборов
позиций, а там наши обычные трафареты покрывают чип Генотека плохо:
v3 — 30,7 %, v5 — 91,1 %, этот — 98,5 % (по аутосомам и X; митохондрия
исключена, см. test_template_covers_only_autosomes_and_x).

Ловушка, ради которой этот файл тестов существует отдельно: переносы
строк. У нового формата они CRLF, как у v5, но код в двух местах
проверял `format_version == "v5"` жёстко — и любой новый формат молча
получал LF, после чего validate_output() объявлял собственный же файл
невалидным.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from template.assembler import (
    CRLF_FORMATS, assemble_final, line_ending_for, validate_output,
)
from template.skeleton import SkeletonRow, extract_skeleton

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = PROJECT_ROOT / "samples" / "template_genotek.txt"
needs_template = pytest.mark.skipif(
    not TEMPLATE.is_file(), reason="template_genotek.txt не установлен",
)


# ---------------------------------------------------------------------------
# Переносы строк
# ---------------------------------------------------------------------------
def test_line_endings_are_driven_by_a_list_not_by_a_hardcoded_v5():
    assert line_ending_for("v3") == "\n"
    assert line_ending_for("v5") == "\r\n"
    assert line_ending_for("genotek") == "\r\n", (
        "у трафарета Генотека оформление v5-шное — CRLF"
    )
    assert {"v5", "genotek"} <= CRLF_FORMATS


def _rows():
    def row(rsid, chrom, pos):
        return SkeletonRow(rsid=rsid, chrom=chrom, pos=pos,
                           raw_line=f"{rsid}\t{chrom}\t{pos}\t--\r\n")
    return [row("rs1", "1", 100), row("rs2", "1", 200), row("rsm", "MT", 10)]


def test_genotek_format_round_trips_through_validation(tmp_path):
    """
    Главная регрессия: собранный в формате genotek файл обязан проходить
    собственную же проверку. С прежней жёсткой проверкой на "v5" он
    получал LF и падал на «CRLF/LF не соответствует формату».
    """
    rows = _rows()
    template = tmp_path / "template_genotek.txt"
    template.write_text(
        "# rsid\tchromosome\tposition\tgenotype\r\n"
        + "".join(f"{r.rsid}\t{r.chrom}\t{r.pos}\t--\r\n" for r in rows),
        encoding="utf-8", newline="",
    )
    out = tmp_path / "out.txt"
    assemble_final(rows, {"1_100": "AG"}, out, format_version="genotek",
                   template_path=template)

    raw = out.read_bytes()
    assert raw.count(b"\r\n") == len(raw.splitlines()), "все строки должны быть CRLF"

    v = validate_output(out, template, "genotek")
    assert v.is_valid, v.errors
    assert v.structure_identical


# ---------------------------------------------------------------------------
# Сам файл трафарета
# ---------------------------------------------------------------------------
@needs_template
def test_template_parses_and_is_well_formed():
    rows = extract_skeleton(TEMPLATE, autosomes_only=False)
    assert len(rows) == 621_566

    order = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "MT": 25}
    keys = [(order[r.chrom], r.pos) for r in rows]
    assert keys == sorted(keys), "трафарет должен быть отсортирован"
    assert len({(r.chrom, r.pos) for r in rows}) == len(rows), "позиции уникальны"
    assert len({r.rsid for r in rows}) == len(rows), "rsID уникальны"
    assert all(r.rsid.startswith("rs") for r in rows), (
        "записи без rsID отброшены при сборке трафарета"
    )


@needs_template
def test_template_covers_only_autosomes_and_x():
    """
    Ни Y, ни MT в трафарете быть не должно — по разным причинам.

    Y: у чипа Генотека её нет вовсе (в их VCF только chr1..chr22, chrX,
    chrM), и мы её в итоговом файле всё равно очищаем.

    MT: координаты митохондрии у Генотека — в системе hg19/Yoruba
    (NC_001807), а не rCRS (NC_012920), на которой сидят 23andMe и все
    потребительские чипы. Сверка по rsID с чипом Genera дала 13
    совпадений из 842 при систематическом сдвиге +1/+2; контроль
    Genera против template_v5 на тех же данных — 173 совпадения из 174.
    Попасть в такие строки не может ни один источник, поэтому MT
    исключена из трафарета целиком.
    """
    chroms = {r.chrom for r in extract_skeleton(TEMPLATE, autosomes_only=False)}
    assert "Y" not in chroms
    assert "MT" not in chroms
    assert "X" in chroms
    assert chroms == {str(i) for i in range(1, 23)} | {"X"}


@needs_template
def test_template_file_is_crlf():
    raw = TEMPLATE.read_bytes()
    assert b"\r\n" in raw
    assert raw.count(b"\n") == raw.count(b"\r\n"), "смешанных переносов быть не должно"


@needs_template
def test_template_is_registered_in_the_app():
    """
    Файл мало положить в samples/ — программа должна знать про него по
    ключу формата, иначе он выбирается только руками через «Обзор».
    """
    # gui/app.py тянет tkinter/customtkinter — на машинах без GUI-стека
    # (в т.ч. в этом контейнере) тест пропускается, но на Windows и в CI
    # он реально проверяет регистрацию формата.
    pytest.importorskip("tkinter")
    pytest.importorskip("customtkinter")
    import gui.app as app  # noqa: PLC0415

    assert app.SAMPLE_TEMPLATE_NAMES["genotek"] == "template_genotek.txt"
    assert app.FORMAT_LABELS["genotek"].startswith("genotek"), (
        "ключ формата опознаётся по началу подписи"
    )
    assert app._find_sample_template("genotek") is not None
