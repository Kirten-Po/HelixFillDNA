"""
tests/test_panel_advisor.py

Проверяем главное свойство модуля: панель выбирается по СОСТАВУ ЧИПА,
а не по качеству файла. И конкретно — что решение принимается по редкому
хвосту, а доля позиций вне 1000G его только усиливает.

Почему это важно закрепить тестом: на реальных данных две метрики
расходятся в разные стороны (FTDNA — редкий хвост 13,4 %, вне 1000G
2,7 %; AncestryDNA — 3,3 % и 9,9 %), и если бы доля отсутствующих
работала самостоятельным триггером, AncestryDNA молча уезжал бы на
TopMed с лифтовером без всякой пользы.
"""
from __future__ import annotations

import gzip

import pytest

from core import panel_advisor as pa


def _make_donor_dir(tmp_path, afs, chip_positions=None):
    """
    Мини-кэш доноров: один kgp_sub_1.vcf.gz с заданными EUR_AF и список
    позиций чипа рядом (тот самый знаменатель "18 837 из 724 937").
    """
    d = tmp_path / "donors" / "ftdna" / "hrc"
    d.mkdir(parents=True)
    lines = [
        "##fileformat=VCFv4.1",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    for i, af in enumerate(afs):
        info = "AC=1;AN=40" if af is None else f"AC=1;AN=40;EUR_AF={af}"
        lines.append(f"1\t{1000 + i}\trs{i}\tA\tG\t100\tPASS\t{info}\tGT\t0|0")
    with gzip.open(d / "kgp_sub_1.vcf.gz", "wt", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    total = chip_positions if chip_positions is not None else len(afs)
    (d / "ftdna_pos.txt").write_text(
        "\n".join(f"1\t{1000 + i}" for i in range(total)) + "\n", encoding="utf-8"
    )
    return d


def test_rare_tail_counted_as_maf_not_af(tmp_path):
    """
    EUR_AF=0.998 — это редкий маркер (минорная частота 0,2 %), а не
    сплошь распространённый. Считать надо минорную частоту.
    """
    d = _make_donor_dir(tmp_path, [0.998, 0.5, 0.4, 0.6])
    comp = pa.analyse_chip(d, bcftools_path=None)
    assert comp.matched == 4
    assert comp.rare == 1


def test_missing_from_1000g_uses_chip_denominator(tmp_path):
    d = _make_donor_dir(tmp_path, [0.5] * 90, chip_positions=100)
    comp = pa.analyse_chip(d, bcftools_path=None)
    assert comp.chip_positions == 100
    assert comp.matched == 90
    assert comp.missing == 10
    assert comp.missing_pct == pytest.approx(10.0)


def test_heavy_rare_tail_recommends_topmed(tmp_path):
    afs = [0.001] * 20 + [0.4] * 80          # 20 % редких
    comp = pa.analyse_chip(_make_donor_dir(tmp_path, afs), bcftools_path=None)
    rec = pa.recommend_panel(comp)
    assert rec.panel == pa.PANEL_TOPMED
    assert "MAF" in rec.headline
    assert rec.counter_reasons, "цену за TopMed (лифтовер) надо назвать явно"


def test_light_rare_tail_recommends_hrc_even_with_many_missing(tmp_path):
    """Случай AncestryDNA: мало редких, много отсутствующих в 1000G."""
    afs = [0.001] * 3 + [0.4] * 97           # 3 % редких
    comp = pa.analyse_chip(_make_donor_dir(tmp_path, afs, chip_positions=112),
                           bcftools_path=None)
    assert comp.missing_pct > pa.MISSING_PCT_THRESHOLD
    rec = pa.recommend_panel(comp)
    assert rec.panel == pa.PANEL_HRC, (
        "доля отсутствующих в 1000G — вспомогательный сигнал, а не триггер"
    )
    assert any("отсутствуют в 1000G" in c for c in rec.counter_reasons), (
        "но умолчать о ней тоже нельзя"
    )


def test_non_european_sample_forces_topmed(tmp_path):
    """
    Неевропейское происхождение перевешивает состав чипа: подвыборка HRC
    почти сплошь европейская.
    """
    afs = [0.4] * 100
    comp = pa.analyse_chip(_make_donor_dir(tmp_path, afs), bcftools_path=None)
    rec = pa.recommend_panel(comp, sample_is_european=False)
    assert rec.panel == pa.PANEL_TOPMED


def test_missing_af_is_not_counted_as_rare(tmp_path):
    """Запись без EUR_AF — неизвестность, а не редкость."""
    comp = pa.analyse_chip(_make_donor_dir(tmp_path, [None] * 10),
                           bcftools_path=None)
    assert comp.rare == 0
    assert comp.matched == 10


def test_no_cache_raises(tmp_path):
    empty = tmp_path / "nothing"
    empty.mkdir()
    with pytest.raises(pa.PanelAdvisorError):
        pa.analyse_chip(empty, bcftools_path=None)


def test_find_donor_cache_picks_any_panel(tmp_path):
    """
    Частоты 1000 Genomes от выбора панели не зависят, поэтому состав чипа
    считается по любому уже скачанному кэшу — иначе получилась бы курица
    и яйцо: чтобы посоветовать панель, надо скачать доноров для панели.
    """
    d = _make_donor_dir(tmp_path, [0.5] * 5)
    found = pa.find_donor_cache(tmp_path / "donors", "ftdna")
    assert found == d
    assert pa.find_donor_cache(tmp_path / "donors", "myheritage") is None


def test_falls_back_when_bcftools_unusable(tmp_path):
    """
    bcftools может отсутствовать, не запуститься или отказаться от файла,
    в заголовке которого нет описания тега EUR_AF. Ни один из этих
    случаев не должен выбрасывать файл из статистики: тогда молча уехал
    бы вниз знаменатель "найдено в 1000G", а за ним и доля отсутствующих
    позиций — отчёт показал бы неправду вместо того, чтобы просто читать
    помедленнее.
    """
    afs = [0.001] * 5 + [0.4] * 15
    d = _make_donor_dir(tmp_path, afs)

    reference = pa.analyse_chip(d, bcftools_path=None)
    fallback = pa.analyse_chip(d, bcftools_path="/несуществующий/bcftools")

    assert fallback.matched == reference.matched == 20
    assert fallback.rare == reference.rare == 5
