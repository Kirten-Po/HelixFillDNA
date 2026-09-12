"""
tests/test_blank_y.py

Y-хромосома в итоговом файле всегда остаётся без вызовов.

Предыстория: Y не импутируется — на сервер уходят только 1-22 и X
(core/pure_python_core.py::UPLOAD_CHROMS), поэтому в собранный файл она
попадает единственным путём, прямыми измерениями чипа. Новые экспорты
MyHeritage такие измерения содержат (в трафарете v5 под Y отведено 3325
строк), и блок Y — он идёт предпоследним, между X и MT, то есть почти в
самом конце файла — внезапно оказывался заполнен.

Два свойства, которые тесты обязаны удержать:
  * строки Y НЕ удаляются, а становятся "--": итоговый файл обязан
    построчно повторять трафарет, это проверяет validate_output();
  * очистка происходит ДО подбора порога Rsq, иначе предсказанная
    заполняемость разойдётся с фактической.
"""
from __future__ import annotations

import pytest

from core import rsq_tuner as rt
from template.assembler import (
    BLANKED_CHROMS, assemble_final, blank_chromosomes, merge_dictionaries,
    validate_output,
)
from template.skeleton import SkeletonRow


def _skeleton():
    """Мини-трафарет с тем же порядком блоков, что у настоящего: ... X, Y, MT."""
    def row(rsid, chrom, pos):
        return SkeletonRow(rsid=rsid, chrom=chrom, pos=pos,
                           raw_line=f"{rsid}\t{chrom}\t{pos}\tAG\n")

    rows = [row(f"rs{i}", "1", 1000 + i) for i in range(5)]
    rows += [row(f"rsx{i}", "X", 2000 + i) for i in range(2)]
    rows += [row(f"rsy{i}", "Y", 3000 + i) for i in range(3)]
    rows += [row(f"rsm{i}", "MT", 40 + i) for i in range(2)]
    return rows


def _all_genotypes(rows):
    return {f"{r.chrom}_{r.pos}": "AG" for r in rows}


def test_only_y_is_removed():
    rows = _skeleton()
    kept, removed = blank_chromosomes(_all_genotypes(rows))
    assert removed == 3
    assert not any(k.startswith("Y_") for k in kept)
    # X и MT трогать нельзя: MT — единственный источник mtDNA-гаплогруппы,
    # X импутируется и является половиной смысла всей затеи.
    assert any(k.startswith("X_") for k in kept)
    assert any(k.startswith("MT_") for k in kept)


def test_chr_prefix_is_handled():
    """
    После обратного лифтовера с панели TopMed ключи проходят через
    несколько рук — префикс "chr" не должен спасать Y от очистки.
    """
    kept, removed = blank_chromosomes({"chrY_3000": "AG", "chr1_100": "AA"})
    assert removed == 1
    assert list(kept) == ["chr1_100"]


def test_empty_and_idempotent():
    assert blank_chromosomes({}) == ({}, 0)
    once, _ = blank_chromosomes(_all_genotypes(_skeleton()))
    twice, removed = blank_chromosomes(once)
    assert removed == 0 and twice == once


def test_y_rows_stay_in_file_as_no_calls(tmp_path):
    """
    Главное свойство: строк в файле столько же, сколько в трафарете, и
    все Y-строки на своих местах со значением "--".
    """
    rows = _skeleton()
    template = tmp_path / "template_v3.txt"
    template.write_text(
        "# rsid\tchromosome\tposition\tgenotype\n"
        + "\n".join(f"{r.rsid}\t{r.chrom}\t{r.pos}\tAG" for r in rows) + "\n",
        encoding="utf-8",
    )

    genotypes, _ = blank_chromosomes(_all_genotypes(rows))
    out = tmp_path / "result.txt"
    assemble_final(rows, genotypes, out, format_version="v3", template_path=template)

    data = [l for l in out.read_text(encoding="utf-8").splitlines()
            if not l.startswith("#")]
    assert len(data) == len(rows), "структура файла обязана повторять трафарет"

    y_lines = [l.split("\t") for l in data if l.split("\t")[1] == "Y"]
    assert len(y_lines) == 3
    assert all(parts[3] == "--" for parts in y_lines)
    # Порядок сохранён: Y по-прежнему между X и MT, а не выброшен в конец.
    chroms = [l.split("\t")[1] for l in data]
    assert chroms.index("Y") > chroms.index("X")
    assert chroms.index("MT") > chroms.index("Y")

    validation = validate_output(out, template, "v3")
    assert validation.structure_identical
    assert validation.call_rate == pytest.approx(100.0 * 9 / 12), (
        "заполнены все строки, кроме трёх Y"
    )


def test_measured_y_does_not_leak_through_merge():
    """
    Импутация Y не даёт (Y не уходит на сервер), поэтому единственный
    источник — measured. Проверяем, что после слияния Y не появляется.
    """
    imputed = {"1_1000": "AA"}
    measured = {"1_1001": "AG", "Y_3000": "TT", "MT_40": "CC"}
    measured, removed = blank_chromosomes(measured)
    assert removed == 1
    merged = merge_dictionaries(imputed, measured)
    assert not any(k.startswith("Y_") for k in merged)
    assert "MT_40" in merged


def test_rsq_curve_matches_file_after_blanking():
    """
    Очистка Y обязана происходить ДО подбора порога: иначе подобранный
    порог обещал бы заполняемость, которой в файле не будет.
    """
    rows = _skeleton()
    skeleton_keys = [f"{r.chrom}_{r.pos}" for r in rows]
    measured = {f"{r.chrom}_{r.pos}" for r in rows if r.chrom in ("1", "Y")}

    измерено_без_y, _ = blank_chromosomes({k: "AG" for k in measured})
    r = rt.tune(skeleton_keys, {}, set(измерено_без_y), target_call_rate=None)
    # 5 аутосомных строк из 12; Y в знаменателе остаётся, в числителе нет.
    assert r.total_rows == len(rows)
    assert r.achieved_call_rate == pytest.approx(100.0 * 5 / len(rows))


def test_blanked_chroms_is_exactly_y():
    """
    Осознанный выбор, а не случайность: MT остаётся (по нему определяется
    материнская гаплогруппа), X остаётся (ради него делалась импутация X).
    """
    assert BLANKED_CHROMS == frozenset({"Y"})
