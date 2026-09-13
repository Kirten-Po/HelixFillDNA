"""
tests/test_rsq_tuner.py

Подбор порога Rsq под целевую заполняемость — прямое решение случая, где
приёмку не прошли из-за недобора 3,4 п.п. Здесь легко сделать вид, что
цель достигнута: достаточно опустить порог в пол. Поэтому тесты в первую
очередь про границы, за которые модуль заходить не должен.
"""
from __future__ import annotations

import pytest

from core import rsq_tuner as rt


def _fixture(n_measured=30, n_imputed=60, n_total=100):
    """Трафарет из n_total строк: часть закрыта измерениями, часть импутацией
    с равномерно распределённым Rsq от 0 до 1."""
    skeleton = [f"1_{i}" for i in range(n_total)]
    measured = {f"1_{i}" for i in range(n_measured)}
    imputed_rsq = {
        f"1_{i}": (i - n_measured) / n_imputed
        for i in range(n_measured, n_measured + n_imputed)
    }
    return skeleton, imputed_rsq, measured


def test_picks_highest_threshold_meeting_target():
    """
    "Самый высокий порог, дающий цель" — принципиально: цель в том, чтобы
    пройти приёмку, а не набить файл вызовами похуже.
    """
    skeleton, imputed_rsq, measured = _fixture()
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=80.0)
    assert r.target_reached
    assert r.achieved_call_rate >= 80.0
    # Шаг выше цель уже не выполняется — значит взят действительно максимум.
    higher = [p for p in r.curve if p.threshold > r.chosen_threshold + 1e-9]
    assert all(p.call_rate < 80.0 for p in higher)


def test_never_goes_below_floor():
    skeleton, imputed_rsq, measured = _fixture()
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=99.9)
    assert not r.target_reached, "недостижимую цель нельзя объявлять достигнутой"
    assert r.chosen_threshold >= rt.RSQ_FLOOR


def test_unreachable_target_reports_max():
    skeleton, imputed_rsq, measured = _fixture()
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=99.0)
    assert r.achieved_call_rate == pytest.approx(r.max_call_rate)
    text = rt.format_result(r)
    assert "НЕ достигнута" in text
    assert "шум" in text, "надо объяснить, почему порог не опущен ещё ниже"


def test_measured_positions_are_threshold_independent():
    """
    Реальные измерения чипа перекрывают импутацию (merge_dictionaries) и
    заполнены при любом пороге — иначе кривая была бы занижена.
    """
    skeleton, imputed_rsq, measured = _fixture(n_measured=40, n_imputed=0,
                                               n_total=100)
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=None)
    assert r.achieved_call_rate == pytest.approx(40.0)
    assert all(p.call_rate == pytest.approx(40.0) for p in r.curve)


def test_price_is_reported():
    """Заполняемость покупается за качество, и это должно быть написано."""
    skeleton, imputed_rsq, measured = _fixture()
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=85.0)
    assert r.below_standard_used > 0
    assert r.gained_pp > 0
    assert 0.0 < r.median_rsq_used <= 1.0
    text = rt.format_result(r)
    assert "Чем заплачено" in text


def test_no_target_keeps_standard_threshold():
    skeleton, imputed_rsq, measured = _fixture()
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=None,
                baseline=0.30)
    assert r.chosen_threshold == pytest.approx(0.30, abs=rt.RSQ_STEP)
    assert r.target_reached


def test_filter_by_threshold_matches_curve():
    """
    Кривая — предсказание, filter_by_threshold — то, что реально уйдёт в
    сборку. Они обязаны сходиться, иначе подобранный порог даст не ту
    заполняемость, которую пообещали.
    """
    skeleton, imputed_rsq, measured = _fixture()
    genotypes = {k: "AG" for k in imputed_rsq}
    r = rt.tune(skeleton, imputed_rsq, measured, target_call_rate=80.0)
    kept = rt.filter_by_threshold(genotypes, imputed_rsq, r.chosen_threshold)
    merged = set(kept) | measured
    predicted = 100.0 * sum(1 for k in skeleton if k in merged) / len(skeleton)
    assert predicted == pytest.approx(r.achieved_call_rate)


def test_duplicate_skeleton_rows_counted_twice():
    """
    В трафарете 23andMe встречаются повторяющиеся позиции, и
    validate_output() считает СТРОКИ, а не уникальные ключи. Кривая
    должна считать так же, иначе предсказание разойдётся с проверкой.
    """
    skeleton = ["1_1", "1_1", "1_2", "1_3"]
    r = rt.tune(skeleton, {"1_2": 0.9}, {"1_1"}, target_call_rate=None)
    assert r.total_rows == 4
    assert r.achieved_call_rate == pytest.approx(75.0)


def test_empty_input_does_not_crash():
    r = rt.tune([], {}, set(), target_call_rate=90.0)
    assert r.achieved_call_rate == 0.0
    assert not r.target_reached


# ---------------------------------------------------------------------------
# Метрика GP: тот же алгоритм, другие пол, точка отсчёта и подпись
# ---------------------------------------------------------------------------
def test_tune_takes_floor_and_baseline_from_metric():
    """
    Пол и точка отсчёта не задаются вручную: для GP они другие (0.50/0.90
    против 0.10/0.30), и передавать их на каждом вызове означало бы рано
    или поздно подобрать порог по чужой шкале.
    """
    from core import rsq_tuner

    keys = [f"1_{i}" for i in range(100)]
    # Половина позиций — уверенные вызовы, половина — размазанные.
    gp = {f"1_{i}": (0.99 if i % 2 == 0 else 0.62) for i in range(100)}

    r = rsq_tuner.tune(keys, gp, set(), target_call_rate=None, metric="gp")
    assert r.metric == "gp"
    assert r.baseline_threshold == rsq_tuner.GP_STANDARD
    # При пороге 0.90 проходит ровно половина.
    assert abs(r.baseline_call_rate - 50.0) < 1e-6
    # Пол сетки — GP_FLOOR, а не RSQ_FLOOR: на 0.50 проходят все.
    assert abs(r.max_call_rate - 100.0) < 1e-6
    assert r.curve[0].threshold == rsq_tuner.GP_FLOOR


def test_tune_reaches_target_by_lowering_gp_threshold():
    from core import rsq_tuner

    keys = [f"1_{i}" for i in range(100)]
    gp = {f"1_{i}": (0.99 if i % 2 == 0 else 0.62) for i in range(100)}

    r = rsq_tuner.tune(keys, gp, set(), target_call_rate=90.0, metric="gp")
    assert r.target_reached
    assert r.chosen_threshold <= 0.62 + 1e-9
    assert r.achieved_call_rate >= 90.0


def test_format_result_labels_the_metric():
    """
    Подпись обязана называть метрику: 0,30 по Rsq и 0,90 по GP стоят в
    разных шкалах, и отчёт без имени метрики нечитаем.
    """
    from core import rsq_tuner

    keys = [f"1_{i}" for i in range(10)]
    gp = {k: 0.95 for k in keys}
    text = rsq_tuner.format_result(
        rsq_tuner.tune(keys, gp, set(), metric="gp")
    )
    assert "max(GP)" in text
    assert "ПОДБОР ПОРОГА max(GP)" in text

    rsq = {k: 0.8 for k in keys}
    text_rsq = rsq_tuner.format_result(
        rsq_tuner.tune(keys, rsq, set(), metric="rsq")
    )
    assert "Rsq" in text_rsq and "max(GP)" not in text_rsq


def test_metrics_row_records_which_metric():
    from core import rsq_tuner

    keys = [f"1_{i}" for i in range(10)]
    row = rsq_tuner.metrics_row(
        rsq_tuner.tune(keys, {k: 0.95 for k in keys}, set(), metric="gp")
    )
    assert row["quality_metric"] == "gp"
    assert row["rsq_chosen"] == rsq_tuner.GP_STANDARD
