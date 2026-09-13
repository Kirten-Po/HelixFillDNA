"""
core/rsq_tuner.py
Подбор порога качества под ЦЕЛЕВУЮ заполняемость итогового файла.

Метрик две — Rsq и max(GP), см. METRIC_SETTINGS и параметр metric в
tune(). Алгоритм у них общий: значения сортируются один раз, и кривая
"порог -> заполняемость" считается бинарным поиском. Различаются только
пол, точка отсчёта и подпись — имя модуля историческое.

Зачем это раньше автовыбора панели
----------------------------------
Отказ приёмки (случай GAEVA) был не из-за панели, а из-за недобора 3,4
процентных пункта по заполняемости. Смена панели — это новая закачка
доноров, новая очередь на Michigan Imputation Server и лифтовер; подбор
порога — это пересчёт по УЖЕ СКАЧАННЫМ дозам, стоящий секунды и не
требующий повторного задания на MIS. Эффект прямой и проверяется на уже
имеющихся данных.

Как это работает
----------------
Дозы читаются ОДИН раз, с низким порогом-полом (RSQ_FLOOR), причём
template/assembler.py::load_imputed_genotypes() по параметру rsq_out
попутно отдаёт Rsq каждой принятой позиции. Дальше зависимость
"порог -> заполняемость" считается арифметикой, без единого повторного
запуска bcftools: позиции трафарета, закрытые РЕАЛЬНЫМИ измерениями
чипа, от порога не зависят вовсе, а импутированные отсортированы по Rsq
один раз и считаются двоичным поиском.

Что здесь принципиально
-----------------------
Заполняемость покупается за качество, и это должно быть написано прямо,
а не спрятано. Поэтому TuningResult несёт не только выбранный порог, но
и цену: сколько позиций прошло с Rsq ниже стандартного 0,30, какая у
использованных позиций медиана Rsq и на сколько пунктов заполняемость
выросла относительно стандартного порога. Приложение НЕ занижает порог
молча и не занижает его ниже RSQ_FLOOR ни при каких целях: если цель
недостижима, оно так и говорит, а не подгоняет цифру любой ценой.
"""
from __future__ import annotations

import bisect
import logging
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

logger = logging.getLogger(__name__)

#: Ниже этого порога подбор не опускается ни при какой цели. Rsq 0,10 —
#: это уже почти шум: генотип в такой позиции определён немногим лучше,
#: чем подбрасыванием монетки со смещением к частой аллели. Заполнять
#: файл такими вызовами ради процента в отчёте — обман приёмки, а не
#: решение задачи.
RSQ_FLOOR = 0.10

#: Стандартный порог MIS — точка отсчёта, относительно которой считается
#: "сколько качества заплачено".
RSQ_STANDARD = 0.30

#: Шаг сетки перебора.
RSQ_STEP = 0.01

# ---------------------------------------------------------------------------
# Вторая метрика качества: max(GP)
#
# Кривая "порог -> заполняемость" в build_curve() и tune() не знает, ЧТО за
# число ей дали: она сортирует значения и режет их бинарным поиском. Поэтому
# перевод на GP — это не новый алгоритм, а другие пол, точка отсчёта и
# подпись. Важно лишь, чтобы словарь качества, приходящий из
# template/assembler.py (параметр rsq_out), содержал ОДНУ метрику: шкалы у
# Rsq и GP разные, и смешанный словарь дал бы бессмысленный порог.
# ---------------------------------------------------------------------------

#: Ниже этого порога подбор по GP не опускается. 0,50 — точка, где самый
#: вероятный генотип перестаёт быть вероятнее всех остальных вместе взятых.
#: Ровно тот же смысл, что у RSQ_FLOOR: дальше уже не вызовы, а шум.
GP_FLOOR = 0.50

#: Стандартный порог GP — точка отсчёта, относительно которой считается
#: "сколько качества заплачено". Выбран по форме кривой "порог ->
#: заполняемость", замеренной на реальном прогоне (FTDNA -> трафарет
#: genotek, весь геном, 503 донора; чистый эффект перехода с Rsq 0,30):
#:
#:      0,99  +0,48 п.п.        0,85  +2,32 п.п.
#:      0,95  +1,76 п.п.        0,80  +2,42 п.п.
#:      0,90  +2,15 п.п.        0,70  +2,54 п.п.
#:
#: ⚠ Сами величины верны только для ТОГО прогона: размер выигрыша сильно
#: зависит от числа доноров и от покрытия трафарета чипом (см. докстринг
#: template/assembler.py::load_imputed_genotypes()). Здесь важна не
#: величина, а ФОРМА: колено около 0,85-0,90, ниже каждая следующая
#: ступень добавляет всё меньше и состоит из всё менее уверенных
#: вызовов. 0,90 — консервативный край колена; он же принят в локальном
#: семейном пайплайне (Beagle, gp=true).
GP_STANDARD = 0.90

#: Пол, точка отсчёта и подпись для каждой метрики. Всё, чем метрики
#: различаются в этом модуле.
METRIC_SETTINGS = {
    "rsq": (RSQ_FLOOR, RSQ_STANDARD, "Rsq"),
    "gp": (GP_FLOOR, GP_STANDARD, "max(GP)"),
}


def settings_for(metric: str) -> tuple[float, float, str]:
    """(пол, точка отсчёта, подпись) для метрики; неизвестная -> Rsq."""
    return METRIC_SETTINGS.get(metric, METRIC_SETTINGS["rsq"])


@dataclass(frozen=True)
class CurvePoint:
    threshold: float
    call_rate: float          # % заполненных строк трафарета
    imputed_used: int         # сколько строк закрыто импутацией
    measured_used: int        # сколько закрыто реальными измерениями


@dataclass
class TuningResult:
    target_call_rate: Optional[float]
    chosen_threshold: float
    achieved_call_rate: float
    baseline_threshold: float = RSQ_STANDARD
    baseline_call_rate: float = 0.0
    target_reached: bool = True
    max_call_rate: float = 0.0        # что даёт сам пол RSQ_FLOOR
    below_standard_used: int = 0      # позиции с Rsq < RSQ_STANDARD в итоге
    median_rsq_used: float = 0.0
    total_rows: int = 0
    curve: list[CurvePoint] = field(default_factory=list)
    #: Какой метрикой подбирали — "rsq" или "gp". Нужна и для подписей, и
    #: для runs_metrics.csv: без неё столбец с порогом нечитаем, потому
    #: что 0,30 по Rsq и 0,90 по GP стоят в разных шкалах.
    metric: str = "rsq"

    @property
    def gained_pp(self) -> float:
        """На сколько процентных пунктов выросла заполняемость."""
        return self.achieved_call_rate - self.baseline_call_rate


def _median(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    n = len(values)
    mid = n // 2
    if n % 2:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2.0


def build_curve(
    skeleton_keys: Iterable[str],
    imputed_rsq: dict[str, float],
    measured_keys: set[str],
    *,
    floor: float = RSQ_FLOOR,
    step: float = RSQ_STEP,
) -> tuple[list[CurvePoint], list[float], int, int]:
    """
    Возвращает (кривая, отсортированные Rsq импутированных строк,
    закрыто измерениями, всего строк трафарета).

    skeleton_keys — ключи "<хромосома>_<позиция>" ВСЕХ строк трафарета, в
    том же виде и в том же количестве, как их считает
    template/assembler.py::assemble_final() (включая повторы позиций,
    унаследованные из самого трафарета) — иначе предсказанная
    заполняемость разойдётся с той, что потом посчитает validate_output().

    measured_keys имеют приоритет над импутацией (так делает
    merge_dictionaries), поэтому строка, закрытая измерением, заполнена
    при ЛЮБОМ пороге и в сортируемый массив Rsq не попадает.
    """
    measured_used = 0
    total_rows = 0
    rsq_values: list[float] = []
    for key in skeleton_keys:
        total_rows += 1
        if key in measured_keys:
            measured_used += 1
            continue
        rsq = imputed_rsq.get(key)
        if rsq is not None:
            rsq_values.append(rsq)
    rsq_values.sort()

    curve: list[CurvePoint] = []
    if total_rows:
        # Сетка идёт от пола вверх; порог 0.99 включительно.
        t = floor
        while t <= 0.99 + 1e-9:
            # Сколько импутированных строк переживёт порог t: все, у кого
            # Rsq >= t. Массив отсортирован, поэтому это одно бинарное
            # деление, а не проход по миллиону значений на каждый порог.
            idx = bisect.bisect_left(rsq_values, t)
            imputed_used = len(rsq_values) - idx
            call_rate = 100.0 * (measured_used + imputed_used) / total_rows
            curve.append(CurvePoint(round(t, 4), call_rate,
                                    imputed_used, measured_used))
            t += step
    return curve, rsq_values, measured_used, total_rows


def tune(
    skeleton_keys: Iterable[str],
    imputed_rsq: dict[str, float],
    measured_keys: set[str],
    *,
    target_call_rate: Optional[float] = None,
    floor: Optional[float] = None,
    baseline: Optional[float] = None,
    metric: str = "rsq",
) -> TuningResult:
    """
    Подбирает САМЫЙ ВЫСОКИЙ порог, при котором заполняемость не ниже
    целевой. "Самый высокий" — принципиально: цель в том, чтобы пройти
    приёмку, а не в том, чтобы набить файл вызовами похуже; лишнее
    занижение порога ухудшает файл без всякой пользы.

    target_call_rate=None — режим "только отчёт": порог остаётся
    baseline, но кривая и цена всё равно посчитаны и их можно показать.

    metric — в какой шкале лежат значения imputed_rsq ("rsq" или "gp").
    floor и baseline по умолчанию берутся из неё; передавать их явно
    нужно только чтобы переопределить.
    """
    metric_floor, metric_baseline, _ = settings_for(metric)
    if floor is None:
        floor = metric_floor
    if baseline is None:
        baseline = metric_baseline

    curve, rsq_values, measured_used, total_rows = build_curve(
        skeleton_keys, imputed_rsq, measured_keys, floor=floor,
    )
    if not curve or not total_rows:
        return TuningResult(
            target_call_rate=target_call_rate,
            chosen_threshold=baseline,
            achieved_call_rate=0.0,
            baseline_threshold=baseline,
            target_reached=target_call_rate is None,
            metric=metric,
        )

    def _at(threshold: float) -> CurvePoint:
        """Точка кривой, ближайшая снизу к заданному порогу."""
        best = curve[0]
        for point in curve:
            if point.threshold <= threshold + 1e-9:
                best = point
            else:
                break
        return best

    baseline_point = _at(baseline)
    max_point = curve[0]  # пол сетки даёт максимальную заполняемость

    if target_call_rate is None:
        chosen = baseline_point
        reached = True
    else:
        candidates = [p for p in curve if p.call_rate >= target_call_rate]
        if candidates:
            chosen = max(candidates, key=lambda p: p.threshold)
            reached = True
        else:
            # Цель недостижима даже на полу — берём пол и говорим об этом
            # прямо, а не делаем вид, что подобрали.
            chosen = max_point
            reached = False

    idx = bisect.bisect_left(rsq_values, chosen.threshold)
    used = rsq_values[idx:]
    below_standard = max(0, bisect.bisect_left(rsq_values, baseline) - idx)

    return TuningResult(
        target_call_rate=target_call_rate,
        chosen_threshold=chosen.threshold,
        achieved_call_rate=chosen.call_rate,
        baseline_threshold=baseline,
        baseline_call_rate=baseline_point.call_rate,
        target_reached=reached,
        max_call_rate=max_point.call_rate,
        below_standard_used=below_standard,
        median_rsq_used=_median(used),
        total_rows=total_rows,
        curve=curve,
        metric=metric,
    )


def filter_by_threshold(genotypes: dict[str, str], imputed_rsq: dict[str, float],
                        threshold: float) -> dict[str, str]:
    """
    Оставляет только те импутированные генотипы, чей Rsq не ниже порога.
    Позиции без Rsq (их в rsq_out не бывает, но на всякий случай) остаются:
    отсутствие качества трактуется как "фильтровать нечего", ровно как в
    template/assembler.py.
    """
    return {
        key: gt for key, gt in genotypes.items()
        if imputed_rsq.get(key, 1.0) >= threshold
    }


def format_result(r: TuningResult) -> str:
    metric_floor, _, label = settings_for(r.metric)
    lines = [
        "=" * 70,
        f"ПОДБОР ПОРОГА {label} ПОД ЦЕЛЕВУЮ ЗАПОЛНЯЕМОСТЬ",
        "=" * 70,
    ]
    if r.target_call_rate is None:
        lines.append(f"Цель не задана — порог оставлен стандартным "
                     f"({r.baseline_threshold:.2f}).")
    elif r.target_reached:
        lines.append(f"Цель: не ниже {r.target_call_rate:.2f}% — достигнута.")
    else:
        lines.append(
            f"Цель: не ниже {r.target_call_rate:.2f}% — НЕ достигнута. "
            f"Максимум, что дают эти дозы при пороге {metric_floor:.2f}, — "
            f"{r.max_call_rate:.2f}%. Ниже {metric_floor:.2f} порог не "
            f"опускается сознательно: там уже не генотипы, а шум."
        )
    lines += [
        "",
        f"Выбранный порог {label}:  {r.chosen_threshold:.2f}",
        f"Заполняемость:            {r.achieved_call_rate:.2f}%",
        f"Для сравнения, при {r.baseline_threshold:.2f}:  "
        f"{r.baseline_call_rate:.2f}%  "
        f"({r.gained_pp:+.2f} п.п.)",
        "",
        "Чем заплачено:",
        f"  позиций с {label} ниже {r.baseline_threshold:.2f} в итоговом файле: "
        f"{r.below_standard_used:,}".replace(",", " "),
        f"  медиана {label} использованных импутированных позиций: "
        f"{r.median_rsq_used:.3f}",
        "=" * 70,
    ]
    return "\n".join(lines)


def metrics_row(r: Optional[TuningResult]) -> dict:
    if r is None:
        return {}
    return {
        # ⚠ Столбцы оставлены с прежними именами, чтобы не рвать историю
        # runs_metrics.csv. Но 0,30 по Rsq и 0,90 по GP — разные шкалы,
        # поэтому рядом обязательно пишется, чем именно подбирали.
        "quality_metric": r.metric,
        "rsq_target_call_rate": (round(r.target_call_rate, 2)
                                 if r.target_call_rate is not None else ""),
        "rsq_chosen": round(r.chosen_threshold, 3),
        "rsq_target_reached": int(bool(r.target_reached)),
        "call_rate_at_standard_rsq": round(r.baseline_call_rate, 4),
        "call_rate_max_possible": round(r.max_call_rate, 4),
        "positions_below_standard_rsq": r.below_standard_used,
        "median_rsq_used": round(r.median_rsq_used, 4),
    }
