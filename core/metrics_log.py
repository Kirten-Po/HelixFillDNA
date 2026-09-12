"""
core/metrics_log.py
Одна строка на запуск в runs_metrics.csv — то, из чего потом можно
откалибровать пороги, которые сейчас приходится угадывать.

Мотив прямой. Пороги пре-флайта (полоса гетерозиготности), порог редкого
хвоста для выбора панели и сам порог приёмки заказчика мы сейчас можем
только зажать между двумя известными точками: 88 % (отклонили) и 91,3 %
(приняли). По двум файлам границу не вывести. Поэтому каждый запуск
дописывает сюда полный набор: метрики исходника, панель, порог Rsq,
финальный call rate — и колонку `verdict`, которую пользователь
заполняет РУКАМИ, когда получает ответ приёмки. После десятка запусков
граница нарисуется сама.

Формат — CSV, а не JSON и не база: файл нужно открывать в Excel, строить
по нему точечную диаграмму и дописывать вердикт в последнюю колонку.
Ничего из этого JSON не умеет.

Колонки не фиксированы навечно: если в будущей версии появится новая
метрика, файл со старой шапкой автоматически переписывается с
объединённой шапкой (см. append_run()) — старые строки сохраняются, в
новых колонках у них пусто. Терять историю из-за добавленной колонки
недопустимо: ради накопления истории всё и затевалось.
"""
from __future__ import annotations

import csv
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

METRICS_FILENAME = "runs_metrics.csv"

#: Порядок колонок по умолчанию — от "что это был за запуск" к "что
#: получилось" и дальше к ручному вердикту, чтобы файл читался слева
#: направо как история одного прогона.
BASE_COLUMNS: tuple[str, ...] = (
    "timestamp",
    "run_name",
    "app_version",
    "source",
    "panel",
    "panel_recommended",
    "format",
    # --- исходник (core/preflight.py) ---
    "markers",
    "call_rate",
    "autosomal_het_pct",
    "palindromic_het_pct",
    "x_nonpar_calls",
    "x_het_pct",
    "y_markers",
    "y_call_rate",
    "mt_markers",
    "sex_by_x",
    "sex_by_y",
    "duplicate_positions",
    "duplicate_rsids",
    "unsorted_positions",
    "malformed_rows",
    "indel_calls",
    "build",
    "build_match_pct",
    "preflight_level",
    # --- состав чипа (core/panel_advisor.py) ---
    "chip_positions",
    "chip_matched_1000g",
    "chip_missing_1000g_pct",
    "chip_rare_eur_pct",
    "chip_monomorphic_eur_pct",
    # --- парсинг (ParseResult) ---
    "parsed_variants",
    "both_non_ref_pct",
    "het_self_complementary_pct",
    "lift_failed",
    # --- подбор порога (core/rsq_tuner.py) ---
    "rsq_target_call_rate",
    "rsq_chosen",
    "rsq_target_reached",
    "call_rate_at_standard_rsq",
    "call_rate_max_possible",
    "positions_below_standard_rsq",
    "median_rsq_used",
    # --- результат ---
    "final_call_rate",
    "final_file",
    # --- заполняется РУКАМИ после ответа приёмки ---
    "verdict",
    "notes",
)


def metrics_path(root: Path) -> Path:
    return Path(root) / METRICS_FILENAME


def _read_existing(path: Path) -> tuple[list[str], list[dict]]:
    if not path.is_file():
        return [], []
    try:
        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            header = list(reader.fieldnames or [])
            rows = [dict(r) for r in reader]
        return header, rows
    except (OSError, csv.Error) as e:
        logger.warning("Не удалось прочитать %s (%s) — файл будет создан заново", path, e)
        return [], []


def append_run(root: Path, row: dict) -> Optional[Path]:
    """
    Дописывает строку. Возвращает путь к файлу или None, если записать не
    удалось (запуск из read-only папки и т.п.).

    Накопление истории — вспомогательная задача, поэтому НИ ОДНА ошибка
    здесь не должна ронять пайплайн: всё гасится и логируется.
    """
    try:
        path = metrics_path(root)
        path.parent.mkdir(parents=True, exist_ok=True)
        clean = {k: ("" if v is None else v) for k, v in row.items()}
        clean.setdefault("timestamp", datetime.now().isoformat(timespec="seconds"))

        header, rows = _read_existing(path)
        if not header:
            # Новый файл: полная шапка, даже если часть колонок в этом
            # запуске пуста — так в Excel сразу видно, что вообще
            # собирается, и куда вписывать вердикт.
            header = list(BASE_COLUMNS)
            for key in clean:
                if key not in header:
                    header.append(key)
            with open(path, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
                writer.writeheader()
                writer.writerow(clean)
            return path

        new_keys = [k for k in clean if k not in header]
        if new_keys:
            # Шапка изменилась (обновилась версия приложения) — переписываем
            # файл целиком с объединённой шапкой, сохраняя все старые строки.
            header = header + new_keys
            rows.append(clean)
            with open(path, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
                writer.writeheader()
                for r in rows:
                    writer.writerow(r)
            return path

        with open(path, "a", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
            writer.writerow(clean)
        return path
    except Exception as e:  # noqa: BLE001 — лог метрик не критичен
        logger.warning("Не удалось записать метрики запуска: %s", e)
        return None


def parse_result_row(result) -> dict:
    """Метрики парсинга из ParseResult — то, чего нет у пре-флайта."""
    if result is None:
        return {}
    try:
        return {
            "parsed_variants": len(result.variants),
            "both_non_ref_pct": round(result.both_non_ref_pct, 4),
            "het_self_complementary_pct": round(
                result.het_self_complementary_pct, 4
            ),
            "lift_failed": getattr(result, "lift_failed", 0),
        }
    except Exception:  # noqa: BLE001
        return {}
