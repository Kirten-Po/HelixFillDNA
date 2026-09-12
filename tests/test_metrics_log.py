"""
tests/test_metrics_log.py

runs_metrics.csv существует ради одной вещи — накопить историю, по
которой можно будет откалибровать пороги вместо того, чтобы зажимать их
между двумя известными точками (88 % отклонили, 91,3 % приняли).
Поэтому главное проверяемое свойство — история НЕ ТЕРЯЕТСЯ: ни при
добавлении новой колонки в будущей версии, ни при битом файле, ни при
ошибке записи.
"""
from __future__ import annotations

import csv

from core import metrics_log as ml


def _rows(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def test_creates_file_with_full_header(tmp_path):
    ml.append_run(tmp_path, {"run_name": "a", "final_call_rate": 91.3})
    path = ml.metrics_path(tmp_path)
    rows = _rows(path)
    assert len(rows) == 1
    assert rows[0]["run_name"] == "a"
    # Колонка вердикта должна быть в файле сразу — её заполняют руками,
    # и человек должен видеть, куда писать.
    assert "verdict" in rows[0]
    assert rows[0]["timestamp"]


def test_appends_without_rewriting(tmp_path):
    ml.append_run(tmp_path, {"run_name": "a"})
    ml.append_run(tmp_path, {"run_name": "b"})
    rows = _rows(ml.metrics_path(tmp_path))
    assert [r["run_name"] for r in rows] == ["a", "b"]


def test_new_column_preserves_old_rows(tmp_path):
    """
    Появилась новая метрика — старые строки обязаны остаться на месте, с
    пустым значением в новой колонке. Потерять историю из-за добавленной
    колонки недопустимо: ради истории всё и затевалось.
    """
    ml.append_run(tmp_path, {"run_name": "старый"})
    ml.append_run(tmp_path, {"run_name": "новый", "метрика_из_будущего": 42})
    rows = _rows(ml.metrics_path(tmp_path))
    assert len(rows) == 2
    assert rows[0]["run_name"] == "старый"
    assert rows[0]["метрика_из_будущего"] == ""
    assert rows[1]["метрика_из_будущего"] == "42"


def test_manual_verdict_survives_next_run(tmp_path):
    """Вердикт приёмки вписывается руками — следующий запуск его не трёт."""
    ml.append_run(tmp_path, {"run_name": "a"})
    path = ml.metrics_path(tmp_path)
    rows = _rows(path)
    header = list(rows[0].keys())
    rows[0]["verdict"] = "принято"
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)

    ml.append_run(tmp_path, {"run_name": "b", "метрика_из_будущего": 1})
    rows = _rows(path)
    assert rows[0]["verdict"] == "принято"


def test_write_failure_is_not_fatal(tmp_path, monkeypatch):
    """Лог метрик вспомогательный — он не имеет права уронить пайплайн."""
    def boom(*a, **kw):
        raise OSError("диск только для чтения")
    monkeypatch.setattr("builtins.open", boom)
    assert ml.append_run(tmp_path, {"run_name": "a"}) is None


def test_broken_existing_file_is_recreated(tmp_path):
    path = ml.metrics_path(tmp_path)
    path.write_bytes(b"\x00\x01\x02 not a csv")
    assert ml.append_run(tmp_path, {"run_name": "a"}) is not None
    assert _rows(path)[0]["run_name"] == "a"


def test_parse_result_row_survives_odd_object():
    class Weird:
        pass
    assert ml.parse_result_row(Weird()) == {}
    assert ml.parse_result_row(None) == {}
