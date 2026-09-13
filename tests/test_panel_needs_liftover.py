"""
Признак «панель требует переноса координат» не должен зависеть от того,
какая панель выбрана по умолчанию.

Разбор реального отказа (прогон ATLAS-TOPMED). В пяти местах кода стояло
условие вида `panel != DEFAULT_PANEL`, и читалось оно как «панель не HRC,
значит нужен лифтовер» — потому что DEFAULT_PANEL был "hrc". В 1.3.2
DEFAULT_PANEL стал "topmed", и все пять условий молча инвертировались:

  * для TopMed лифтовер перестал выполняться (Этап 7 не переносил ни
    panel_pos вперёд, ни результат обратно);
  * для HRC он, наоборот, стал запрашиваться там, где не нужен.

Цена: скелет трафарета остаётся в GRCh37, а дозы и измерения — в GRCh38.
Позиции одного локуса в разных сборках почти никогда не совпадают числом,
поэтому в итоговый файл попадали только случайные совпадения — 1,15 %
заполнения вместо ~98 %. Ни одной ошибки в лог при этом не попало:
словарь генотипов собирался правильно (573 534 записи), просто его ключи
были из другой системы координат.

Поэтому тест проверяет не текущее значение, а ИНВАРИАНТ: ответ зависит
только от сборки генома панели и не меняется при подмене DEFAULT_PANEL.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import main  # noqa: E402


def test_hrc_needs_no_liftover():
    """HRC уже в GRCh37 — переносить нечего."""
    assert main._panel_needs_liftover("hrc") is False


def test_topmed_needs_liftover():
    """TopMed в GRCh38 — без переноса результат не совпадёт со скелетом."""
    assert main._panel_needs_liftover("topmed") is True


@pytest.mark.parametrize("default_panel", ["hrc", "topmed"])
def test_answer_does_not_depend_on_default_panel(monkeypatch, default_panel):
    """
    ГЛАВНАЯ проверка: смена панели по умолчанию не должна менять ответ
    ни для одной панели. Ровно это сломалось в 1.3.2.
    """
    monkeypatch.setattr(main, "DEFAULT_PANEL", default_panel)
    assert main._panel_needs_liftover("hrc") is False
    assert main._panel_needs_liftover("topmed") is True


def test_every_panel_agrees_with_its_genome_build():
    """
    Инвариант на весь реестр панелей: лифтовер нужен тогда и только тогда,
    когда сборка панели отличается от сборки скелета. Новая панель,
    добавленная в REFERENCE_PANELS, попадает под проверку автоматически.
    """
    for name, cfg in main.REFERENCE_PANELS.items():
        expected = cfg["genome_build"] != main.SKELETON_BUILD
        assert main._panel_needs_liftover(name) is expected, (
            f"панель {name!r} объявлена в сборке {cfg['genome_build']!r}, "
            f"а _panel_needs_liftover() вернул {main._panel_needs_liftover(name)}"
        )


def test_panel_needing_liftover_has_chain_urls():
    """
    Если панель требует переноса, у неё должны быть chain-файлы в обе
    стороны: вперёд (чип -> сборка панели, Этап 1) и обратно (результат
    -> GRCh37, Этап 7). Панель без обратного chain дала бы ту же тихую
    порчу, что и пропущенное условие.
    """
    for name, cfg in main.REFERENCE_PANELS.items():
        if not main._panel_needs_liftover(name):
            continue
        assert cfg.get("liftover_chain_url") or cfg.get("liftover_chain_urls"), (
            f"панель {name!r} в чужой сборке, но chain-файл не задан"
        )


def test_default_panel_is_a_known_panel():
    """Мелочь, но дешёвая: опечатка в DEFAULT_PANEL ломает весь пайплайн."""
    assert main.DEFAULT_PANEL in main.REFERENCE_PANELS
