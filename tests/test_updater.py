"""
tests/test_updater.py

Проверка обновлений тихая по замыслу: при любой проблеме она возвращает
None и приложение просто запускается. Обратная сторона этого — сломаться
она может незаметно: пользователь никогда не увидит предложения
обновиться и никогда не пожалуется, потому что ничего и не показывалось.
Поэтому логика сравнения версий и разбора ответа GitHub закрыта тестами,
а не проверяется руками.
"""
from __future__ import annotations

import json

import pytest

from core import updater


# ---------------------------------------------------------------------------
# Сравнение версий
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("candidate,current,expected", [
    ("1.3.0", "1.2.5", True),
    ("v1.3.0", "1.2.5", True),        # тег с ведущей "v"
    ("1.2.5", "1.2.5", False),        # та же самая версия — не новее
    ("1.2.4", "1.2.5", False),
    ("1.2.10", "1.2.9", True),        # не строковое сравнение: 10 > 9
    ("2.0.0", "1.99.99", True),
    ("1.3", "1.2.5", True),           # укороченная запись
    ("1.3.0", "1.3.0-beta1", True),   # финал новее своего предрелиза
    ("1.3.0-beta1", "1.3.0", False),  # и обратно — не новее
    ("1.3.0-beta2", "1.3.0-beta1", True),
    ("latest", "1.2.5", False),       # неразбираемый тег -> не предлагаем
    ("", "1.2.5", False),
    ("1.3.0", "непонятно", False),    # неизвестна точка отсчёта
])
def test_is_newer(candidate, current, expected):
    assert updater.is_newer(candidate, current) is expected


def test_parse_version_prerelease():
    assert updater.parse_version("v1.3.0-rc.2") == ((1, 3, 0), ("rc", 2))


# ---------------------------------------------------------------------------
# Разбор ответа GitHub
# ---------------------------------------------------------------------------
def _payload(**over):
    base = {
        "tag_name": "v1.3.0",
        "name": "HelixFillDNA 1.3.0",
        "body": "Что изменилось...",
        "html_url": "https://github.com/o/r/releases/tag/v1.3.0",
        "draft": False,
        "prerelease": False,
        "assets": [
            {"name": "notes.txt", "browser_download_url": "u1", "size": 10},
            {"name": "HelixFillDNA-Setup-1.3.0.exe",
             "browser_download_url": "u2", "size": 52428800},
        ],
    }
    base.update(over)
    return base


def test_parse_release_json_picks_exe_asset():
    info = updater.parse_release_json(_payload())
    assert info is not None
    assert info.version == "1.3.0"
    assert info.tag == "v1.3.0"
    assert info.asset_name == "HelixFillDNA-Setup-1.3.0.exe"
    assert info.asset_url == "u2"
    assert round(info.size_mb) == 50


def test_parse_release_json_rejects_draft_and_prerelease():
    assert updater.parse_release_json(_payload(draft=True)) is None
    assert updater.parse_release_json(_payload(prerelease=True)) is None


def test_parse_release_json_rejects_unparseable_tag():
    assert updater.parse_release_json(_payload(tag_name="nightly")) is None
    assert updater.parse_release_json(_payload(tag_name="")) is None


def test_parse_release_json_survives_missing_assets():
    info = updater.parse_release_json(_payload(assets=[]))
    assert info is not None and info.asset_url == ""


# ---------------------------------------------------------------------------
# Настройки пользователя в ui_state.json
# ---------------------------------------------------------------------------
def test_settings_default_enabled(tmp_path):
    state = tmp_path / "ui_state.json"
    settings = updater.load_update_settings(state)
    assert settings["enabled"] is True
    assert settings["skip_version"] == ""


def test_save_settings_keeps_other_keys(tmp_path):
    """
    ui_state.json общий с настройкой режима интерфейса — раздел обновлений
    не должен затирать settings_mode, иначе пользователь после первой же
    проверки обновлений получал бы сброшенный режим.
    """
    state = tmp_path / "ui_state.json"
    state.write_text(json.dumps({"settings_mode": "Продвинутые настройки"}),
                     encoding="utf-8")
    updater.disable_update_checks(state)
    data = json.loads(state.read_text(encoding="utf-8"))
    assert data["settings_mode"] == "Продвинутые настройки"
    assert data["updates"]["enabled"] is False


def test_save_settings_survives_broken_json(tmp_path):
    state = tmp_path / "ui_state.json"
    state.write_text("{это не json", encoding="utf-8")
    updater.skip_version(state, "1.3.0")
    assert updater.load_update_settings(state)["skip_version"] == "1.3.0"


# ---------------------------------------------------------------------------
# check_for_update: решения пользователя
# ---------------------------------------------------------------------------
def _stub_fetch(monkeypatch, version="1.3.0"):
    info = updater.parse_release_json(_payload(tag_name=f"v{version}"))
    monkeypatch.setattr(updater, "fetch_latest_release", lambda timeout=0: info)
    return info


def test_check_offers_newer(monkeypatch, tmp_path):
    _stub_fetch(monkeypatch)
    got = updater.check_for_update("1.2.5", tmp_path / "s.json")
    assert got is not None and got.version == "1.3.0"


def test_check_silent_when_same_version(monkeypatch, tmp_path):
    _stub_fetch(monkeypatch, "1.2.5")
    assert updater.check_for_update("1.2.5", tmp_path / "s.json") is None


def test_check_respects_skip_version(monkeypatch, tmp_path):
    state = tmp_path / "s.json"
    _stub_fetch(monkeypatch)
    updater.skip_version(state, "1.3.0")
    assert updater.check_for_update("1.2.5", state) is None
    # ...но следующая версия всё равно предлагается
    _stub_fetch(monkeypatch, "1.4.0")
    assert updater.check_for_update("1.2.5", state) is not None


def test_check_respects_disabled(monkeypatch, tmp_path):
    state = tmp_path / "s.json"
    _stub_fetch(monkeypatch)
    updater.disable_update_checks(state)
    assert updater.check_for_update("1.2.5", state) is None


def test_force_ignores_disabled_and_skip(monkeypatch, tmp_path):
    """
    Кнопка "Проверить обновления" — явный вопрос пользователя, и на него
    надо ответить, даже если он когда-то нажал "не напоминать". При этом
    сама настройка не меняется.
    """
    state = tmp_path / "s.json"
    _stub_fetch(monkeypatch)
    updater.disable_update_checks(state)
    updater.skip_version(state, "1.3.0")
    assert updater.check_for_update("1.2.5", state, force=True) is not None
    assert updater.load_update_settings(state)["enabled"] is False


def test_check_never_raises(monkeypatch, tmp_path):
    """Падение сети не должно мешать запуску приложения."""
    def boom(timeout=0):
        raise updater.UpdateCheckError("нет сети")
    monkeypatch.setattr(updater, "fetch_latest_release", boom)
    assert updater.check_for_update("1.2.5", tmp_path / "s.json") is None

    def worse(timeout=0):
        raise RuntimeError("что-то совсем неожиданное")
    monkeypatch.setattr(updater, "fetch_latest_release", worse)
    assert updater.check_for_update("1.2.5", tmp_path / "s.json") is None
