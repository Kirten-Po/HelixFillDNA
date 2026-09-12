"""
core/updater.py
Проверка обновлений на GitHub при запуске приложения.

Зачем отдельный модуль, а не пара строк в gui/app.py:
  * логика сравнения версий и разбора ответа GitHub тестируется без GUI
    (tests/test_updater.py) — иначе её пришлось бы проверять руками, а
    ошибка здесь тихая: приложение просто никогда не предложит обновиться;
  * тем же кодом пользуется CLI (main.py --check-updates), а не только окно.

Договорённости, которые важно не потерять:
  * проверка НИКОГДА не мешает работе. Любая ошибка (нет сети, GitHub
    ответил 403 из-за лимита анонимных запросов, ответ не разобрался)
    гасится и возвращается None — приложение просто запускается молча;
  * запрос идёт в фоновом потоке с коротким таймаутом, чтобы окно не
    зависало на старте на медленной сети;
  * решение пользователя "не напоминать" сохраняется в ui_state.json
    рядом с exe — там же, где живёт settings_mode, а не в реестре/AppData.
    Различаются ДВА решения: "пропустить конкретную версию" (skip_version)
    и "не проверять вообще" (enabled=False);
  * приложение НИЧЕГО не скачивает и не устанавливает само — только
    открывает страницу релиза в браузере. Автоподмена собственного exe на
    Windows требует перезапуска через промежуточный процесс, прав на
    запись в Program Files и подписи; для инструмента, который и так
    ставится установщиком, это лишний класс проблем.
"""
from __future__ import annotations

import contextlib
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Репозиторий, где публикуются релизы. Держится здесь, а не в gui/app.py,
# чтобы CLI и GUI смотрели в одно и то же место.
GITHUB_OWNER = "Kirten-Po"
GITHUB_REPO = "HelixFillDNA"
RELEASES_API_URL = (
    f"https://api.github.com/repos/{GITHUB_OWNER}/{GITHUB_REPO}/releases/latest"
)
RELEASES_PAGE_URL = f"https://github.com/{GITHUB_OWNER}/{GITHUB_REPO}/releases/latest"

# Короткий таймаут: проверка обновлений — украшение, а не работа. Лучше
# молча не проверить, чем задержать старт.
DEFAULT_TIMEOUT = 6.0

# Ключ в ui_state.json, под которым живут настройки обновлений.
STATE_KEY = "updates"


class UpdateCheckError(RuntimeError):
    """Проверка не удалась. Наружу из check_for_update() не выходит."""


# ---------------------------------------------------------------------------
# Сравнение версий
# ---------------------------------------------------------------------------
_VERSION_RE = re.compile(
    r"^\s*v?(\d+)(?:\.(\d+))?(?:\.(\d+))?(?:[-+.]?([0-9A-Za-z.\-]+))?\s*$"
)


def parse_version(text: str) -> Optional[tuple[tuple[int, int, int], tuple]]:
    """
    "v1.3.0" -> ((1, 3, 0), ()), "1.3.0-beta2" -> ((1, 3, 0), ("beta", 2)).

    Возвращает None, если строка вообще не похожа на версию — тогда
    сравнивать нечего и обновление не предлагается (лучше не предложить,
    чем предложить "обновиться" на тег вида "latest").

    Предрелизы считаются СТАРШЕ ничего и МЛАДШЕ финального релиза той же
    тройки (semver): 1.3.0-beta2 < 1.3.0. Это нужно, чтобы человек,
    поставивший бету, получил предложение перейти на финал, а не наоборот.
    """
    if not text:
        return None
    m = _VERSION_RE.match(str(text))
    if not m:
        return None
    major, minor, patch, pre = m.groups()
    core = (int(major), int(minor or 0), int(patch or 0))
    if not pre:
        return core, ()
    parts: list = []
    for chunk in re.split(r"[.\-]", pre):
        if not chunk:
            continue
        parts.append(int(chunk) if chunk.isdigit() else chunk)
    return core, tuple(parts)


def _pre_key(pre: tuple) -> tuple:
    """
    Ключ сортировки предрелизной части. Пустая часть (финальный релиз)
    должна быть БОЛЬШЕ любой непустой — отсюда ведущий флаг 1/0.
    Внутри: числа сравниваются как числа и считаются младше строк
    (semver), поэтому каждый элемент разворачивается в пару.
    """
    if not pre:
        return (1,)
    return (0,) + tuple((0, p, "") if isinstance(p, int) else (1, 0, p) for p in pre)


def is_newer(candidate: str, current: str) -> bool:
    """
    Строго ли candidate новее current. Неразбираемая candidate — False
    (не предлагаем обновление на непонятный тег); неразбираемая current —
    тоже False (не знаем, от чего считать).
    """
    a = parse_version(candidate)
    b = parse_version(current)
    if a is None or b is None:
        return False
    return (a[0], _pre_key(a[1])) > (b[0], _pre_key(b[1]))


# ---------------------------------------------------------------------------
# Ответ GitHub
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ReleaseInfo:
    version: str          # "1.3.0" — тег без ведущей "v"
    tag: str              # "v1.3.0" — как он есть на GitHub
    name: str             # заголовок релиза (может совпадать с тегом)
    notes: str            # тело релиза (markdown), может быть пустым
    page_url: str         # страница релиза для браузера
    asset_url: str = ""   # прямая ссылка на установщик, если он приложен
    asset_name: str = ""
    asset_size: int = 0

    @property
    def size_mb(self) -> float:
        return self.asset_size / 1024 ** 2 if self.asset_size else 0.0


def parse_release_json(payload: dict) -> Optional[ReleaseInfo]:
    """
    Разбирает ответ GitHub Releases API. Отдельная функция — чтобы её
    можно было проверить тестом на сохранённом JSON, без сети.

    Черновики и предрелизы, помеченные на самом GitHub (draft/prerelease),
    отбрасываются: /releases/latest их и так не отдаёт, но эта же функция
    используется и для разбора списка релизов.
    """
    if not isinstance(payload, dict):
        return None
    if payload.get("draft") or payload.get("prerelease"):
        return None
    tag = str(payload.get("tag_name") or "").strip()
    if not tag:
        return None
    parsed = parse_version(tag)
    if parsed is None:
        return None
    core, pre = parsed
    version = ".".join(str(x) for x in core)
    if pre:
        version += "-" + ".".join(str(x) for x in pre)

    asset_url = asset_name = ""
    asset_size = 0
    for asset in payload.get("assets") or []:
        name = str(asset.get("name") or "")
        if name.lower().endswith(".exe"):
            asset_url = str(asset.get("browser_download_url") or "")
            asset_name = name
            with contextlib.suppress(TypeError, ValueError):
                asset_size = int(asset.get("size") or 0)
            break

    return ReleaseInfo(
        version=version,
        tag=tag,
        name=str(payload.get("name") or tag),
        notes=str(payload.get("body") or "").strip(),
        page_url=str(payload.get("html_url") or RELEASES_PAGE_URL),
        asset_url=asset_url,
        asset_name=asset_name,
        asset_size=asset_size,
    )


def fetch_latest_release(timeout: float = DEFAULT_TIMEOUT) -> ReleaseInfo:
    """
    Спрашивает у GitHub последний релиз. Бросает UpdateCheckError при
    любой проблеме — вызывающий check_for_update() её гасит.

    context= обязателен: без него на Windows с битым хранилищем корневых
    сертификатов urlopen падает ещё до подключения (см.
    core/network_utils.py::make_ssl_context()).
    """
    import urllib.error
    import urllib.request

    from core.network_utils import make_ssl_context

    req = urllib.request.Request(
        RELEASES_API_URL,
        headers={
            # GitHub отвечает 403 на запросы без User-Agent.
            "User-Agent": f"{GITHUB_REPO}-updater",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout,
                                    context=make_ssl_context()) as response:
            raw = response.read(1024 * 1024)
    except urllib.error.HTTPError as e:
        # 403 — обычный лимит анонимных запросов к API (60/час на IP),
        # 404 — в репозитории ещё нет ни одного релиза. И то и другое не
        # повод шуметь в интерфейсе.
        raise UpdateCheckError(f"GitHub ответил HTTP {e.code}") from e
    except Exception as e:  # noqa: BLE001 — сеть/SSL/что угодно
        raise UpdateCheckError(str(e)) from e

    try:
        payload = json.loads(raw.decode("utf-8", errors="replace"))
    except ValueError as e:
        raise UpdateCheckError("ответ GitHub не разобрался как JSON") from e

    info = parse_release_json(payload)
    if info is None:
        raise UpdateCheckError("в ответе GitHub нет пригодного тега релиза")
    return info


# ---------------------------------------------------------------------------
# Настройки пользователя (ui_state.json)
# ---------------------------------------------------------------------------
def _read_state(state_file: Path) -> dict:
    try:
        data = json.loads(Path(state_file).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def load_update_settings(state_file: Path) -> dict:
    """
    {"enabled": bool, "skip_version": str} — с умолчаниями.
    Проверка включена по умолчанию: пользователь, который её не трогал,
    должен узнавать о новых версиях.
    """
    section = _read_state(state_file).get(STATE_KEY)
    if not isinstance(section, dict):
        section = {}
    return {
        "enabled": bool(section.get("enabled", True)),
        "skip_version": str(section.get("skip_version") or ""),
        "last_seen": str(section.get("last_seen") or ""),
    }


def save_update_settings(state_file: Path, **fields) -> None:
    """
    Дописывает переданные ключи в раздел updates, не трогая остальной
    ui_state.json (там же лежит settings_mode). Ошибки записи гасятся:
    запуск из read-only папки не должен ронять приложение.
    """
    state_file = Path(state_file)
    data = _read_state(state_file)
    section = data.get(STATE_KEY)
    if not isinstance(section, dict):
        section = {}
    section.update(fields)
    data[STATE_KEY] = section
    try:
        state_file.parent.mkdir(parents=True, exist_ok=True)
        state_file.write_text(
            json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    except OSError as e:
        logger.debug("Не удалось сохранить настройки обновлений: %s", e)


def disable_update_checks(state_file: Path) -> None:
    """Кнопка «Больше не напоминать» — выключает проверку целиком."""
    save_update_settings(state_file, enabled=False)


def skip_version(state_file: Path, version: str) -> None:
    """«Пропустить эту версию» — про неё больше не спрашиваем, про следующую спросим."""
    save_update_settings(state_file, skip_version=str(version))


# ---------------------------------------------------------------------------
# Главная точка входа
# ---------------------------------------------------------------------------
def check_for_update(
    current_version: str,
    state_file: Path,
    *,
    force: bool = False,
    timeout: float = DEFAULT_TIMEOUT,
) -> Optional[ReleaseInfo]:
    """
    Возвращает ReleaseInfo, если ЕСТЬ что предложить пользователю, иначе
    None. Никогда не бросает исключений — вызывается из фонового потока
    при старте, где падение некому поймать.

    force=True (кнопка "Проверить обновления" в настройках) игнорирует и
    выключенную проверку, и пропущенную версию: человек спросил явно —
    надо ответить, даже если он когда-то нажал "не напоминать". При этом
    сама настройка не меняется.
    """
    settings = load_update_settings(state_file)
    if not force and not settings["enabled"]:
        return None

    try:
        info = fetch_latest_release(timeout=timeout)
    except UpdateCheckError as e:
        logger.info("Проверка обновлений не удалась (%s) — продолжаю без неё.", e)
        return None
    except Exception as e:  # noqa: BLE001 — страховка: старт важнее проверки
        logger.info("Проверка обновлений не удалась (%s) — продолжаю без неё.", e)
        return None

    if not is_newer(info.version, current_version):
        return None
    if not force and settings["skip_version"] == info.version:
        return None
    return info
