"""
core/network_utils.py
Единая точка правды для настройки сетевого окружения, от которого зависит
удалённый доступ bcftools/libcurl к HTTPS-источникам (1000 Genomes S3/EBI) —
используется _probe_bcftools_remote_support()/process_chromosome_remote() в
download_donors.py.

=============================================================================
ПРЕДЫСТОРИЯ (промт "Диагностика + устойчивая настройка CA-сертификатов"):
=============================================================================

При ручной диагностике на Windows-машине пользователя `bcftools view -h
"https://1000genomes.s3.amazonaws.com/..."` падал с SSL-ошибкой — сборка
bcftools/libcurl на Windows не использует системное хранилище сертификатов
Windows и ищет файл `cacert.pem` явно (через переменные окружения
CURL_CA_BUNDLE/SSL_CERT_FILE). Это чинится вручную:

    curl.exe -o bin\\cacert.pem https://curl.se/ca/cacert.pem
    $env:CURL_CA_BUNDLE = "...\\bin\\cacert.pem"
    $env:SSL_CERT_FILE  = "...\\bin\\cacert.pem"

...но это переменные окружения ТЕКУЩЕЙ консоли — при следующем запуске из
нового окна PowerShell, из ярлыка или из GUI-процесса они не установлены, и
проблема возвращается. Раньше приложение молча полагалось на то, что
пользователь сам настроил окружение перед запуском — теперь это делает сам
процесс, один раз при старте (см. ensure_ca_bundle()/ensure_network_ready()).

Второй найденный при той же диагностике баг: в bin/ (папка бинарников
htslib, указываемая через --bin-dir) лежит СВОЙ curl.exe (в комплекте со
сборкой bcftools/htslib), и когда main.py/gui/app.py добавляют bin_dir в
начало PATH (чтобы найти bgzip/tabix/bcftools по имени в некоторых местах),
это заодно подменяет системный curl.exe из System32 на бандловый — а у
бандлового зашит relative-путь к сертификатам вида
"<корень установки>/etc/ssl/certs/ca-bundle.crt", которого на диске
пользователя нет. Результат — SSL-ошибка даже после того, как
CURL_CA_BUNDLE выставлен для bcftools, потому что ломается совсем другой
процесс (сам curl.exe, вызываемый из download_donors.py как fallback
полного скачивания).

Итоговая цепочка диагностики была: SSL-ошибка -> нашли и поставили
cacert.pem -> SSL-ошибка сохранилась -> нашли конфликт bin/curl.exe ->
переименовали в curl.exe.disabled -> SSL пропала, но всплыла "No such file
or directory" (URL с суффиксом v5b, которого нет на S3-зеркале, — не
связано с сетью, обычная работа фолбэка суффиксов) -> с правильным
суффиксом v5a запрос успешно вернул VCF-заголовок.

Модуль ниже автоматизирует обе части фикса (сертификаты + конфликт
curl.exe), чтобы пользователю больше не нужно было ничего делать вручную ни
в одной новой сессии/ярлыке/GUI-запуске.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

CA_BUNDLE_FILENAME = "cacert.pem"
CA_BUNDLE_URL = "https://curl.se/ca/cacert.pem"
# Разумный минимум размера актуального набора сертификатов curl.se —
# защита от того, что скачался обрезанный/битый файл (HTML страница с
# ошибкой вместо .pem, редирект на страницу-заглушку и т.п.). Актуальный
# файл на момент написания — около 200+ КБ, берём заведомо заниженный
# порог, чтобы не ловить ложные срабатывания на будущих чуть меньших
# версиях набора, но отсечь совсем испорченные скачивания.
CA_BUNDLE_MIN_SIZE = 50 * 1024  # 50 КБ

# Имена curl, которые ищем в bin_dir как потенциально конфликтующие.
_CURL_BIN_NAMES = ("curl.exe", "curl")


# ---------------------------------------------------------------------------
# Часть 2.1 — CA-сертификаты для libcurl
# ---------------------------------------------------------------------------
def _download_ca_bundle(dest: Path) -> None:
    """
    Скачивает актуальный набор CA-сертификатов с curl.se с докачкой
    (по аналогии с main.py::_download_with_resume()) — используется
    urllib, а не системный curl, чтобы не зависеть от того самого
    curl.exe, чью сертификатную проблему мы, возможно, ещё не починили
    (курица и яйцо: нельзя чинить сертификаты curl.exe через сам curl.exe,
    если у него именно эта проблема).
    """
    import urllib.request
    import urllib.error

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_dest = dest.with_suffix(dest.suffix + ".part")
    existing_size = tmp_dest.stat().st_size if tmp_dest.exists() else 0

    headers = {"User-Agent": "Mozilla/5.0"}
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"

    req = urllib.request.Request(CA_BUNDLE_URL, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            resumed = getattr(response, "status", 200) == 206
            mode = "ab" if (resumed and existing_size > 0) else "wb"
            if not resumed:
                existing_size = 0
            with open(tmp_dest, mode) as f:
                while True:
                    chunk = response.read(65536)
                    if not chunk:
                        break
                    f.write(chunk)
    except urllib.error.HTTPError as e:
        if e.code == 416:
            pass  # уже скачано полностью
        else:
            tmp_dest.unlink(missing_ok=True)
            raise RuntimeError(
                f"Не удалось скачать {CA_BUNDLE_URL}: HTTP {e.code} {e.reason}"
            ) from e
    except Exception as e:
        tmp_dest.unlink(missing_ok=True)
        raise RuntimeError(f"Не удалось скачать {CA_BUNDLE_URL}: {e}") from e

    if not tmp_dest.exists() or tmp_dest.stat().st_size < CA_BUNDLE_MIN_SIZE:
        size = tmp_dest.stat().st_size if tmp_dest.exists() else 0
        tmp_dest.unlink(missing_ok=True)
        raise RuntimeError(
            f"Скачанный {CA_BUNDLE_FILENAME} подозрительно мал ({size} байт, "
            f"ожидалось >= {CA_BUNDLE_MIN_SIZE} байт) — вероятно, скачалась "
            f"страница с ошибкой вместо файла сертификатов."
        )
    tmp_dest.replace(dest)


def ensure_ca_bundle(bin_dir: Optional[Path]) -> Optional[Path]:
    """
    Гарантирует наличие набора CA-сертификатов (bin_dir/cacert.pem) и
    выставляет CURL_CA_BUNDLE/SSL_CERT_FILE В ТЕКУЩЕМ ПРОЦЕССЕ Python
    (os.environ) — так это работает одинаково при запуске из любого
    терминала, ярлыка или GUI, без требования к пользователю вручную
    выставлять переменные окружения перед каждым запуском.

    Если bin_dir не задан — сертификаты негде хранить постоянно
    (используем рабочую директорию как fallback) и функция всё равно
    пытается их обеспечить, чтобы удалённый доступ имел шанс сработать.

    Если скачивание не удалось (нет интернета, сеть недоступна) — НЕ
    бросает исключение наружу: пишет предупреждение в лог и возвращает
    None. Вызывающий код (ensure_network_ready) должен относиться к этому
    как к сигналу "удалённый доступ, скорее всего, не заработает" — а не
    как к фатальной ошибке: _probe_bcftools_remote_support() и так сама
    проверяет реальную работоспособность и откатывается на полное
    скачивание при неудаче.

    Устанавливает CURL_CA_BUNDLE/SSL_CERT_FILE даже если сертификаты уже
    были прописаны иначе на уровне системы/пользователя (переопределяет
    их на заведомо рабочий файл для этого процесса) — это самый надёжный
    вариант: он не полагается на то, что случайно настроено в окружении
    снаружи, и не может конфликтовать с чужими системными путями.
    """
    target_dir = Path(bin_dir) if bin_dir else Path.cwd()
    ca_path = target_dir / CA_BUNDLE_FILENAME

    if ca_path.exists() and ca_path.stat().st_size >= CA_BUNDLE_MIN_SIZE:
        logger.info("✓ CA-сертификаты найдены: %s", ca_path)
    else:
        if ca_path.exists():
            logger.warning(
                "⚠ Файл %s повреждён или неполон — перекачиваю", ca_path
            )
            ca_path.unlink(missing_ok=True)
        logger.info(
            "ℹ CA-сертификаты для libcurl не найдены — скачиваю с %s "
            "(нужно для удалённого доступа bcftools к HTTPS без ручной "
            "настройки $env:CURL_CA_BUNDLE)", CA_BUNDLE_URL,
        )
        try:
            _download_ca_bundle(ca_path)
        except RuntimeError as e:
            logger.warning(
                "⚠ Не удалось получить CA-сертификаты (%s) — удалённый "
                "доступ bcftools к HTTPS, скорее всего, не заработает "
                "(SSL-ошибка), но пайплайн продолжит работу через полное "
                "скачивание хромосом (fallback).", e,
            )
            return None
        logger.info("✓ CA-сертификаты скачаны: %s", ca_path)

    os.environ["CURL_CA_BUNDLE"] = str(ca_path)
    os.environ["SSL_CERT_FILE"] = str(ca_path)
    return ca_path


# ---------------------------------------------------------------------------
# Часть 2.2 — конфликтующий curl.exe в bin_dir
# ---------------------------------------------------------------------------
def find_conflicting_bin_curl(bin_dir: Optional[Path]) -> Optional[Path]:
    """
    Возвращает путь к curl(.exe) внутри bin_dir, если он там есть, иначе
    None. На Linux/macOS bin/curl обычно не бандлится вместе со сборкой
    bcftools/htslib (в отличие от Windows-дистрибутивов) — там эта функция
    в норме просто ничего не найдёт, no-op.
    """
    if not bin_dir:
        return None
    bin_dir = Path(bin_dir)
    for name in _CURL_BIN_NAMES:
        candidate = bin_dir / name
        if candidate.is_file():
            return candidate
    return None


def which_curl_ignoring_dir(exclude_dir: Optional[Path]) -> Optional[str]:
    """
    Аналог shutil.which("curl"), но временно убирающий exclude_dir (обычно
    --bin-dir) из PATH на время поиска — гарантирует, что найдётся
    системный curl (System32 на Windows, /usr/bin на Linux/macOS), а не
    потенциально битый бандловый curl.exe из папки бинарников htslib (см.
    докстринг модуля: относительный путь к сертификатам внутри бандлового
    curl.exe ломает SSL даже после того, как CURL_CA_BUNDLE выставлен для
    bcftools — это два РАЗНЫХ процесса).

    Используется вместо голого shutil.which("curl") везде, где
    download_donors.py ищет системный curl как способ полного скачивания
    (download_with_curl / _downloader_chain).
    """
    if not exclude_dir:
        return shutil.which("curl")

    exclude_str = str(Path(exclude_dir).resolve())
    original_path = os.environ.get("PATH", "")
    filtered_entries = [
        p for p in original_path.split(os.pathsep)
        if p and Path(p).resolve() != Path(exclude_str)
    ]
    filtered_path = os.pathsep.join(filtered_entries)

    old_path = os.environ.get("PATH")
    try:
        os.environ["PATH"] = filtered_path
        return shutil.which("curl")
    finally:
        if old_path is not None:
            os.environ["PATH"] = old_path


def warn_if_conflicting_curl(bin_dir: Optional[Path]) -> None:
    """
    Логирует явное предупреждение, если в bin_dir лежит свой curl(.exe) —
    не переименовывает и не удаляет файл автоматически (это было бы
    неожиданным побочным эффектом для пользователя, трогающим файлы вне
    рабочих директорий приложения), только предупреждает и объясняет,
    что для скачивания донорских файлов и cacert.pem приложение всё
    равно будет использовать системный curl в обход этого файла (см.
    which_curl_ignoring_dir()).
    """
    conflicting = find_conflicting_bin_curl(bin_dir)
    if conflicting is None:
        return
    logger.warning(
        "⚠ В папке бинарников (%s) обнаружен собственный %s — сборки "
        "curl, идущие в комплекте с htslib/bcftools на Windows, часто "
        "используют относительный путь к сертификатам и падают с SSL-"
        "ошибкой, даже если CURL_CA_BUNDLE настроен верно для bcftools "
        "(это два разных процесса). Приложение будет ИГНОРИРОВАТЬ этот "
        "curl при поиске системного curl для скачивания (см. "
        "which_curl_ignoring_dir) — переименовывать/удалять файл не "
        "требуется, но при ручных экспериментах в консоли имейте в виду, "
        "что 'curl' в PATH может резолвиться в этот, потенциально "
        "нерабочий, экземпляр, если бинарники добавлены в начало PATH.",
        conflicting.parent, conflicting.name,
    )


# ---------------------------------------------------------------------------
# Единая точка инициализации (Часть 2.3)
# ---------------------------------------------------------------------------
def ensure_network_ready(bin_dir: Optional[Path]) -> None:
    """
    Единая точка входа: вызывается один раз при старте (CLI main(),
    GUI при инициализации HTSLIB, CLI download_donors.py) — до любого
    использования bcftools с удалённым URL (перед
    _probe_bcftools_remote_support()).

    Делает оба фикса из докстринга модуля:
      1. ensure_ca_bundle(bin_dir)  — сертификаты для libcurl (bcftools).
      2. warn_if_conflicting_curl(bin_dir) — предупреждение о конфликте
         curl.exe (сам обход конфликта — на уровне
         which_curl_ignoring_dir(), используемого в download_donors.py
         вместо голого shutil.which("curl")).

    Никогда не бросает исключение — при любой проблеме (нет интернета,
    странная файловая система) пишет предупреждение в лог и продолжает:
    это подготовка окружения для УСКОРЕНИЯ (Часть 1 предыдущего промта),
    а не критический путь — при неудаче пайплайн просто откатывается на
    полное скачивание хромосом, как и раньше.
    """
    try:
        ensure_ca_bundle(bin_dir)
    except Exception as e:  # noqa: BLE001 — подготовка окружения не должна ронять пайплайн
        logger.warning("⚠ ensure_ca_bundle завершилась с ошибкой (%s) — продолжаю без неё.", e)

    try:
        warn_if_conflicting_curl(bin_dir)
    except Exception as e:  # noqa: BLE001
        logger.warning("⚠ warn_if_conflicting_curl завершилась с ошибкой (%s) — продолжаю без неё.", e)
