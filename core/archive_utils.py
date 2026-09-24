"""
core/archive_utils.py
Единый модуль распаковки ZIP-архивов результатов Michigan Imputation Server.

Раньше логика распаковки (поиск 7z.exe, санитайзинг пароля, вызов
subprocess) была продублирована в main.py (_extract_zip) и в
mis_adapter.py (MISAdapter._extract_one). Теперь оба места вызывают
функции отсюда — единая точка правды (Задача 6).

Почему пароль передаётся как bytes, а не str (Задача 4):
Пароли MIS иногда содержат спецсимволы '(', ')', '&'. На Windows
subprocess.run с списком аргументов (shell=False, как у нас) сам по себе
не подвержен классическому shell-инъекция/экранированию `&`/`()` —
CreateProcess получает аргументы напрямую. Тем не менее 7-Zip читает
пароль как CLI-аргумент в кодировке консоли, и при копировании пароля
из HTML-письма в него нередко попadают невидимые unicode-символы
(NBSP, zero-width space, BOM) — вот это и ломает ввод. Поэтому пароль
сначала прогоняется через _sanitize_password(), а сравнение/подстановка
делается на уровне bytes (UTF-8), чтобы не терять контроль над тем,
какие байты реально уходят во внешний процесс.
"""
from __future__ import annotations

import logging
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Пробелы + типичные невидимые unicode-символы, которые часто попадают
# в буфер обмена при копировании пароля из HTML-письма.
_INVISIBLE_PASSWORD_CHARS = re.compile(r"[\u00a0\u200b\u200c\u200d\ufeff\s]")


class ArchiveExtractionError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Пароль
# ---------------------------------------------------------------------------
def sanitize_password_text(raw: str) -> str:
    """Убирает пробелы и невидимые unicode-символы, возвращает строку."""
    return _INVISIBLE_PASSWORD_CHARS.sub("", raw)


def sanitize_password_bytes(raw: str) -> bytes:
    """
    Нормализует пароль и возвращает его как bytes (UTF-8).
    Используется там, где важно контролировать точные байты, уходящие
    во внешний процесс (7z) или в zipfile.extractall(pwd=...).
    """
    cleaned = sanitize_password_text(raw)
    return cleaned.encode("utf-8")


# ---------------------------------------------------------------------------
# Поиск 7z.exe
# ---------------------------------------------------------------------------
def find_7z(explicit_path: Optional[str] = None) -> Optional[str]:
    if explicit_path:
        p = Path(explicit_path)
        if p.is_file():
            return str(p)
    found = shutil.which("7z") or shutil.which("7z.exe")
    if found:
        return found
    for p in (
        Path(r"C:\Program Files\7-Zip\7z.exe"),
        Path(r"C:\Program Files (x86)\7-Zip\7z.exe"),
    ):
        if p.exists():
            return str(p)
    return None


# ---------------------------------------------------------------------------
# Распаковка одного архива
# ---------------------------------------------------------------------------
def extract_zip(
    zip_path: Path,
    dest_dir: Path,
    password: str,
    sevenzip_path: Optional[str] = None,
) -> None:
    """
    Распаковывает один ZIP-архив в dest_dir с указанным паролем.

    Приоритет:
      1. 7-Zip (обязателен для AES-256, который используют архивы MIS —
         встроенный zipfile его не поддерживает и либо падает, либо в
         некоторых версиях Python молча создаёт битые файлы даже при
         верном пароле).
      2. Встроенный zipfile — только запасной путь, с явным
         предупреждением о возможных проблемах с AES-256.

    Пароль передаётся в 7z как отдельный байтовый аргумент (encode/decode
    через UTF-8 без потерь), чтобы спецсимволы вроде '(', ')', '&' не
    портились при передаче.
    """
    zip_path = Path(zip_path)
    dest_dir = Path(dest_dir)
    if not zip_path.exists():
        raise ArchiveExtractionError(f"ZIP файл не найден: {zip_path}")
    dest_dir.mkdir(parents=True, exist_ok=True)

    pwd_bytes = sanitize_password_bytes(password) if password else b""
    # Если очистка ничего не изменила, всё равно используем очищенную
    # версию — она либо совпадает с исходной, либо чинит проблему.
    pwd_str = pwd_bytes.decode("utf-8")

    seven_zip = find_7z(sevenzip_path)
    if seven_zip:
        # -p<пароль> собирается как один аргумент списка — subprocess
        # с shell=False передаёт его в CreateProcess/execve без участия
        # shell-парсера, так что скобки/амперсанды не интерпретируются.
        #
        # -sccUTF-8: 7-Zip пишет в консоль в OEM-кодировке (cp866 на
        # русской Windows), а text=True декодировал вывод как cp1251 — и
        # пользователь видел вместо «Недостаточно места на диске» набор
        # «KI a®6B BB3-®». Просим у 7z UTF-8 и так же его и читаем.
        cmd = [seven_zip, "x", str(zip_path), f"-o{dest_dir}", f"-p{pwd_str}",
               "-y", "-sccUTF-8"]
        result = subprocess.run(
            cmd, capture_output=True, text=True, encoding="utf-8",
            errors="replace", shell=False,
        )
        if result.returncode != 0:
            details = (result.stderr or "").strip() or (result.stdout or "").strip()
            message = (
                f"7z завершился с кодом {result.returncode} при распаковке "
                f"{zip_path.name}:\n{details}"
            )
            if _is_disk_full(details):
                message += "\n\n" + _disk_full_hint(dest_dir)
            raise ArchiveExtractionError(message)
        return

    logger.warning(
        "7z.exe не найден — распаковываю через встроенный zipfile. "
        "Архивы MIS используют AES-256, который zipfile не поддерживает: "
        "если распаковка упадёт или файлы окажутся битыми, установите "
        "7-Zip (https://7-zip.org) и укажите путь к нему."
    )
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(path=dest_dir, pwd=pwd_bytes)
    except OSError as e:
        if _is_disk_full(str(e)) or getattr(e, "errno", None) == 28:
            raise ArchiveExtractionError(
                f"Ошибка при распаковке {zip_path.name}: {e}\n\n"
                + _disk_full_hint(dest_dir)
            ) from e
        raise
    except (zipfile.BadZipFile, RuntimeError, NotImplementedError) as e:
        raise ArchiveExtractionError(
            f"Ошибка при распаковке {zip_path.name}: {e}. Скорее всего "
            f"архив зашифрован AES-256 — установите 7-Zip и повторите."
        ) from e


# ---------------------------------------------------------------------------
# Свободное место на диске
# ---------------------------------------------------------------------------
#: Как 7-Zip и ОС сообщают о нехватке места: английский вывод 7z, русский
#: текст Windows-ошибки ERROR_DISK_FULL (7z подставляет его как есть) и
#: POSIX ENOSPC.
_DISK_FULL_MARKERS = (
    "not enough space", "недостаточно места", "no space left",
    "disk full", "диск заполнен",
)

#: Запас сверх объёма распакованных файлов: tabix-индексы, временные
#: файлы Этапа 7 и сам итоговый файл.
DISK_SPACE_MARGIN = 2 * 1024 ** 3


def _is_disk_full(text: str) -> bool:
    low = (text or "").lower()
    return any(marker in low for marker in _DISK_FULL_MARKERS)


def _gb(n: int) -> str:
    return f"{n / 1024 ** 3:.1f} ГБ"


def _disk_full_hint(dest_dir: Path) -> str:
    try:
        free = _gb(shutil.disk_usage(dest_dir).free)
    except OSError:
        free = "не удалось определить"
    return (
        f"Причина — на диске кончилось место (свободно сейчас: {free}). "
        f"Результаты импутации распаковываются в {dest_dir}; распакованные "
        f"дозы занимают десятки гигабайт. Освободите место на этом диске "
        f"(или перенесите папку программы на диск побольше) и повторите "
        f"Этап 7 — уже скачанные архивы заново не скачиваются.\n"
        f"Файл хромосомы, на которой оборвалась распаковка, может остаться "
        f"недописанным: повторная распаковка перезапишет его."
    )


def required_extract_space(zip_paths: list[Path], dest_dir: Path) -> int:
    """
    Сколько байт ДОБАВИТСЯ на диске после распаковки архивов.

    Размеры берутся из оглавления ZIP — оно не шифруется даже у архивов
    с AES-256, поэтому пароль не нужен. Файлы, которые уже лежат в
    dest_dir под тем же именем, будут перезаписаны: их текущий размер
    вычитается. Нечитаемый архив пропускается — о нём внятно скажет
    сама распаковка.
    """
    dest_dir = Path(dest_dir)
    need = 0
    for zip_path in zip_paths:
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        continue
                    target = dest_dir / info.filename
                    existing = target.stat().st_size if target.is_file() else 0
                    need += max(0, info.file_size - existing)
        except (OSError, zipfile.BadZipFile):
            continue
    return need


def check_free_space(zip_paths: list[Path], dest_dir: Path) -> None:
    """
    Проверяет ДО распаковки, что архивы поместятся на диск.

    Без этой проверки 7-Zip успевал распаковать часть хромосом и падал
    посреди следующей с «Cannot set length for output file», а на диске
    оставался недописанный chrN.dose.vcf.gz.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    need = required_extract_space(zip_paths, dest_dir)
    try:
        free = shutil.disk_usage(dest_dir).free
    except OSError:
        return
    if need + DISK_SPACE_MARGIN > free:
        raise ArchiveExtractionError(
            f"Недостаточно места на диске для распаковки результатов: нужно "
            f"около {_gb(need + DISK_SPACE_MARGIN)} (из них {_gb(need)} — сами "
            f"файлы), свободно {_gb(free)}.\n"
            f"Папка распаковки: {dest_dir}\n"
            f"Освободите место на этом диске (или перенесите папку программы "
            f"на диск побольше) и повторите Этап 7 — уже скачанные архивы "
            f"заново не скачиваются."
        )


# ---------------------------------------------------------------------------
# Распаковка списка архивов (с автоматическим повтором на "очищенном" пароле)
# ---------------------------------------------------------------------------
def extract_all(
    zip_paths: list[Path],
    dest_dir: Path,
    password: str,
    sevenzip_path: Optional[str] = None,
) -> None:
    """
    Распаковывает все переданные архивы. Если распаковка с исходным
    паролем не удалась и его "очистка" (sanitize) реально что-то поменяла,
    делает вторую попытку с очищенным паролем.
    """
    if not zip_paths:
        raise ArchiveExtractionError("Список ZIP-архивов пуст")

    check_free_space(zip_paths, dest_dir)

    failed: list[tuple[Path, str]] = []
    for zip_path in zip_paths:
        try:
            extract_zip(zip_path, dest_dir, password, sevenzip_path)
        except ArchiveExtractionError as e:
            sanitized = sanitize_password_text(password)
            # Нехватка места — не про пароль: повтор с «очищенным» паролем
            # упал бы так же, только позже и с запутывающим сообщением.
            if sanitized != password and not _is_disk_full(str(e)):
                logger.warning(
                    "Распаковка %s с исходным паролем не удалась (%s). "
                    "Пробую с очищенным от пробелов/невидимых символов паролем.",
                    zip_path.name, e,
                )
                try:
                    extract_zip(zip_path, dest_dir, sanitized, sevenzip_path)
                    continue
                except ArchiveExtractionError as e2:
                    failed.append((zip_path, str(e2)))
                    continue
            failed.append((zip_path, str(e)))

    if failed:
        details = "\n".join(f"  - {p.name}: {err}" for p, err in failed)
        raise ArchiveExtractionError(
            f"Не удалось распаковать {len(failed)} из {len(zip_paths)} архивов:\n{details}"
        )
    logger.info("Все %d архивов распакованы в: %s", len(zip_paths), dest_dir)


# ---------------------------------------------------------------------------
# Диагностика пароля (Задача 7) — используется кнопкой "🔍 Диагностика" в GUI
# ---------------------------------------------------------------------------
def diagnose_password(
    password: str,
    test_archive: Path,
    sevenzip_path: Optional[str] = None,
) -> list[str]:
    """
    Возвращает список диагностических строк для показа пользователю:
      - длина пароля / наличие невидимых символов
      - наличие 7z.exe
      - пробная распаковка bin/_password_test.7z (если он существует)
    Ничего не бросает — все проблемы отражены в тексте сообщений.
    """
    lines: list[str] = [f"Длина пароля: {len(password)} символов"]

    if password:
        sanitized = sanitize_password_text(password)
        if sanitized != password:
            lines.append(
                f"⚠ Обнаружены пробелы/невидимые символы! "
                f"После очистки длина: {len(sanitized)}"
            )
        else:
            lines.append("✓ Невидимых символов не обнаружено")
        lines.append(f"Первые 3 символа: {password[:3]}***")
    else:
        lines.append("⚠ Поле пароля пустое")

    seven_zip = find_7z(sevenzip_path)
    lines.append(
        f"✓ 7z.exe найден: {seven_zip}" if seven_zip
        else "❌ 7z.exe НЕ найден! Установите: https://7-zip.org"
    )

    test_archive = Path(test_archive)
    if test_archive.exists() and seven_zip and password:
        cmd = [seven_zip, "t", str(test_archive), f"-p{sanitize_password_text(password)}", "-y"]
        res = subprocess.run(cmd, capture_output=True, text=True, shell=False)
        if res.returncode == 0:
            lines.append("✓ Тестовый архив распакован успешно — пароль рабочий")
        else:
            lines.append(
                f"❌ Тестовый архив НЕ распакован (код {res.returncode}). "
                f"stderr: {res.stderr.strip()[:300]}"
            )
    else:
        lines.append("ℹ Тестовый архив bin/_password_test.7z не найден — проверка пропущена")

    return lines
