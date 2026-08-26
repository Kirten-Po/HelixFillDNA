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
        cmd = [seven_zip, "x", str(zip_path), f"-o{dest_dir}", f"-p{pwd_str}", "-y"]
        result = subprocess.run(cmd, capture_output=True, text=True, shell=False)
        if result.returncode != 0:
            raise ArchiveExtractionError(
                f"7z завершился с кодом {result.returncode} при распаковке "
                f"{zip_path.name}:\n{result.stderr.strip()}"
            )
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
    except (zipfile.BadZipFile, RuntimeError, NotImplementedError) as e:
        raise ArchiveExtractionError(
            f"Ошибка при распаковке {zip_path.name}: {e}. Скорее всего "
            f"архив зашифрован AES-256 — установите 7-Zip и повторите."
        ) from e


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

    failed: list[tuple[Path, str]] = []
    for zip_path in zip_paths:
        try:
            extract_zip(zip_path, dest_dir, password, sevenzip_path)
        except ArchiveExtractionError as e:
            sanitized = sanitize_password_text(password)
            if sanitized != password:
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
