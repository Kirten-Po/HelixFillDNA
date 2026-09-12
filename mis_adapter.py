"""
mis_adapter.py
Адаптер для работы с Michigan Imputation Server.
Подготовка 22 VCF-файлов для загрузки, открытие браузера,
скачивание и распаковка результатов.
"""
from __future__ import annotations
import atexit
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import threading
import time
import webbrowser
import zipfile
from pathlib import Path
from typing import Callable, Optional

from core.pure_python_core import UPLOAD_CHROMS
from core.archive_utils import (
    extract_zip, extract_all, find_7z, ArchiveExtractionError,
)

logger = logging.getLogger(__name__)


# ===========================================================================
# Скачивание архивов результатов MIS
#
# Разбор реального отказа (жалоба "chr_16 качается по 11 КБ/с, то
# появляется, то пропадает; потом [WinError 32] файл занят другим
# процессом"):
#
#   1. curl запускался как `curl -sL <url> -o <финальное имя>`:
#      * `-s` глушит сообщения об ошибках, поэтому до пользователя
#        долетало только "returned non-zero exit status 56" — без единого
#        слова о том, что 56 означает обрыв приёма данных;
#      * без `-f` HTTP-ошибка (404, 500, страница логина) не считается
#        ошибкой: curl пишет тело ответа в файл и возвращает 0;
#      * без `--speed-limit` соединение, просевшее до 11 КБ/с, живёт
#        вечно — 250 МБ на такой скорости это больше шести часов, и со
#        стороны это выглядит как "зависло";
#      * без `-C -` каждый обрыв начинал файл заново, с нуля;
#      * запись сразу в финальное имя означает, что оборванная закачка
#        остаётся на диске под именем настоящего архива.
#   2. Дочерний curl.exe не убивался при закрытии окна: поток GUI —
#      daemon, Python при выходе его снимает, но ДОЧЕРНИЙ ПРОЦЕСС на
#      Windows продолжает жить и писать в chr_16.zip. Отсюда и "две
#      хромосомы качаются одновременно" (скачивание здесь строго
#      последовательное — второй писатель был осиротевшим curl от
#      прошлого запуска), и [WinError 32] на следующем запуске: файл
#      занят живым процессом, а `dest.unlink()` не был защищён и ронял
#      весь Шаг 3.
#   3. Файлы, уже лежащие в папке, пропускались по признаку "валидный
#      ZIP" — без всякой привязки к ЗАДАНИЮ. Перезапустив задание на MIS
#      (новый job id, новый пароль), пользователь получал старые архивы,
#      которые честно проходили проверку целостности, пропускались как
#      "уже скачанные" — и падали на распаковке с "Wrong password" по
#      всем 23 сразу.
#
# Ниже — фиксы всех трёх пунктов.
# ===========================================================================

#: Расшифровка кодов возврата curl. Голое "exit status 56" не говорит
#: пользователю ничего, а причина отказа у этих кодов разная настолько,
#: что и действия разные: 56 — чинить связь, 22 — ссылка устарела, 60 —
#: сертификаты.
_CURL_EXIT_HINTS: dict[int, str] = {
    6: "не удалось определить адрес сервера (DNS)",
    7: "не удалось подключиться к серверу",
    18: "передача оборвалась, файл получен не целиком",
    22: "сервер вернул HTTP-ошибку — обычно это истёкшая ссылка "
        "(результаты на MIS хранятся ~3 дня) или неверный адрес",
    23: "не удалось записать файл на диск (нет места или нет прав)",
    28: "истёк таймаут: соединение слишком долго молчало или скорость "
        "надолго упала почти до нуля",
    33: "сервер не поддерживает докачку с середины файла",
    35: "ошибка установки защищённого соединения (TLS)",
    52: "сервер ответил пустым ответом",
    55: "не удалось отправить данные в сеть",
    56: "обрыв приёма данных: соединение разорвано на середине закачки",
    60: "сертификат сервера не проверен — обычно это антивирус с "
        "перехватом HTTPS или повреждённое хранилище сертификатов Windows",
}


def curl_error_text(code: int, stderr: str = "") -> str:
    """Человекочитаемое описание отказа curl вместо голого кода."""
    hint = _CURL_EXIT_HINTS.get(code)
    tail = (stderr or "").strip().splitlines()
    detail = tail[-1].strip() if tail else ""
    parts = [f"curl завершился с кодом {code}"]
    if hint:
        parts.append(hint)
    if detail:
        parts.append(detail)
    return ": ".join(parts[:2]) + (f" ({detail})" if detail and hint else "")


# ---------------------------------------------------------------------------
# Реестр запущенных curl-процессов
#
# Без него дочерний curl.exe переживает закрытие окна и продолжает писать
# в файл результатов — именно так возникает [WinError 32] "файл занят
# другим процессом" на СЛЕДУЮЩЕМ запуске и вторая "качающаяся" строка в
# панели прогресса. subprocess.run() тут не годится вовсе: он не отдаёт
# наружу объект процесса, а значит убить его некому.
# ---------------------------------------------------------------------------
_ACTIVE_CURLS: set[subprocess.Popen] = set()
_ACTIVE_CURLS_LOCK = threading.Lock()


def _register_curl(proc: subprocess.Popen) -> None:
    with _ACTIVE_CURLS_LOCK:
        _ACTIVE_CURLS.add(proc)


def _unregister_curl(proc: subprocess.Popen) -> None:
    with _ACTIVE_CURLS_LOCK:
        _ACTIVE_CURLS.discard(proc)


def kill_active_curls() -> int:
    """
    Убивает все запущенные этим процессом закачки. Вызывается на выходе
    (atexit) и при отмене пользователем. Возвращает число убитых.
    """
    with _ACTIVE_CURLS_LOCK:
        procs = list(_ACTIVE_CURLS)
    killed = 0
    for proc in procs:
        if proc.poll() is not None:
            continue
        try:
            proc.kill()
            killed += 1
        except OSError:
            pass
    return killed


atexit.register(kill_active_curls)


class DownloadCancelled(RuntimeError):
    """Пользователь отменил скачивание результатов."""


class FileLockedError(RuntimeError):
    """Файл держит другой процесс — почти всегда осиротевший curl.exe."""


def _locked_file_message(path: Path, err: Exception) -> str:
    return (
        f"Файл {path.name} занят другим процессом и его нельзя ни удалить, "
        f"ни перезаписать ({err}).\n\n"
        f"Почти наверняка это curl.exe, оставшийся от предыдущего запуска: "
        f"он продолжает медленно докачивать этот же архив, даже если окно "
        f"программы было закрыто.\n\n"
        f"Что сделать: откройте Диспетчер задач (Ctrl+Shift+Esc), вкладка "
        f"«Подробности», снимите все процессы curl.exe — и нажмите "
        f"«Скачать результаты» ещё раз. Эта версия программы больше не "
        f"оставляет таких процессов после себя."
    )


def _safe_unlink(path: Path) -> None:
    """
    Удаление с внятной ошибкой вместо голого [WinError 32].

    Раньше `dest.unlink(missing_ok=True)` на занятом файле выбрасывал
    PermissionError, который никто не ловил, и весь Шаг 3 падал с
    сообщением "[WinError 32] Процесс не может получить доступ к файлу" —
    без единого намёка на то, ЧТО за процесс и что с этим делать.
    """
    try:
        path.unlink(missing_ok=True)
    except PermissionError as e:
        raise FileLockedError(_locked_file_message(path, e)) from e
    except OSError as e:
        raise MISAdapterError(f"Не удалось удалить {path}: {e}") from e


# ---------------------------------------------------------------------------
# Манифест скачивания: какой архив от какого задания
# ---------------------------------------------------------------------------
DOWNLOAD_MANIFEST_NAME = "download_manifest.json"

#: Сколько раз подряд пытаться скачать ОДИН файл, прежде чем признать
#: неудачу. Раньше цикл повторов был `while True` без счётчика: на
#: нестабильном канале пользователь мог до бесконечности жать «Да» в
#: диалоге «повторить скачивание этого файла?», каждый раз начиная файл с
#: нуля. Теперь попытки считаются, а докачка идёт с места обрыва.
MAX_FILE_ATTEMPTS = 3


def job_id_from_urls(urls: list[str]) -> str:
    """
    Отпечаток задания MIS — по набору ссылок на архивы. Ссылки содержат
    идентификатор задания, поэтому два разных запуска на сервере дают
    разный отпечаток, а повторный запуск программы для ТОГО ЖЕ задания —
    тот же самый.
    """
    h = hashlib.sha1()
    for url in sorted(urls):
        h.update(url.encode("utf-8", errors="replace"))
        h.update(b"\n")
    return h.hexdigest()[:16]


def read_manifest(results_dir: Path) -> dict:
    try:
        data = json.loads((Path(results_dir) / DOWNLOAD_MANIFEST_NAME)
                          .read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def write_manifest(results_dir: Path, manifest: dict) -> None:
    try:
        (Path(results_dir) / DOWNLOAD_MANIFEST_NAME).write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8",
        )
    except OSError as e:
        logger.debug("Не удалось записать манифест скачивания: %s", e)


# ---------------------------------------------------------------------------
# Один файл: скачивание с докачкой, отменой и защитой от залипания
# ---------------------------------------------------------------------------
#: Ниже этой скорости (байт/с), продержавшейся _STALL_SECONDS, соединение
#: считается залипшим и обрывается — чтобы сработала докачка с другого
#: соединения, а не тянулись сутки. Живой пример: chr_16.zip шёл на
#: 11 КБ/с; 250 МБ на такой скорости — больше шести часов, и со стороны
#: это неотличимо от зависшей программы.
_MIN_SPEED_BYTES = 20 * 1024
_STALL_SECONDS = 60
#: Пауза между опросами процесса — на неё же завязана реакция на отмену.
_POLL_SECONDS = 0.5


def _curl_base_args() -> list[str]:
    return [
        "curl",
        # -s глушит прогресс-бар, -S ВОЗВРАЩАЕТ сообщения об ошибках
        # (раньше стоял голый -s, и до пользователя долетал только код
        # возврата), -f делает HTTP-ошибку настоящей ошибкой, а не
        # молча скачанной страницей с текстом ошибки вместо архива.
        "-sSfL",
        "--connect-timeout", "30",
        "--speed-limit", str(_MIN_SPEED_BYTES),
        "--speed-time", str(_STALL_SECONDS),
        "--retry", "3",
        "--retry-delay", "3",
        "--retry-connrefused",
    ]


def _run_curl(args: list[str], cancel_check: Optional[Callable[[], bool]] = None,
              ) -> tuple[int, str]:
    """
    Запускает curl как Popen (а не subprocess.run) по двум причинам:
    процесс надо уметь УБИТЬ по отмене и на выходе из приложения, и надо
    успевать реагировать на отмену, пока идёт многоминутная закачка.

    Возвращает (код возврата, stderr).
    """
    try:
        proc = subprocess.Popen(
            # stdout в DEVNULL: с "-o файл" curl в него ничего не пишет, а
            # неопустошаемый PIPE — это потенциальный дедлок на большом
            # выводе. stderr нужен: в нём текст ошибки, ради которого и
            # добавлен флаг -S.
            args, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
            text=True, errors="replace",
        )
    except FileNotFoundError as e:
        raise MISAdapterError(
            "curl не найден в системе. На Windows 10+ он встроен — "
            "проверьте PATH или обновите Windows."
        ) from e

    _register_curl(proc)
    try:
        while True:
            try:
                proc.wait(timeout=_POLL_SECONDS)
                break
            except subprocess.TimeoutExpired:
                pass
            if cancel_check is not None and cancel_check():
                proc.kill()
                proc.wait(timeout=10)
                raise DownloadCancelled("Скачивание отменено пользователем")
        stderr = proc.stderr.read() if proc.stderr else ""
        return proc.returncode, stderr
    finally:
        _unregister_curl(proc)
        if proc.stderr is not None:
            try:
                proc.stderr.close()
            except OSError:
                pass


def remote_size(url: str) -> Optional[int]:
    """
    Размер файла на сервере (HEAD-запрос) или None, если сервер не
    ответил/не сообщил длину. Нужен, чтобы отличить полностью скачанный
    архив от оборванного, когда манифеста ещё нет (файлы остались от
    предыдущей версии программы).
    """
    try:
        proc = subprocess.run(
            [*_curl_base_args(), "-I", url],
            capture_output=True, text=True, errors="replace", timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    size: Optional[int] = None
    for line in proc.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            try:
                size = int(line.split(":", 1)[1].strip())
            except ValueError:
                continue
    return size


def download_one(url: str, dest: Path,
                 cancel_check: Optional[Callable[[], bool]] = None) -> None:
    """
    Скачивает один архив. Бросает MISAdapterError с человекочитаемым
    текстом при неудаче, DownloadCancelled при отмене, FileLockedError
    если файл держит чужой процесс.

    Пишет в <имя>.part и переименовывает только после проверки
    целостности. Раньше curl писал сразу в финальное имя, и оборванная
    закачка оставалась на диске под именем настоящего архива.

    Докачка (-C -) включена: обрыв на 210-м мегабайте из 250 не должен
    означать «начать заново». Если сервер докачку не поддерживает (код
    33), .part удаляется и файл качается с нуля — один раз, без
    зацикливания.
    """
    dest = Path(dest)
    part = dest.with_name(dest.name + ".part")

    for resume in (True, False):
        args = [*_curl_base_args()]
        if resume and part.exists() and part.stat().st_size > 0:
            args += ["-C", "-"]
            logger.info("Докачиваю %s с %.1f МБ", dest.name,
                        part.stat().st_size / 1024 ** 2)
        args += [url, "-o", str(part)]

        try:
            code, stderr = _run_curl(args, cancel_check)
        except PermissionError as e:
            raise FileLockedError(_locked_file_message(part, e)) from e

        if code == 0:
            break
        if code == 33 and resume:
            # Сервер не умеет отдавать кусок с середины — начинаем заново.
            logger.info("Сервер не поддерживает докачку %s — качаю заново",
                        dest.name)
            _safe_unlink(part)
            continue
        raise MISAdapterError(
            f"Не удалось скачать {dest.name}: {curl_error_text(code, stderr)}"
        )

    if not part.exists() or part.stat().st_size == 0:
        _safe_unlink(part)
        raise MISAdapterError(
            f"Файл {dest.name} скачался пустым или не скачался вовсе."
        )
    if not _is_valid_zip(part):
        size_mb = part.stat().st_size / 1024 ** 2
        _safe_unlink(part)
        raise MISAdapterError(
            f"Файл {dest.name} скачался ({size_mb:.1f} МБ), но не является "
            f"целым ZIP-архивом — закачка оборвалась либо вместо архива "
            f"пришла страница с ошибкой."
        )

    # os.replace атомарен и перезаписывает существующий файл — но на
    # Windows падает с WinError 32, если цель держит чужой процесс.
    try:
        os.replace(part, dest)
    except PermissionError as e:
        raise FileLockedError(_locked_file_message(dest, e)) from e
    except OSError as e:
        raise MISAdapterError(f"Не удалось сохранить {dest.name}: {e}") from e


# ---------------------------------------------------------------------------
# Промт "проверять уже скачанные файлы на соответствие ссылке / не
# сломанные": лёгкая проверка целостности ZIP-архива без полной
# распаковки — используется и для файлов, уже найденных на диске (перед
# тем как их пропустить), и сразу после свежего скачивания (перед тем
# как считать файл успешно скачанным).
# ---------------------------------------------------------------------------
def _is_valid_zip(path: Path) -> bool:
    """
    Проверяет, что файл — структурно корректный ZIP-архив: читает
    центральную директорию (zipfile.ZipFile(...).namelist()) и убеждается,
    что в архиве есть хотя бы одна запись.

    Это НЕ полная проверка CRC каждого файла внутри архива (для этого
    пришлось бы распаковать архив целиком, что дорого и не нужно на
    данном этапе — реальная целостность содержимого всё равно
    перепроверяется на распаковке, core/archive_utils.py) — а быстрая
    структурная проверка: битый файл (оборванная докачка, HTML-страница
    с ошибкой вместо архива, повреждение на диске) не пройдёт даже
    чтение центральной директории и здесь будет надёжно отловлен.

    Возвращает False для любой проблемы (файл не существует, пустой, не
    ZIP вовсе, повреждённая центральная директория, архив без единой
    записи) — никогда не бросает исключение наружу.
    """
    try:
        if not path.exists() or path.stat().st_size == 0:
            return False
        with zipfile.ZipFile(path, "r") as zf:
            names = zf.namelist()
        return len(names) > 0
    except (zipfile.BadZipFile, OSError, EOFError):
        return False



class MISAdapterError(RuntimeError):
    """Ошибка при работе с MIS."""


class MISAdapter:
    """
    Адаптер для взаимодействия с Michigan Imputation Server.

    Основные методы:
    - prepare_upload_files(): разбивка merged VCF на 23 файла (1-22 + X)
    - open_upload_page(): открытие страницы загрузки в браузере
    - download_results(): скачивание результатов по curl-ссылке
    - extract_all_results(): распаковка ВСЕХ ZIP-архивов результатов
    """

    MIS_URL = "https://imputationserver.sph.umich.edu"

    def __init__(
        self,
        upload_dir: Path,
        results_dir: Path,
        bcftools_path: Optional[str] = None,
        sevenzip_path: Optional[str] = None,
    ):
        """
        upload_dir: папка для подготовки файлов загрузки (upload/)
        results_dir: папка для сохранения результатов (rerun_results/)
        bcftools_path: явный путь к bcftools.exe (из --bin-dir приложения).
            Если не задан, ищем "bcftools" в системном PATH — этого может
            не хватить на Windows, если бинарник лежит в отдельной папке
            бандла, а не добавлен в PATH.
        sevenzip_path: явный путь к 7z.exe. Архивы MIS зашифрованы AES-256,
            который встроенный zipfile НЕ поддерживает (падает даже с
            верным паролем) — 7-Zip нужен как основной инструмент
            распаковки, а не опциональное ускорение.
        """
        self.upload_dir = Path(upload_dir)
        self.results_dir = Path(results_dir)
        self.upload_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.bcftools_path = bcftools_path or shutil.which("bcftools") or "bcftools"
        # Поиск 7z.exe теперь в одном месте — core/archive_utils.py
        # (Задача 6: раньше эта логика была продублирована и здесь, и в main.py).
        self.sevenzip_path = find_7z(sevenzip_path)

    def prepare_upload_files(
        self,
        merged_vcf: Path,
        chromosomes=None,
    ) -> list[Path]:
        """
        Разбивает merged VCF на отдельные файлы по хромосомам.

        merged_vcf: путь к объединённому VCF (batch_merged.vcf.gz)
        chromosomes: перечень хромосом (по умолчанию UPLOAD_CHROMS —
            1-22 + X). X добавлена вместе с поддержкой импутации
            X-хромосомы: Michigan Imputation Server сам делит присланный
            chrX.vcf.gz на PAR1/nonPAR/PAR2 и возвращает результат одним
            файлом, отдельной подготовки с нашей стороны не требуется.

        Возвращает список созданных файлов.
        """
        if chromosomes is None:
            chromosomes = UPLOAD_CHROMS
        merged_vcf = Path(merged_vcf)
        if not merged_vcf.exists():
            raise MISAdapterError(f"VCF файл не найден: {merged_vcf}")

        output_files: list[Path] = []

        for chrom in chromosomes:
            out_file = self.upload_dir / f"chr{chrom}.vcf.gz"

            cmd = [
                self.bcftools_path, "view",
                str(merged_vcf),
                # --targets (в отличие от --regions) не требует индекса
                # (.tbi/.csi) входного файла — merged_vcf у нас собирается
                # чисто-питоновским кодом и индекса не имеет.
                "--targets", str(chrom),
                "-Oz", "-o", str(out_file),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                output_files.append(out_file)
                logger.info("Создан файл: %s", out_file)
            except FileNotFoundError as e:
                raise MISAdapterError(
                    f"bcftools не найден ({self.bcftools_path!r}). "
                    f"Передайте bcftools_path в MISAdapter() или добавьте bcftools в PATH."
                ) from e
            except subprocess.CalledProcessError as e:
                stderr = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else e.stderr
                logger.error("Ошибка при создании chr%s: %s", chrom, stderr)
                raise MISAdapterError(f"Не удалось создать chr{chrom}: {stderr}") from e

        logger.info("Подготовлено %d файлов для загрузки на MIS", len(output_files))
        return output_files

    def open_upload_page(self) -> None:
        """Открывает страницу загрузки MIS в браузере."""
        logger.info("Открываю страницу Michigan Imputation Server...")
        webbrowser.open(self.MIS_URL)
        print("\n" + "="*70)
        print("ИНСТРУКЦИЯ ПО ЗАГРУЗКЕ:")
        print("="*70)
        print("1. Зарегистрируйтесь/войдите на сайте")
        print("2. Нажмите 'Genotype Imputation' → 'RUN'")
        print("3. В поле 'Name' введите любое название (например: genotek)")
        print("4. В 'Reference Panel' выберите: HRC r1.1 2016 (GRCh37/hg19)")
        print("5. Нажмите 'Select Files' и загрузите ВСЕ 23 файла (1-22 + X) из папки:")
        print(f"   {self.upload_dir.absolute()}")
        print("6. Поставьте галочки в необходимых параметрах")
        print("7. Нажмите 'Start Imputation'")
        print("8. Через 10-40 минут придёт письмо со ссылкой")
        print("="*70 + "\n")

    def _existing_file_verdict(self, dest: Path, file_url: str,
                               known: dict, foreign_job: bool = False) -> str:
        """
        Что делать с файлом, который уже лежит в папке результатов:
        "keep" — точно наш и целый, "adopt" — происхождение неизвестно, но
        размер сошёлся с сервером, иначе — короткая причина перекачать.

        Три уровня доверия, от дешёвого к дорогому:
          1. запись в манифесте с той же ссылкой и тем же размером —
             верим без единого сетевого запроса;
          2. манифеста нет (файлы от предыдущей версии программы) —
             спрашиваем у сервера длину и сверяем. Один HEAD-запрос
             дешевле, чем перекачивать гигабайты, и надёжнее, чем
             прежнее «валидный ZIP, значит наш»;
          3. всё остальное — перекачать.
        """
        if foreign_job:
            # Манифест прямо говорит, что папка осталась от ДРУГОГО задания
            # MIS. Спрашивать у сервера длину бессмысленно и опасно: размеры
            # архивов разных заданий вполне могут совпасть, и файл был бы
            # принят как «уже скачанный» — ровно тот баг, из-за которого
            # распаковка падала с «Wrong password» по всем 23 архивам.
            return "архив от другого задания MIS"

        try:
            size = dest.stat().st_size
        except OSError as e:
            return f"файл недоступен ({e})"
        if size == 0:
            return "нулевой размер"

        record = known.get(dest.name)
        if record:
            if record.get("url") != file_url:
                return "скачан по другой ссылке (другое задание MIS)"
            if record.get("size") != size:
                return "размер не совпал с записанным при скачивании"
            if not _is_valid_zip(dest):
                return "повреждён (не читается как ZIP)"
            return "keep"

        # Записи нет — происхождение файла неизвестно. Проверяем по длине
        # на сервере, а не по одному факту «это валидный ZIP»: архив от
        # ПРОШЛОГО задания MIS тоже валидный ZIP, но пароль к нему уже
        # другой, и распаковка упадёт по всем 23 архивам сразу.
        if not _is_valid_zip(dest):
            return "повреждён (не читается как ZIP)"
        expected = remote_size(file_url)
        if expected is None:
            return "не удалось сверить размер с сервером"
        if expected != size:
            return (f"размер не совпал с сервером "
                    f"({size} на диске против {expected})")
        return "adopt"

    def download_results(
        self,
        curl_command: str,
        on_file_error: Optional[Callable[[str, str], bool]] = None,
        cancel_check: Optional[Callable[[], bool]] = None,
    ) -> list[Path]:
        """
        Скачивает результаты по curl-команде из письма — БЕЗ bash.

        Старая версия выполняла `curl ... | bash` через shell=True. На
        Windows shell — это cmd.exe, где bash либо отсутствует, либо
        является bash из WSL, и тогда архивы скачивались внутрь файловой
        системы WSL, а не в results_dir. Поэтому ZIP-файлы «не находились».

        Теперь работаем напрямую:
        1. Извлекаем URL из команды вида `curl -sL <url> | bash`.
        2. Скачиваем bash-скрипт по этому URL как обычный текст.
        3. Вытаскиваем из скрипта ссылки на ZIP-архивы (по одному на
           хромосому плюс служебные, обычно ~22-27 файлов).
        4. Качаем каждый архив напрямую в results_dir через curl.

        curl_command: полная команда из письма, например:
            curl -sL https://imputationserver.sph.umich.edu/get/... | bash
            Также принимается просто голый URL.

        Промт "проверять уже скачанные файлы + предлагать повтор при
        ошибке + проверка целостности":
          - Перед скачиванием КАЖДОГО файла проверяется, нет ли его уже
            на диске (results_dir/<filename>). Если файл есть и
            структурно цел (валидный ZIP, см. _is_valid_zip()) — повторное
            скачивание пропускается. Если файл есть, но сломан
            (оборванная докачка с прошлого раза, битый архив) — он
            удаляется и качается заново, как будто его не было.
          - Сразу после скачивания нового файла он ТОЖЕ проверяется через
            _is_valid_zip() — недостаточно того, что curl вернул код 0 и
            файл непустой: content может оказаться HTML-страницей с
            ошибкой или оборванным потоком. Невалидный результат
            трактуется как ошибка скачивания этого файла (файл удаляется).
          - Ошибка скачивания ОДНОГО файла больше не прерывает всё
            скачивание немедленно — остальные файлы всё равно
            докачиваются. Если передан on_file_error(filename, error) ->
            bool, при неудаче конкретного файла он вызывается, и, если
            возвращает True, попытка для ЭТОГО ЖЕ файла повторяется
            (цикл длится, пока on_file_error не вернёт False или файл не
            скачается успешно и не пройдёт проверку целостности). Если
            on_file_error не передан (или вернул False) — файл
            добавляется в список неудавшихся, а обработка продолжается
            со следующего файла.
          - После обработки ВСЕХ файлов, если остались неудавшиеся,
            бросается MISAdapterError со списком всех проблемных файлов
            сразу (а не только первого, как раньше) — уже успешно
            скачанные и проверенные файлы при этом остаются на диске и
            будут пропущены при повторном вызове.

        Возвращает список путей ко ВСЕМ скачанным ZIP-файлам, реально
        присутствующим в results_dir на момент успешного завершения.
        """
        if not curl_command.strip():
            raise MISAdapterError("curl команда пуста")

        url_match = re.search(r"https?://[^\s|'\"]+", curl_command)
        if not url_match:
            raise MISAdapterError(
                f"Не удалось найти URL в команде: {curl_command!r}"
            )
        script_url = url_match.group(0)

        logger.info("Получаю скрипт скачивания: %s", script_url)
        try:
            res = subprocess.run(
                ["curl", "-sSfL", "--connect-timeout", "30", script_url],
                capture_output=True, text=True, errors="replace", check=True,
            )
        except FileNotFoundError as e:
            raise MISAdapterError(
                "curl не найден в системе. На Windows 10+ он встроен — "
                "проверьте PATH или обновите Windows."
            ) from e
        except subprocess.CalledProcessError as e:
            raise MISAdapterError(
                f"Не удалось получить скрипт скачивания: "
                f"{curl_error_text(e.returncode, e.stderr or '')}.\n"
                f"Чаще всего это истёкшая ссылка — результаты на MIS "
                f"хранятся около 3 дней после письма."
            ) from e

        script = res.stdout
        if not script or script.lstrip().lower().startswith("<"):
            raise MISAdapterError(
                "Вместо скрипта скачивания получена HTML-страница. "
                "Скорее всего срок действия ссылки истёк (результаты на MIS "
                "хранятся ~3 дня) — запустите задание на сервере заново."
            )

        # Формат-агностичный парсинг: скрипты разных серверов семейства
        # Cloudgene/eMIS оформлены по-разному.
        #   Michigan (imputationserver.sph.umich.edu): строки вида
        #     curl -sL https://.../chr_1.zip -o chr_1.zip
        #   BioDataCatalyst (imputation.biodatacatalyst.nhlbi.nih.gov):
        #     подтверждено реальным письмом — ссылки БЕЗ схемы вообще
        #     (imputation.biodatacatalyst.../share/results/...zip, а не
        #     https://imputation...), и без флага -O/-o. Раньше оба regex
        #     жёстко требовали https?:// в начале — из-за этого не
        #     срабатывал ни строгий формат (Style A), ни запасной
        #     (Style B), хотя нужные `curl ... -o chr_1.zip` строки в
        #     скрипте были. Схема теперь опциональна в обоих; "домен"
        #     обязателен (буквы/цифры/точки/дефисы + TLD из ≥2 букв),
        #     чтобы не подхватить случайные голые имена файлов вида
        #     "chr_1.zip" без пути.
        _DOMAIN = r'[a-zA-Z0-9](?:[a-zA-Z0-9.-]*[a-zA-Z0-9])?\.[a-zA-Z]{2,}'
        downloads = re.findall(
            r'(?:curl|wget)[^\n]*?'
            r'((?:https?://)?' + _DOMAIN + r'/[^\s"\']*\.zip)(?=[\s"\']|$)[^\n]*?'
            r'(?:-o|-O)\s*["\']?([^\s"\']+)',
            script,
        )
        if not downloads:
            # Style B (BioDataCatalyst и подобные): bare "wget <url>" или
            # голый URL (со схемой или без) без -o/-O, по одному на строку.
            zip_urls = list(dict.fromkeys(
                re.findall(
                    r'(?:https?://)?' + _DOMAIN + r'/[^\s"\']*\.zip(?=[\s"\']|$)',
                    script,
                )
            ))
            downloads = [
                (u, u.split("/")[-1].split("?")[0] or f"archive_{i}.zip")
                for i, u in enumerate(zip_urls)
            ]

        if not downloads:
            # Диагностика вместо «слепой» ошибки: показываем начало
            # реально полученного ответа — сразу видно, истекла ли
            # ссылка (HTML/логин-страница) или сервер использует ещё
            # какой-то третий формат скрипта.
            preview = script[:500].replace("\n", " | ")
            raise MISAdapterError(
                "В скрипте MIS не найдено ссылок на ZIP-архивы. Возможно, "
                "задание ещё не завершено, ссылка устарела, либо сервер "
                "вернул ответ неожиданного формата.\n"
                f"Начало полученного ответа: {preview}"
            )

        logger.info("В скрипте найдено %d архивов, начинаю скачивание...", len(downloads))

        # --- К какому ЗАДАНИЮ относятся файлы, уже лежащие в папке -------
        # Раньше файл пропускался как «уже скачанный» по одному признаку:
        # это валидный ZIP. Задание при этом не проверялось никак. Стоило
        # перезапустить задание на MIS (новый job id и, главное, НОВЫЙ
        # пароль), как архивы от прошлого задания честно проходили
        # проверку целостности, пропускались — и распаковка падала с
        # «Wrong password» сразу по всем 23 архивам, хотя пароль был верный
        # для нового задания, а файлы лежали от старого.
        job_id = job_id_from_urls([u for u, _ in downloads])
        manifest = read_manifest(self.results_dir)
        known: dict = manifest.get("files") or {}
        if manifest.get("job") and manifest["job"] != job_id:
            logger.warning(
                "⚠ В папке результатов лежат архивы от ДРУГОГО задания MIS "
                "(%s вместо %s) — они будут перекачаны. Пароль из письма "
                "подходит только к своему заданию, поэтому оставлять их "
                "нельзя.", manifest["job"], job_id,
            )
            foreign_job = True
            known = {}
            for stale in self.results_dir.glob("*.zip.part"):
                _safe_unlink(stale)
        else:
            foreign_job = False
        manifest = {"job": job_id, "files": known}

        already_present = 0
        redownloaded_broken = 0
        adopted = 0
        failed: list[tuple[str, str]] = []

        for file_url, filename in downloads:
            if cancel_check is not None and cancel_check():
                raise DownloadCancelled("Скачивание отменено пользователем")

            # Нормализация: если regex поймал ссылку без схемы (как у
            # BioDataCatalyst), curl без -L/схемы не поймёт, куда стучаться —
            # сайты MIS/eMIS-семейства HTTPS-only, поэтому дополняем сами.
            if not re.match(r'^https?://', file_url):
                file_url = "https://" + file_url
            dest = self.results_dir / filename

            if dest.exists():
                verdict = self._existing_file_verdict(
                    dest, file_url, known, foreign_job=foreign_job,
                )
                if verdict == "keep":
                    logger.info("✓ %s уже скачан и цел — пропускаю", filename)
                    already_present += 1
                    known[filename] = {"url": file_url,
                                       "size": dest.stat().st_size}
                    write_manifest(self.results_dir, manifest)
                    continue
                if verdict == "adopt":
                    logger.info(
                        "✓ %s уже лежит на диске, размер совпал с сервером — "
                        "принимаю как скачанный", filename,
                    )
                    adopted += 1
                    known[filename] = {"url": file_url,
                                       "size": dest.stat().st_size}
                    write_manifest(self.results_dir, manifest)
                    continue
                logger.warning(
                    "⚠ %s на диске не подходит (%s) — качаю заново",
                    filename, verdict,
                )
                _safe_unlink(dest)
                redownloaded_broken += 1

            attempt = 0
            while True:
                attempt += 1
                logger.info("Скачиваю: %s (попытка %d)", filename, attempt)
                try:
                    download_one(file_url, dest, cancel_check=cancel_check)
                    error_msg = None
                except (DownloadCancelled, FileLockedError):
                    # Отмена и занятый чужим процессом файл — не тот случай,
                    # где уместно предлагать «повторить этот файл»: повтор
                    # упрётся в то же самое. Пробрасываем наружу.
                    raise
                except MISAdapterError as e:
                    error_msg = str(e)

                if error_msg is None:
                    known[filename] = {"url": file_url,
                                       "size": dest.stat().st_size}
                    write_manifest(self.results_dir, manifest)
                    break

                logger.warning("⚠ %s", error_msg)
                if attempt >= MAX_FILE_ATTEMPTS:
                    error_msg += (
                        f"\n(исчерпаны {MAX_FILE_ATTEMPTS} попытки — "
                        f"уже скачанная часть сохранена, следующий запуск "
                        f"продолжит с этого места)"
                    )
                    failed.append((filename, error_msg))
                    break
                if on_file_error is not None and on_file_error(filename, error_msg):
                    continue
                failed.append((filename, error_msg))
                break

        write_manifest(self.results_dir, manifest)

        if already_present or adopted:
            logger.info(
                "✓ Пропущено как уже скачанное: %d (из них принято по "
                "совпадению размера с сервером: %d) из %d файлов",
                already_present + adopted, adopted, len(downloads),
            )
        if redownloaded_broken:
            logger.info(
                "ℹ Перекачано файлов, не прошедших проверку (битые или от "
                "другого задания): %d", redownloaded_broken,
            )

        zip_files = sorted(self.results_dir.glob("*.zip"))

        if failed:
            details = "\n".join(f"  - {name}: {err}" for name, err in failed)
            raise MISAdapterError(
                f"Не удалось скачать {len(failed)} из {len(downloads)} файлов:\n{details}\n\n"
                f"Уже успешно скачанные и проверенные файлы сохранены в "
                f"{self.results_dir} — повторный запуск пропустит их и "
                f"продолжит проблемные С МЕСТА ОБРЫВА, а не с нуля."
            )

        if not zip_files:
            raise MISAdapterError("ZIP файлы с результатами не найдены")

        logger.info("Скачано %d ZIP-архивов в: %s", len(zip_files), self.results_dir)
        return zip_files

    def extract_all_results(self, zip_paths: list[Path], password: str) -> None:
        """
        Распаковывает ВСЕ переданные ZIP-архивы результатов.

        Задача 6: логика поиска 7z.exe, санитайзинга пароля и самого вызова
        subprocess/zipfile теперь целиком живёт в core/archive_utils.py —
        здесь остаётся только адаптация MISAdapterError под интерфейс класса.
        """
        if not zip_paths:
            raise MISAdapterError("Список ZIP-архивов пуст")
        try:
            extract_all(
                [Path(p) for p in zip_paths],
                self.results_dir,
                password,
                sevenzip_path=self.sevenzip_path,
            )
        except ArchiveExtractionError as e:
            raise MISAdapterError(str(e)) from e
        logger.info("Все %d архивов результатов распакованы в: %s", len(zip_paths), self.results_dir)

    def extract_results(self, zip_path: Path, password: str) -> None:
        """Распаковывает ОДИН ZIP (оставлено для обратной совместимости).
        Для полного набора результатов MIS используйте extract_all_results()
        со списком из download_results() — иначе распакуется только одна
        хромосома из 23 (1-22 + X)."""
        self._extract_one(Path(zip_path), password)

    def _extract_one(self, zip_path: Path, password: str) -> None:
        try:
            extract_zip(Path(zip_path), self.results_dir, password, sevenzip_path=self.sevenzip_path)
        except ArchiveExtractionError as e:
            raise MISAdapterError(str(e)) from e

    def verify_results(self) -> dict[str, int]:
        """
        Проверяет наличие и целостность результатов.

        Возвращает словарь с количеством файлов:
        {
            "dose_vcf": 23,   # chr*.dose.vcf.gz (1-22 + X)
            "info": 23,        # chr*.info.gz
        }
        """
        dose_files = list(self.results_dir.glob("chr*.dose.vcf.gz"))
        info_files = list(self.results_dir.glob("chr*.info.gz"))

        result = {
            "dose_vcf": len(dose_files),
            "info": len(info_files),
        }

        expected = len(UPLOAD_CHROMS)
        if len(dose_files) != expected or len(info_files) != expected:
            # Не ошибка: X может отсутствовать (задание было отправлено
            # только с аутосомами, или панель/сервер не вернули X).
            logger.warning(
                "Ожидается %d файла каждого типа (1-22 + X), найдено: "
                "dose=%d, info=%d",
                expected, len(dose_files), len(info_files),
            )
        else:
            logger.info("Проверка результатов: OK (%d dose + %d info)",
                        expected, expected)

        return result


def main():
    """CLI для тестирования адаптера."""
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    parser = argparse.ArgumentParser(description="MIS Adapter")
    parser.add_argument(
        "--upload-dir",
        type=Path,
        default=Path("upload"),
        help="Папка для файлов загрузки",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("rerun_results"),
        help="Папка для результатов",
    )
    parser.add_argument(
        "--merged-vcf",
        type=Path,
        help="Путь к merged VCF для разбивки",
    )
    parser.add_argument(
        "--curl",
        type=str,
        help="curl команда для скачивания результатов",
    )
    parser.add_argument(
        "--password",
        type=str,
        help="Пароль для распаковки ZIP",
    )
    parser.add_argument(
        "--bin-dir",
        type=Path,
        default=None,
        help="Папка с bcftools.exe/7z.exe (бандл-бинарники приложения)",
    )

    args = parser.parse_args()

    bcftools_path = None
    sevenzip_path = None
    if args.bin_dir:
        cand = args.bin_dir / "bcftools.exe"
        if cand.is_file():
            bcftools_path = str(cand)
        cand7z = args.bin_dir / "7z.exe"
        if cand7z.is_file():
            sevenzip_path = str(cand7z)

    adapter = MISAdapter(
        args.upload_dir, args.results_dir,
        bcftools_path=bcftools_path, sevenzip_path=sevenzip_path,
    )

    if args.merged_vcf:
        print("=== Подготовка файлов для MIS ===")
        files = adapter.prepare_upload_files(args.merged_vcf)
        print(f"✓ Создано {len(files)} файлов в {args.upload_dir}")
        adapter.open_upload_page()

    if args.curl:
        print("=== Скачивание результатов ===")
        zip_paths = adapter.download_results(args.curl)
        print(f"✓ Скачано архивов: {len(zip_paths)}")

        if args.password:
            print("=== Распаковка ===")
            adapter.extract_all_results(zip_paths, args.password)
            print("✓ Распаковано")

    print("=== Проверка результатов ===")
    stats = adapter.verify_results()
    print(f"Dose VCF: {stats['dose_vcf']} файлов")
    print(f"Info: {stats['info']} файлов")


if __name__ == "__main__":
    main()
