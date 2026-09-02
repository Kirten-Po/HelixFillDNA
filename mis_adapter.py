"""
mis_adapter.py
Адаптер для работы с Michigan Imputation Server.
Подготовка 22 VCF-файлов для загрузки, открытие браузера,
скачивание и распаковка результатов.
"""
from __future__ import annotations
import logging
import re
import shutil
import subprocess
import webbrowser
import zipfile
from pathlib import Path
from typing import Callable, Optional

from core.pure_python_core import UPLOAD_CHROMS
from core.archive_utils import (
    extract_zip, extract_all, find_7z, ArchiveExtractionError,
)

logger = logging.getLogger(__name__)


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

    def download_results(
        self,
        curl_command: str,
        on_file_error: Optional[Callable[[str, str], bool]] = None,
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
                ["curl", "-sL", script_url],
                capture_output=True, text=True, check=True,
            )
        except FileNotFoundError as e:
            raise MISAdapterError(
                "curl не найден в системе. На Windows 10+ он встроен — "
                "проверьте PATH или обновите Windows."
            ) from e
        except subprocess.CalledProcessError as e:
            raise MISAdapterError(
                f"Не удалось получить скрипт скачивания (код {e.returncode}). "
                f"Возможно, срок действия ссылки истёк (результаты на MIS "
                f"хранятся ~3 дня после письма)."
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
        already_present = 0
        redownloaded_broken = 0
        failed: list[tuple[str, str]] = []

        for file_url, filename in downloads:
            # Нормализация: если regex поймал ссылку без схемы (как у
            # BioDataCatalyst), curl без -L/схемы не поймёт, куда стучаться —
            # сайты MIS/eMIS-семейства HTTPS-only, поэтому дополняем сами.
            if not re.match(r'^https?://', file_url):
                file_url = "https://" + file_url
            dest = self.results_dir / filename

            # Проверка "уже скачан и цел" — если файл уже лежит на диске
            # и проходит проверку целостности (валидный ZIP, не только
            # "непустой"), повторно его не качаем. Битый файл (оборванная
            # докачка/HTML-страница с ошибкой с прошлого раза) удаляется
            # и качается заново, как будто его не было.
            if dest.exists():
                if dest.stat().st_size > 0 and _is_valid_zip(dest):
                    logger.info("✓ %s уже скачан и цел (%d байт) — пропускаю", filename, dest.stat().st_size)
                    already_present += 1
                    continue
                logger.warning(
                    "⚠ %s уже есть на диске, но повреждён (не проходит проверку "
                    "целостности ZIP) — удаляю и качаю заново", filename,
                )
                dest.unlink(missing_ok=True)
                redownloaded_broken += 1

            while True:
                logger.info("Скачиваю: %s", filename)
                error_msg: Optional[str] = None
                try:
                    subprocess.run(
                        ["curl", "-sL", file_url, "-o", str(dest)],
                        cwd=str(self.results_dir), check=True,
                    )
                    if not dest.exists() or dest.stat().st_size == 0:
                        error_msg = f"Файл {filename} скачался пустым или не скачался вовсе."
                    elif not _is_valid_zip(dest):
                        error_msg = (
                            f"Файл {filename} скачался, но не проходит проверку "
                            f"целостности ZIP (похоже, скачался HTML-страницей с "
                            f"ошибкой или обрывом соединения)."
                        )
                except subprocess.CalledProcessError as e:
                    error_msg = f"Ошибка при скачивании {filename}: {e}"

                if error_msg is None:
                    break  # успех — переходим к следующему файлу

                dest.unlink(missing_ok=True)
                logger.warning("⚠ %s", error_msg)

                if on_file_error is not None and on_file_error(filename, error_msg):
                    # Пользователь (или вызывающий код) попросил повторить
                    # попытку именно для этого файла — качаем его ещё раз.
                    continue

                failed.append((filename, error_msg))
                break

        if already_present:
            logger.info(
                "✓ Уже было скачано ранее и прошло проверку целостности "
                "(пропущено): %d из %d файлов",
                already_present, len(downloads),
            )
        if redownloaded_broken:
            logger.info(
                "ℹ Обнаружено и перекачано повреждённых файлов с прошлого "
                "раза: %d", redownloaded_broken,
            )

        zip_files = sorted(self.results_dir.glob("*.zip"))

        if failed:
            details = "\n".join(f"  - {name}: {err}" for name, err in failed)
            raise MISAdapterError(
                f"Не удалось скачать {len(failed)} из {len(downloads)} файлов:\n{details}\n\n"
                f"Уже успешно скачанные и проверенные файлы сохранены в "
                f"{self.results_dir} — повторный запуск пропустит их и "
                f"попробует докачать только проблемные."
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
