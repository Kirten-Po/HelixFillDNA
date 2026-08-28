"""
gui/app.py
Графический интерфейс конвертера FTDNA/MyHeritage/VCF → 23andMe (Генотек).
Построен на CustomTkinter.

Изменения в этой версии (по промту "Комплексное улучшение GUI"):
 - Задача 1: горячие клавиши Ctrl+C/V/X/A больше не дублируют вставку —
   обработчики возвращают "break", гасящий встроенный биндинг виджета.
 - Задача 2: диагностика пароля и вся логика распаковки переиспользуют
   core/archive_utils.py вместо собственной копии.
 - Задача 6: прогресс-бар этапов 1-6 двигается плавно внутри каждого этапа
   через _set_subprogress(), а не только между этапами.
 - Задача 7: прогресс-бар этапа 7 показывает скачивание (0-50%) и сборку
   (50-100%); в stage_lbl — краткое сообщение, полные детали — в лог.
 - Задача A/B (устранение Invalid alleles + раздельные доноры по
   источникам): здесь БОЛЬШЕ НЕТ вызова pipeline._save_chip_signature()
   сразу после парсинга (Этап 1). Раньше сигнатура чипа перезаписывалась
   ДО проверки кэша доноров, из-за чего сравнение "cached == chip_signature"
   в _check_donors() было тривиально истинным всегда, и устаревший/чужой
   кэш доноров (donors/) принимался как валидный. Теперь _check_donors()
   — тонкая обёртка над pipeline.check_donor_cache(), которая только
   ЧИТАЕТ и СРАВНИВАЕТ существующую сигнатуру в donors/<source>/<panel>/,
   без права её перезаписывать; запись сигнатуры происходит только в
   download_donors.py после свежего скачивания.
 - Задача 1 ("Автопредложение скачать доноров через GUI"): _check_donors()
   на Этапе 3 в _run_stages_1_6() теперь обёрнут в _ensure_donors(), который
   при RuntimeError от check_donor_cache() (доноры отсутствуют/устарели)
   показывает диалог "Скачать автоматически?" и, при согласии, качает
   доноров прямо в этом же фоновом потоке через
   download_donors.download_donors_for_chip() (v11, без модульных глобалов,
   с инвалидацией устаревшего кэша при смене чипа). Диалог из фонового
   потока вызывается ТОЛЬКО через self.after(0, ...) + threading.Event —
   messagebox нельзя дёргать напрямую не из главного потока Tkinter.
   Добавлена кнопка "Отменить скачивание доноров" (self.cancel_donor_btn),
   активная только во время скачивания, привязанная к
   self._cancel_donor_download (threading.Event) — download_donors.py
   реально прерывает текущий subprocess (curl) через terminate()/kill(),
   а не ждёт окончания текущей хромосомы. CLI (main.py::main()) не
   затронут — RuntimeError из check_donor_cache() там по-прежнему
   завершает процесс как раньше.
 - Задача 2 ("Детекция несоответствия источника и файла"): в
   _run_stages_1_6() сразу после проверки, что csv_path существует, и ДО
   "[0/7] Проверка референсного генома" (ensure_reference_genome() может
   качать/проверять несколько ГБ, а источнику 'vcf' референс не нужен
   вовсе) вызывается pipeline.detect_source_from_file(csv_path). Если
   уверенность >= 0.8 и определённый источник отличается от выбранного в
   выпадающем списке — показывается thread-safe диалог с тремя вариантами
   (_prompt_source_mismatch(): "Сменить источник" / "Продолжить с
   выбранным" / "Отмена"), реализованный как CTkToplevel (messagebox
   поддерживает только 2 исхода). При "Сменить источник" обновляется и
   выпадающий список (визуально, через self.after), и локальная
   переменная `source`, которой дальше пользуется весь остальной код
   _run_stages_1_6() — повторный парсинг файла не требуется, потому что
   детекция вызывается строго до Этапа 1. Попутно _check_donors() получил
   необязательный параметр `source`: раньше он всегда читался из
   self.source_dd через self._get_source_key(), а после появления
   переключения источника в середине _run_stages_1_6() это создавало
   гонку (обновление виджета уходит через self.after(0, ...) асинхронно
   в главный поток) — теперь _ensure_donors() передаёт source явно.
 - Шаг 1 промта "Поддержка выбора референсной панели HRC / TopMed r3":
   добавлен выпадающий список self.panel_dd ("Референсная панель:") —
   HRC r1.1 2016 (GRCh37/hg19) / TOPMed r3 (GRCh38/hg38), по умолчанию
   HRC (поведение без выбора не меняется). Выбранная панель прокидывается
   явным параметром panel через pipeline.ensure_reference_genome(),
   pipeline._donor_source_dir() (через save_pos_fn/positions_cache_dir),
   _check_donors()/_ensure_donors() — по аналогии с тем, как уже
   прокидывается source. Текст инструкции на вкладке "Запуск" (Reference
   Panel: ...) обновляется динамически при смене выбора в self.panel_dd.
   ⚠ Реальный лифтовер координат GRCh37<->GRCh38 для TopMed на ЭТОМ ШАГЕ
   (Шаг 1) ещё не был реализован — выбор TopMed переключал пути на диске
   и скачивание референса, но практическое использование до завершения
   следующих шагов плана не рекомендовалось.
   Промт "встроить лифтовер HRC/TopMed в gui/app.py": это ограничение
   снято — лифтовер (core/liftover.py::ChainLiftover) реализован и
   подключён и в _run_stages_1_6() (форвард), и в _run_stage_7()
   (обратно), см. соответствующие комментарии там и в _on_panel_changed().
   Оставшееся реальное ограничение — source='vcf' лифтовер не
   поддерживает (см. pipeline._supports_liftover()).
"""
from __future__ import annotations

import contextlib
import io
import json
import logging
import os
import pickle
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
import tkinter as tk
import traceback
import webbrowser
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox
from urllib.parse import quote

import customtkinter as ctk

# === Корень проекта ===
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

# Версия приложения — из version.py в корне проекта (см. его докстринг:
# та же строка сверяется с git-тегом при сборке релиза). Импорт возможен
# только ПОСЛЕ sys.path.insert() выше.
from version import __version__

# === Фикс "--- Logging error ---" / AttributeError: 'NoneType' object has
# no attribute 'write' в windowed-сборке PyInstaller (console=False) ===
# main.py при импорте (см. "Импорты пайплайна" ниже) вызывает
# logging.basicConfig(...) без явного stream= — это создаёт
# logging.StreamHandler(), который привязывается к sys.stderr МОМЕНТАЛЬНО,
# в момент вызова, и хранит эту привязку (self.stream) до конца работы
# процесса. В windowed-сборке (console=False в HelixFillDNA.spec) у
# процесса нет консоли, поэтому на старте sys.stdout/sys.stderr — None, а
# не поток. Если не подменить их ДО импорта main.py, handler привяжется к
# None и будет падать на каждый logger.info()/warning() с
# "AttributeError: 'NoneType' object has no attribute 'write'" — сам
# logging перехватывает эту ошибку внутри Handler.handleError() (поэтому
# приложение не падает), но выводит трейсбек в лог при каждом вызове.
# Подмена sys.stdout/sys.stderr на LogRedirector ниже (в _run_stages_1_6
# и т.п.) не помогает — она подменяет sys.stderr ПОСЛЕ того, как
# basicConfig() уже создал handler и запомнил старое значение (None).
# Поэтому подменяем None на no-op поток здесь, ДО import main as pipeline.
if sys.stdout is None:
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w", encoding="utf-8")


# ---------------------------------------------------------------------------
# Промт "скачивание должно быть видно в приложении, а не в отдельных окнах
# cmd".
#
# В оконной сборке PyInstaller (console=False) у GUI-процесса нет своей
# консоли. Когда такой процесс запускает КОНСОЛЬНОЕ приложение
# (aria2c.exe, curl.exe, bcftools.exe, tabix.exe, bgzip.exe, 7z.exe),
# Windows создаёт для него НОВОЕ консольное окно: у долгих скачиваний оно
# висит на экране отдельным чёрным окном, у коротких вызовов bcftools —
# мельтешит вспышками на каждый вызов. Прогресс при этом идёт мимо лога
# приложения.
#
# В проекте таких вызовов больше тридцати (subprocess.run/Popen в
# main.py, download_donors.py, core/*, template/*) — проставлять флаг в
# каждом месте пришлось бы руками и легко забыть при следующей правке.
# Поэтому флаг ставится ОДИН раз здесь, в точке входа, до импорта любого
# кода проекта: subprocess.run()/check_output() внутри создают тот же
# Popen, так что патч конструктора покрывает все вызовы разом.
#
# Явно переданный вызывающим кодом creationflags сохраняется (флаг только
# добавляется битовой маской), а на не-Windows этот блок не делает ничего.
# ---------------------------------------------------------------------------
if os.name == "nt":
    import subprocess as _subprocess

    _CREATE_NO_WINDOW = 0x08000000
    _orig_popen_init = _subprocess.Popen.__init__

    def _popen_init_no_window(self, *args, **kwargs):
        kwargs["creationflags"] = kwargs.get("creationflags", 0) | _CREATE_NO_WINDOW
        return _orig_popen_init(self, *args, **kwargs)

    _subprocess.Popen.__init__ = _popen_init_no_window


def _detect_bin_dir() -> Path:
    """
    Возвращает папку с htslib-бинарниками (bcftools/tabix/bgzip).

    Раньше поле "Папка с бинарниками htslib" в GUI жёстко предзаполнялось
    как PROJECT_ROOT / "bin", в расчёте на плоский макет (bin/ рядом с
    exe). В собранном PyInstaller'ом .exe это предположение не всегда
    верно: начиная с PyInstaller 6.0 бинарники по умолчанию попадают в
    PROJECT_ROOT / "_internal" / "bin" (см. HelixFillDNA.spec —
    contents_directory должен возвращать плоский макет, но на практике
    поведение зависит от версии PyInstaller и надёжно проверять его на
    каждой машине разработчика непрактично).

    Вместо того чтобы полагаться на конкретный макет сборки, проверяем
    оба возможных расположения и выбираем то, где реально лежит
    bcftools(.exe) — так поле в GUI всегда предзаполняется верно,
    независимо от того, как именно PyInstaller разложил файлы в этой
    конкретной сборке. Если бинарник не найден ни там, ни там (например,
    запуск из исходников до первой сборки) — возвращаем прежний плоский
    путь как запасной вариант, поведение не меняется.
    """
    exe_name = "bcftools.exe" if os.name == "nt" else "bcftools"
    flat_bin = PROJECT_ROOT / "bin"
    if (flat_bin / exe_name).is_file():
        return flat_bin
    internal_bin = PROJECT_ROOT / "_internal" / "bin"
    if (internal_bin / exe_name).is_file():
        return internal_bin
    return flat_bin


def _find_app_icon() -> Path | None:
    """
    Путь к app_icon.ico для иконки окна (та, что слева от заголовка и на
    панели задач). Ищем в тех же двух местах, что и bin/ и samples/ —
    PyInstaller кладёт файлы либо плоско рядом с exe, либо в _internal/,
    а установщик Inno Setup кладёт иконку в {app} (см. [Files] в
    HelixFillDNA.iss). None, если файла нет — окно просто останется с
    иконкой Tk по умолчанию, падать из-за этого приложение не должно.
    """
    for candidate in (
        PROJECT_ROOT / "app_icon.ico",
        PROJECT_ROOT / "_internal" / "app_icon.ico",
    ):
        if candidate.is_file():
            return candidate
    return None


# ---------------------------------------------------------------------------
# Промт "обычные / продвинутые настройки".
#
# Вкладка "Подготовка" исторически показывала СРАЗУ все параметры пайплайна
# (формат вывода, порог Rsq, нормализация multiallelic, число EUR-доноров,
# кэш сырых хромосом, путь к трафарету и к бинарникам htslib). Для человека,
# который просто хочет "дополнить свой FTDNA-файл", это стена настроек, где
# каждая ошибка стоит часов перекачки доноров.
#
# Поэтому у вкладки теперь два режима:
#   * MODE_SIMPLE   — видны только источник данных, референсная панель и файл
#                     с данными. Всё остальное выставляется автоматически
#                     (см. _apply_simple_presets()), включая трафарет, который
#                     берётся из папки samples/ рядом с программой.
#   * MODE_ADVANCED — прежнее поведение, все настройки видны и правятся руками.
#
# Обычный режим НЕ подменяет значения "на лету" в момент запуска: он
# физически проставляет их в те же самые виджеты, что и раньше. Поэтому весь
# код запуска (_get_format_key(), _get_rsq_threshold(), _get_eur_sample_count()
# и т.д.) читает настройки ровно как прежде и ничего не знает о режимах, а
# пользователь, переключившись в продвинутый режим, видит именно то, что
# будет использовано.
# ---------------------------------------------------------------------------
# Подписи подвкладок на вкладке "Запуск". Нумерация — не этапы пайплайна
# (их 7), а то, что видит пользователь: подготовка файлов, ручная работа на
# сайте MIS, сборка итогового файла. Базовые имена без галочки; актуальные
# (с "✓" у пройденных) живут в App._run_tab_names, см. _set_wizard_step().
RUN_TAB_BASE_NAMES = (
    "1 · Подготовка файлов",
    "2 · Импутация на MIS",
    "3 · Сборка файла",
)

# Сколько ячеек хромосом в ряду карты этапа 3.
DONOR_GRID_COLUMNS = 4

# Идентификатор задания в письме MIS (job-20260828-123456 / -1). Нужен
# только для дружелюбного "✓ Распознано задание ..." — на работу
# скачивания не влияет, поэтому несовпадение не считается ошибкой.
_MIS_JOB_RE = re.compile(r"job-\d{8}-\d{6}(?:-\d+)?", re.IGNORECASE)


def _fmt_duration(seconds: float) -> str:
    """Человеческая длительность: «меньше минуты», «12 мин», «1 ч 20 мин»."""
    seconds = max(0.0, seconds)
    minutes = int(seconds // 60)
    if minutes < 1:
        return "меньше минуты"
    if minutes < 60:
        return f"{minutes} мин"
    hours, rest = divmod(minutes, 60)
    return f"{hours} ч {rest:02d} мин" if rest else f"{hours} ч"



# ---------------------------------------------------------------------------
# Промт "живая карта хромосом на этапе 3".
#
# Скачивание доноров идёт в несколько потоков и часами. Общий счётчик
# «Обработано хромосом: 7/22» не отвечает на главный вопрос пользователя —
# «оно вообще шевелится или зависло?». download_donors.py уже печатает
# подробности по каждой хромосоме (скачивание с процентами, фильтрация,
# индексация, ошибки), просто они тонули в потоке лога. Ниже — разбор этих
# строк в состояние конкретной хромосомы для карты 22 ячеек на вкладке
# "Запуск". Менять сам download_donors.py при этом не потребовалось.
#
# Ранг нужен из-за параллельности: строка «скачивание 41%» может прийти в
# очередь ПОЗЖЕ, чем «готово» (её напечатал другой поток чуть раньше), и без
# ранга ячейка откатывалась бы назад. Состояние может только расти.
# ---------------------------------------------------------------------------
_DONOR_RANK_START = 1
_DONOR_RANK_DOWNLOAD = 2
_DONOR_RANK_FILTER = 3
_DONOR_RANK_DONE = 4
_DONOR_RANK_FAILED = 5  # выше "готово": ошибку не должно перекрывать ничем

_DONOR_CHR_RE = re.compile(r"chr(\d{1,2})\b", re.IGNORECASE)
# Процент может быть дробным (curl --progress-bar печатает "30.7%"), и без
# необязательной дробной части регэксп цеплялся бы за "7%" вместо "30%".
_DONOR_PCT_RE = re.compile(r"(\d{1,3})(?:[.,]\d+)?\s*%")

# Как часто опрашивать размеры файлов доноров на диске, мс.
_DONOR_WATCH_INTERVAL_MS = 1500


def _fmt_size(num_bytes: float) -> str:
    """Размер/скорость по-человечески: «812 МБ», «1.2 ГБ», «340 КБ»."""
    num_bytes = max(0.0, float(num_bytes))
    for unit, limit in (("ГБ", 1024 ** 3), ("МБ", 1024 ** 2), ("КБ", 1024)):
        if num_bytes >= limit:
            value = num_bytes / limit
            return f"{value:.1f} {unit}" if value < 10 else f"{value:.0f} {unit}"
    return f"{num_bytes:.0f} Б"


def _parse_donor_state(msg: str):
    """
    (номер_хромосомы, текст, цвет, ранг) для строки лога скачивания
    доноров, либо None — если строка не про конкретную хромосому.
    """
    match = _DONOR_CHR_RE.search(msg)
    if not match:
        return None
    try:
        chrom = int(match.group(1))
    except ValueError:
        return None
    if not 1 <= chrom <= 22:
        return None

    low = msg.lower()

    if msg.lstrip().startswith(("✗", "❌")) or "не удалось" in low or "ошибка" in low:
        return chrom, "✗ ошибка", "#F44336", _DONOR_RANK_FAILED

    if "уже готов" in low or "отфильтров" in low or "проверено" in low:
        return chrom, "✓ готово", "#4CAF50", _DONOR_RANK_DONE

    if "ALL.chr" in msg:
        pct = _DONOR_PCT_RE.search(msg)
        if pct:
            return (chrom, f"⬇ скачивание {pct.group(1)}%", "#42A5F5",
                    _DONOR_RANK_DOWNLOAD)
        if "индекс" in low:
            return chrom, "⬇ индекс", "#42A5F5", _DONOR_RANK_DOWNLOAD

    if "подвыборк" in low or "chrom донора" in low or "фильтр" in low:
        return chrom, "⚙ фильтрация", "#F9A825", _DONOR_RANK_FILTER

    if "---" in msg:
        kind = "удалённо" if "удалённо" in low else "начата"
        return chrom, f"… {kind}", "gray70", _DONOR_RANK_START

    return None


# ---------------------------------------------------------------------------
# Промт "итоговый файл в отдельной папке".
#
# Раньше собранный файл ложился прямо в output/runs/<запуск>/, вперемешку с
# десятками промежуточных VCF, логов и служебных json — найти его там было
# отдельной задачей. Теперь итоговые файлы ВСЕХ запусков складываются в одну
# папку results/ рядом с программой, а имя файла начинается с названия
# запуска, чтобы они не путались между собой.
# ---------------------------------------------------------------------------
RESULTS_DIR_NAME = "results"


# ---------------------------------------------------------------------------
# Промт "обратная связь автору".
#
# Адрес получателя. ЕДИНСТВЕННОЕ место, где он задан, — меняется здесь.
#
# Программа не отправляет письмо сама, а открывает почтовый клиент
# пользователя заготовленным mailto:. Это сознательный отказ от варианта
# "приложение шлёт письмо через SMTP": логин и пароль от ящика пришлось бы
# положить внутрь распространяемого exe, откуда их извлекает кто угодно за
# пять минут — и рассылает с этого ящика спам, пока его не заблокируют.
# Вариант со сторонним сервисом форм (Formspree и т.п.) такой проблемы не
# создаёт, но требует аккаунта сервиса и делает отправку невидимой для
# пользователя. Через mailto пользователь видит письмо целиком и жмёт
# "Отправить" сам.
# ---------------------------------------------------------------------------
FEEDBACK_EMAIL = "kirten-tempest2026@outlook.com"
FEEDBACK_SUBJECT_PREFIX = "HelixFillDNA"
FEEDBACK_KINDS = ("Ошибка", "Предложение")
FEEDBACK_LOG_LINES = 40

# Windows передаёт mailto: через командную строку, где длина ограничена, и
# слишком длинное письмо может обрезаться на полуслове. Считать надо длину
# ГОТОВОЙ ссылки, а не текста: при percent-кодировании одна кириллическая
# буква превращается в девять символов (%D0%9E).
#
# Порог подобран по реальным замерам, а не «на глаз» — первая версия с
# лимитом 1900 срабатывала ВСЕГДА, даже на пустом письме, и пользователь
# получал «текст не поместился», ничего не написав:
#   пустое описание + техданные          -> ссылка ~2300
#   описание на 500 знаков + техданные   -> ссылка ~5000
#   то же плюс 40 строк лога             -> ссылка ~12900
# 7000 пропускает обычное письмо с содержательным описанием и отправляет
# в файл только то, что действительно огромно, — прежде всего лог.
FEEDBACK_MAILTO_URL_LIMIT = 7000


def _results_dir() -> Path:
    """Папка с итоговыми файлами; создаётся при первом обращении."""
    target = PROJECT_ROOT / RESULTS_DIR_NAME
    try:
        target.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    return target


def _unique_result_path(path: Path) -> Path:
    """
    Путь, который точно не затрёт уже существующий файл: к имени
    добавляется _2, _3 и т.д. Повторная сборка того же запуска (например с
    другим порогом Rsq) не должна молча уничтожать предыдущий результат.
    """
    if not path.exists():
        return path
    stem, suffix = path.stem, path.suffix
    n = 2
    while True:
        candidate = path.with_name(f"{stem}_{n}{suffix}")
        if not candidate.exists():
            return candidate
        n += 1


# Число этапов, которые выполняет кнопка "Запустить этапы 1-6 (до MIS)" —
# делитель шкалы прогресса на вкладке "Запуск".
STAGES_TOTAL = 6

MODE_SIMPLE = "Обычные настройки"
MODE_ADVANCED = "Продвинутые настройки"

# Источник данных -> формат вывода в обычном режиме. FTDNA Family Finder
# исторически соответствует трафарету v3 (LF), MyHeritage — v5 (CRLF).
# Готовый VCF не привязан ни к одному из экспортов, поэтому для него берём
# v3 как более полный по call rate.
#
# AncestryDNA -> v3 по той же причине, что и VCF, но с конкретным
# основанием: чип Ancestry V2.0 пересекается с трафаретом v3 по 473 136
# позициям (49,3 % трафарета) и лишь по 167 684 с v5 (26,1 %) — то есть
# на v3 у него почти втрое больше собственных измерений, которые не
# придётся импутировать.
SIMPLE_FORMAT_BY_SOURCE = {
    "ftdna": "v3", "myheritage": "v5", "ancestry": "v3", "vcf": "v3",
}
SIMPLE_RSQ = "0.30"          # стандартный порог MIS
SIMPLE_EUR_COUNT = 20        # компромисс трафик/качество для обычного режима
SIMPLE_NORMALIZE = True      # нормализовать multiallelic-сайты перед split
SIMPLE_RAW_CACHE = True      # хранить сырые хромосомы 1000 Genomes

# Имена трафаретов в папке samples/ — по одному на формат вывода.
SAMPLE_TEMPLATE_NAMES = {"v3": "template_v3.txt", "v5": "template_v5.txt"}

# Выбранный режим переживает перезапуск программы: маленький JSON рядом с
# exe, а не реестр/AppData — программа и так держит свои данные (donors/,
# output/, reference/) в своей папке.
UI_STATE_FILE = PROJECT_ROOT / "ui_state.json"

SAMPLES_README = """Папка samples — трафареты для сборки итогового файла.

Положите сюда файлы:
    template_v3.txt  — трафарет для источника FTDNA Family Finder (формат v3, LF)
    template_v5.txt  — трафарет для источника MyHeritage (формат v5, CRLF)

Трафарет — это реальный экспорт 23andMe соответствующей версии: программа
берёт из него порядок и набор rsid/позиций, а генотипы подставляет ваши.

В режиме "Обычные настройки" программа сама подставляет нужный трафарет из
этой папки по выбранному источнику данных, поэтому выбирать файл вручную не
нужно. В режиме "Продвинутые настройки" путь к трафарету по-прежнему можно
указать вручную кнопкой "Обзор".
"""


def _samples_candidates() -> list[Path]:
    """
    Возможные расположения папки samples/ — тот же приём, что и в
    _detect_bin_dir(): PyInstaller в зависимости от версии кладёт datas либо
    плоско рядом с exe, либо в _internal/.
    """
    return [PROJECT_ROOT / "samples", PROJECT_ROOT / "_internal" / "samples"]


def _samples_dir() -> Path:
    """
    Папка с трафаретами. Если её нет ни в одном из ожидаемых мест —
    создаёт PROJECT_ROOT/samples и кладёт туда README с инструкцией, куда
    какой трафарет положить (в установщик трафареты входят, но при запуске
    из исходников или после ручного удаления папки её нужно воссоздать).
    Ошибки создания (запуск из read-only каталога) намеренно проглатываются:
    отсутствие папки не должно мешать продвинутому режиму, где путь к
    трафарету указывается вручную.
    """
    for candidate in _samples_candidates():
        if candidate.is_dir():
            return candidate
    target = PROJECT_ROOT / "samples"
    try:
        target.mkdir(parents=True, exist_ok=True)
        readme = target / "README.txt"
        if not readme.exists():
            readme.write_text(SAMPLES_README, encoding="utf-8")
    except OSError:
        pass
    return target


def _find_sample_template(fmt: str) -> Path | None:
    """Путь к трафарету samples/template_<fmt>.txt или None, если его нет."""
    name = SAMPLE_TEMPLATE_NAMES.get(fmt)
    if not name:
        return None
    for base in _samples_candidates():
        candidate = base / name
        if candidate.is_file():
            return candidate
    return None


def _load_ui_mode() -> str:
    """Режим вкладки "Подготовка" из прошлого запуска (по умолчанию — обычный)."""
    try:
        data = json.loads(UI_STATE_FILE.read_text(encoding="utf-8"))
        mode = data.get("settings_mode")
    except (OSError, ValueError):
        return MODE_SIMPLE
    return mode if mode in (MODE_SIMPLE, MODE_ADVANCED) else MODE_SIMPLE


def _save_ui_mode(mode: str) -> None:
    """Сохраняет выбранный режим. Не критично: ошибки записи игнорируются."""
    try:
        data = {}
        if UI_STATE_FILE.is_file():
            with contextlib.suppress(ValueError):
                data = json.loads(UI_STATE_FILE.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            data = {}
        data["settings_mode"] = mode
        UI_STATE_FILE.write_text(
            json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
        )
    except OSError:
        pass


# === Импорты пайплайна ===
import main as pipeline
import download_donors
from adapters.ftdna_v3 import ReferenceGenome
from core import archive_utils
from core import network_utils


# ---------------------------------------------------------------------------
# Промт "простой лог скачивания — одна обновляющаяся строка на файл":
# _pump_progress() в download_donors.py прореживает вывод aria2c/curl по
# времени (не чаще раза в 2с НА ОДИН файл), но при нескольких хромосомах,
# качающихся параллельно (DEFAULT_PARALLEL_CHROMOSOMES), плюс отдельные
# .tbi-индексы — в лог всё равно сыпались десятки перемежающихся строк
# вида "aria2c ALL.chr14...: [...41%...]". Ниже — распознавание таких
# строк по регэкспу и обновление ОДНОЙ и той же строки в логе (через
# именованные метки Tkinter Text, а не insert("end", ...)) вместо
# добавления новой на каждое обновление процента. Ключ — номер хромосомы
# + тип файла (сам VCF донора или его .tbi-индекс), поэтому при
# параллельном скачивании нескольких хромосом каждая держит свою строку.
_CHR_PROGRESS_RE = re.compile(r"ALL\.chr(\d+)\.")
_PERCENT_RE = re.compile(r"(\d{1,3})%")


def _parse_progress_line(msg: str) -> tuple[str, str] | None:
    """
    Возвращает (ключ, текст_для_показа) для строки прогресса скачивания
    конкретной хромосомы (aria2c/curl), либо None, если msg — обычное
    сообщение (этап, ошибка, готово и т.п.), которое нужно просто
    добавить в лог как есть.
    """
    if "ALL.chr" not in msg:
        return None
    chr_match = _CHR_PROGRESS_RE.search(msg)
    pct_match = _PERCENT_RE.search(msg)
    if not chr_match or not pct_match:
        return None
    chrom = chr_match.group(1)
    is_index = ".tbi" in msg
    key = f"chr{chrom}_idx" if is_index else f"chr{chrom}"
    kind = "индекс" if is_index else "донор"
    text = f"⬇ chr{chrom} ({kind}): {pct_match.group(1)}%"
    return key, text


# ---------------------------------------------------------------------------
# Перехват print() → очередь → GUI-лог
# ---------------------------------------------------------------------------
class UserCancelledRun(RuntimeError):
    """
    Промт "Доноры для VCF-источника: понятная отмена + общий кэш сырых
    хромосом", Шаг 3: маркерный класс исключений для ШТАТНЫХ, ожидаемых
    прерываний запуска по решению самого пользователя (отказ скачивать
    доноров, отмена при несоответствии источника и т.п.) — в отличие от
    НЕОЖИДАННЫХ ошибок (сеть, bcftools, повреждённый файл), для которых
    полный traceback реально помогает диагностике.

    _run_stages_1_6() проверяет тип исключения и для UserCancelledRun
    печатает только короткое дружелюбное сообщение без сырого Python
    traceback — раньше даже полностью штатный отказ ("Нет" в диалоге
    скачивания доноров) выглядел в логе как настоящий программный сбой.
    """


class LogRedirector(io.TextIOBase):
    """
    Перенаправляет stdout/stderr в очередь для безопасного вывода в GUI.

    log_file_path (промт "Именованные папки запуска", опционально):
    если задан, каждое сообщение дополнительно дописывается в этот файл
    (обычно <run_dir>/run.log) — logger.*-сообщения так не дублируются
    (см. pipeline.attach_run_log_handler(), у него отдельный
    FileHandler на root logger), а вот все print()-сообщения (и свои, и
    из pipeline/download_donors и т.д., пока действует
    contextlib.redirect_stdout/stderr) — дублируются именно здесь.
    """
    def __init__(self, q: queue.Queue, log_file_path: Path | None = None):
        self.q = q
        self._log_file = None
        if log_file_path is not None:
            try:
                self._log_file = open(log_file_path, "a", encoding="utf-8")
            except OSError:
                self._log_file = None

    def write(self, s):
        if s and s.strip():
            clean = s.replace("\r", "").rstrip("\n")
            if clean:
                self.q.put(clean)
                if self._log_file is not None:
                    try:
                        self._log_file.write(clean + "\n")
                        self._log_file.flush()
                    except OSError:
                        pass
        return len(s or "")

    def flush(self):
        pass

    def close(self):
        if self._log_file is not None:
            try:
                self._log_file.close()
            except OSError:
                pass
            self._log_file = None


# ---------------------------------------------------------------------------
# Горячие клавиши и контекстное меню
# ---------------------------------------------------------------------------
def _is_textbox(widget) -> bool:
    return isinstance(widget, ctk.CTkTextbox) or isinstance(widget, tk.Text)


def _get_clipboard(widget):
    try:
        return widget.clipboard_get()
    except (tk.TclError, AttributeError):
        try:
            return widget.winfo_toplevel().clipboard_get()
        except tk.TclError:
            return None


def _do_copy(widget):
    if _is_textbox(widget):
        try:
            text = widget.get("sel.first", "sel.last")
            widget.clipboard_clear()
            widget.clipboard_append(text)
        except Exception:
            # Нет выделения (TclError) или CTk-обёртка не поддержала
            # какой-то из вызовов (AttributeError и т.п.) — в любом
            # случае просто нечего копировать, тихо выходим.
            pass
    else:
        try:
            text = widget.selection_get()
            if text:
                widget.clipboard_clear()
                widget.clipboard_append(text)
        except Exception:
            pass


def _do_cut(widget):
    if _is_textbox(widget):
        try:
            text = widget.get("sel.first", "sel.last")
            widget.delete("sel.first", "sel.last")
            widget.clipboard_clear()
            widget.clipboard_append(text)
        except Exception:
            pass
    else:
        try:
            text = widget.selection_get()
            if text:
                start = widget.index("sel.first")
                end = widget.index("sel.last")
                widget.delete(start, end)
                widget.clipboard_clear()
                widget.clipboard_append(text)
        except Exception:
            # Раньше здесь ловился только tk.TclError — если CTkEntry не
            # поддерживает .index()/.selection_get() напрямую (AttributeError),
            # исключение улетало необработанным внутрь Tkinter-обработчика
            # события и вырезание молча ничего не делало, без видимой ошибки.
            pass


def _do_paste(widget):
    """
    Задача "Ctrl+C/Ctrl+V для пароля и ссылки не работали": раньше
    финальный widget.insert("insert", text) для НЕ-textbox виджетов
    (CTkEntry — как раз поля пароля и curl-команды) не был обёрнут в
    try/except вовсе, а блок try/except перед ним ловил только
    tk.TclError, не AttributeError. Если CTkEntry.index()/.insert() не
    отработали ожидаемым образом (у CustomTkinter это composite-виджет
    поверх настоящего tk.Entry, и не все стандартные методы Entry
    гарантированно доступны/ведут себя как в голом tkinter), исключение
    улетало необработанным внутрь Tkinter-обработчика события — Ctrl+V
    просто ничего не делал, без единого сообщения об ошибке.

    Теперь: любая проблема на пути "вставить по выделению/курсору"
    перехватывается широко (Exception, а не только TclError), и есть
    надёжный запасной вариант — заменить содержимое поля целиком. Для
    полей пароля/ссылки это и есть основной сценарий использования
    (вставить скопированное целиком), так что запасной вариант не
    ухудшает UX, а гарантирует, что вставка вообще срабатывает.
    """
    text = _get_clipboard(widget)
    if not text:
        return

    if _is_textbox(widget):
        try:
            widget.delete("sel.first", "sel.last")
        except Exception:
            pass
        try:
            widget.insert("insert", text)
        except Exception:
            # Запасной вариант — заменить всё содержимое textbox целиком.
            try:
                widget.delete("1.0", "end")
                widget.insert("1.0", text)
            except Exception:
                pass
    else:
        try:
            start = widget.index("sel.first")
            end = widget.index("sel.last")
            widget.delete(start, end)
            widget.insert("insert", text)
        except Exception:
            # Нет выделения, либо CTkEntry не поддержал .index()/.insert()
            # с этими аргументами — заменяем содержимое поля целиком, это
            # и есть типичный сценарий вставки пароля/ссылки.
            try:
                widget.delete(0, "end")
                widget.insert(0, text)
            except Exception:
                pass


def _do_select_all(widget):
    if _is_textbox(widget):
        try:
            widget.tag_add("sel", "1.0", "end")
            widget.mark_set("insert", "end")
        except Exception:
            pass
    else:
        try:
            widget.select_range(0, "end")
            widget.icursor("end")
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Промт "горячие клавиши работают не до конца".
#
# Раньше Ctrl+C/V/X/A привязывались как <Control-c>, <Control-v> и т.д. На
# Windows Tk подставляет в event.keysym символ ТЕКУЩЕЙ раскладки: при
# русской раскладке нажатие физической клавиши C даёт keysym "es" (или
# Cyrillic_es), а не "c" — событие <Control-c> просто не срабатывает, и
# копирование/вставка «не работали», причём молча и через раз (в
# английской раскладке всё было в порядке, поэтому баг выглядел
# плавающим). Вариант с перечислением кириллических keysym'ов пришлось бы
# расширять под каждую новую раскладку пользователя.
#
# Решение — одна привязка <Control-KeyPress> и разбор по event.keycode,
# то есть по ФИЗИЧЕСКОЙ клавише (на Windows это VK-код: C=67, V=86, X=88,
# A=65), не зависящей от раскладки вообще. На не-Windows keycode другой,
# поэтому там остаётся разбор по keysym — в Linux/macOS сборках этой
# программы всё равно нет, но и ломать их незачем.
#
# Дополнительно поддержаны классические Windows-сочетания Ctrl+Insert /
# Shift+Insert / Shift+Delete: их ждут пользователи старой школы, а стоят
# они три строки.
#
# Каждый обработчик возвращает "break": без этого Tk после нашего кода
# прогоняет ещё и встроенный биндинг того же события — текст вставлялся
# дважды, а Ctrl+A в поле ввода вместо выделения прыгал курсором в начало
# строки (стандартное emacs-поведение tk.Entry).
# ---------------------------------------------------------------------------
_WIN_VK_ACTIONS = {67: "copy", 86: "paste", 88: "cut", 65: "select_all"}
_KEYSYM_ACTIONS = {"c": "copy", "v": "paste", "x": "cut", "a": "select_all"}


def _hotkey_action_name(event) -> str | None:
    """Какое действие запрошено сочетанием с Ctrl — независимо от раскладки."""
    if os.name == "nt":
        name = _WIN_VK_ACTIONS.get(getattr(event, "keycode", None))
        if name is not None:
            return name
    return _KEYSYM_ACTIONS.get((getattr(event, "keysym", "") or "").lower())


def attach_hotkeys(widget):
    """Привязывает Ctrl+C/V/X/A (и Ctrl+Insert/Shift+Insert/Shift+Delete)."""
    actions = {
        "copy": _do_copy,
        "paste": _do_paste,
        "cut": _do_cut,
        "select_all": _do_select_all,
    }

    def on_ctrl_key(event, w=widget):
        name = _hotkey_action_name(event)
        if name is None:
            return None  # прочие сочетания с Ctrl отдаём Tk как есть
        actions[name](w)
        return "break"

    widget.bind("<Control-KeyPress>", on_ctrl_key, add="+")

    def _bind_simple(seq, fn):
        def handler(event, w=widget, f=fn):
            f(w)
            return "break"
        widget.bind(seq, handler, add="+")

    _bind_simple("<Control-Insert>", _do_copy)
    _bind_simple("<Shift-Insert>", _do_paste)
    _bind_simple("<Shift-Delete>", _do_cut)


def attach_context_menu(widget):
    menu = tk.Menu(widget, tearoff=0)
    menu.add_command(label="Вырезать", accelerator="Ctrl+X",
                     command=lambda: _do_cut(widget))
    menu.add_command(label="Копировать", accelerator="Ctrl+C",
                     command=lambda: _do_copy(widget))
    menu.add_command(label="Вставить", accelerator="Ctrl+V",
                     command=lambda: _do_paste(widget))
    menu.add_separator()
    menu.add_command(label="Выделить всё", accelerator="Ctrl+A",
                     command=lambda: _do_select_all(widget))

    def show_menu(event):
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()

    widget.bind("<Button-3>", show_menu, add="+")
    if os.name == "nt":
        # <App> — клавиша "контекстное меню" рядом с правым Ctrl. Такого
        # keysym нет в X11-сборках Tk, там bind() на него бросает TclError.
        widget.bind("<App>", show_menu, add="+")
    return menu


def attach_input_features(widget):
    attach_hotkeys(widget)
    attach_context_menu(widget)


# ---------------------------------------------------------------------------
# Главное окно
# ---------------------------------------------------------------------------
class App(ctk.CTk):
    def __init__(self):
        # ВАЖЕН ПОРЯДОК: тему ставим ДО создания окна.
        #
        # ctk.set_appearance_mode() на Windows перекрашивает заголовок окна
        # через DwmSetWindowAttribute, а чтобы новый цвет отрисовался,
        # customtkinter прячет и заново показывает окно (withdraw + update,
        # см. CTk._windows_set_titlebar_color). Этот приём сбрасывает
        # состояние "развёрнуто на весь экран" обратно в обычное окно —
        # именно поэтому раньше приложение на мгновение открывалось во весь
        # экран и тут же схлопывалось до 1000x750. Пока окна ещё нет,
        # перекрашивать нечего, и никакого withdraw не происходит.
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        super().__init__()

        self.title(f"HelixFillDNA  v{__version__}")

        # Иконка окна (слева от заголовка и на панели задач) вместо
        # стандартной заглушки Tk. iconbitmap() применяется дважды
        # намеренно: customtkinter внутри пересоздаёт часть оформления
        # окна уже после __init__, и на некоторых версиях Windows/Tk
        # первый вызов при этом теряется — повторный через after()
        # ставит иконку окончательно. Оба вызова защищены: отсутствие
        # или повреждение .ico не должно мешать запуску приложения.
        self._apply_window_icon()
        self.after(200, self._apply_window_icon)
        # Стартовый размер сразу задаём по экрану, а не 1000x750: иначе на
        # системах, где state("zoomed") срабатывает не мгновенно (или не
        # срабатывает вовсе), пользователь успевает увидеть, как маленькое
        # окно прыгает в полный экран. Минус в том, что «свернуть в окно»
        # вернёт почти тот же размер — это меньшее зло, чем прыжок при
        # каждом запуске.
        try:
            self.geometry(
                f"{self.winfo_screenwidth()}x{max(400, self.winfo_screenheight() - 70)}+0+0"
            )
        except tk.TclError:
            self.geometry("1000x750")
        self.minsize(850, 650)
        # Открываемся развёрнутыми на весь экран: на вкладке "Запуск" три
        # шага, шкала прогресса и карта из 22 хромосом — в окне 1000x750
        # это сплошная прокрутка.
        #
        # after(0) — а не прямой вызов: customtkinter при первом показе
        # окна ещё раз прячет и показывает его (CTk.mainloop() -> withdraw
        # + deiconify), и разворот, сделанный до этого, потерялся бы.
        # Обработчик after выполняется уже после этой процедуры.
        self.after(0, self._maximize_window)

        self.log_q: queue.Queue = queue.Queue()
        self.running = False
        self._pwd_visible = False
        # Задача 1: сигнал отмены для download_donors.download_donors_for_chip().
        # threading.Event, а не bool-флаг — так cancel_check() из фонового
        # потока и .set()/.clear() из главного потока не гонятся друг с
        # другом без явной блокировки.
        self._cancel_donor_download = threading.Event()

        # Промт "Именованные папки запуска": какой конкретно запуск
        # (output/runs/<run_name>/) сейчас активен — выставляется либо
        # при старте Этапа 1-6 (_on_start), либо при выборе существующего
        # запуска из истории ("Продолжить (Этап 2)"). Этап 7 (_on_mis)
        # без этого не знает, куда класть rerun_results/собранный файл.
        self.current_run_dir: Path | None = None
        self.current_run_name: str | None = None
        self._run_log_handler = None  # logging.Handler, снимается при смене запуска
        self._run_history_map: dict[str, Path] = {}  # подпись в списке -> папка запуска
        # Ключи ("chr14", "chr14_idx", ...), для которых в log_text уже
        # создана обновляемая строка прогресса (см. _upsert_progress_line).
        self._progress_keys: set[str] = set()

        # Промт "сделать вкладку Запуск юзерфрендли": состояние мастера.
        # _run_started_at — момент старта длинной операции (time.monotonic),
        # по нему считается оценка оставшегося времени под шкалой; None,
        # когда ничего не выполняется.
        self._wizard_step = 1
        self._run_started_at: float | None = None
        self._run_details_visible = False
        # Карта состояний донорских хромосом этапа 3: {номер: (ранг, текст)}.
        self._donor_states: dict[int, tuple[int, str]] = {}
        # Наблюдение за папками доноров: какие папки опрашивать, последние
        # увиденные размеры файлов и id запланированного after()-вызова.
        self._donor_watch_dirs: list[Path] = []
        # Что именно сейчас наблюдаем: "donors" (этап 3 Шага 1) или "mis"
        # (скачивание результатов на Шаге 3). От этого зависит, в какую
        # панель писать — сам механизм опроса общий.
        self._watch_target = "donors"
        self._donor_file_sizes: dict[Path, tuple[float, int]] = {}
        self._donor_watch_id = None
        # Путь к последнему собранному итоговому файлу (в results/).
        self._last_result_path: Path | None = None

        # Промт "обычные / продвинутые настройки": папка с трафаретами
        # создаётся программой сама при первом запуске, чтобы в обычном
        # режиме было куда класть template_v3.txt / template_v5.txt.
        _samples_dir()
        _results_dir()

        # Нижняя полоска: версия слева, обратная связь справа. Пакуется
        # ДО tabview и с side="bottom" — иначе tabview с expand=True забрал
        # бы всё место и полоску прижало бы в ноль.
        footer = ctk.CTkFrame(self, fg_color="transparent")
        footer.pack(side="bottom", fill="x", padx=12, pady=(0, 8))
        ctk.CTkLabel(
            footer, text=f"HelixFillDNA v{__version__}", text_color="gray50",
        ).pack(side="left")
        ctk.CTkButton(
            footer, text="✉ Сообщить об ошибке / предложить улучшение", width=340,
            fg_color="transparent", border_width=1,
            command=self._open_feedback_dialog,
        ).pack(side="right")

        self.tabview = ctk.CTkTabview(self)
        self.tabview.pack(fill="both", expand=True, padx=10, pady=10)

        self.tab_settings = self.tabview.add("Подготовка")
        self.tab_run = self.tabview.add("Запуск")
        self.tab_log = self.tabview.add("Лог")

        self._build_settings_tab()
        self._build_run_tab()
        self._build_log_tab()

        self._poll_logs()

    # -----------------------------------------------------------------------
    # Вкладка "Подготовка"
    # -----------------------------------------------------------------------
    def _maximize_window(self):
        """
        Разворачивает окно на весь экран. state("zoomed") работает на
        Windows и macOS; в X11-сборках Tk такого состояния нет — там
        отдельный атрибут "-zoomed", а если и его нет (некоторые оконные
        менеджеры), просто растягиваем окно по размеру экрана. Любая
        неудача не должна мешать запуску: окно останется 1000x750.
        """
        self.after(300, self._verify_maximized)
        try:
            self.state("zoomed")
            return
        except tk.TclError:
            pass
        try:
            self.attributes("-zoomed", True)
            return
        except tk.TclError:
            pass
        self._stretch_to_screen()

    def _stretch_to_screen(self):
        try:
            self.geometry(f"{self.winfo_screenwidth()}x{self.winfo_screenheight()}+0+0")
        except tk.TclError:
            pass

    def _verify_maximized(self):
        """
        state("zoomed") и атрибут "-zoomed" не бросают ошибку, когда просто
        ничего не делают (например, оконный менеджер их не поддерживает) —
        поэтому недостаточно вызвать их и понадеяться. Через треть секунды
        после старта смотрим на фактический размер окна и, если оно так и
        осталось маленьким, растягиваем его руками.
        """
        try:
            too_narrow = self.winfo_width() < self.winfo_screenwidth() * 0.9
            too_short = self.winfo_height() < self.winfo_screenheight() * 0.8
        except tk.TclError:
            return
        if too_narrow or too_short:
            self._stretch_to_screen()

    def _apply_window_icon(self):
        """Ставит app_icon.ico на окно. Молча ничего не делает, если файла
        нет или Tk отказался его читать (не-Windows платформа, битый .ico)."""
        icon = _find_app_icon()
        if icon is None:
            return
        try:
            self.iconbitmap(str(icon))
        except tk.TclError:
            pass

    def _build_settings_tab(self):
        scroll = ctk.CTkScrollableFrame(self.tab_settings)
        scroll.pack(fill="both", expand=True, padx=10, pady=10)

        # --- Промт "обычные / продвинутые настройки": переключатель ------
        ctk.CTkLabel(
            scroll, text="Основные настройки",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 8))
        self.mode_switch = ctk.CTkSegmentedButton(
            scroll, values=[MODE_SIMPLE, MODE_ADVANCED],
            command=self._on_mode_changed, width=340,
        )
        self.mode_switch.set(_load_ui_mode())
        self.mode_switch.pack(anchor="w", pady=(0, 15))

        # ==================================================================
        # Базовый блок — виден в ОБОИХ режимах.
        # ==================================================================
        base = ctk.CTkFrame(scroll, fg_color="transparent")
        base.pack(fill="x")
        self.basic_box = base

        ctk.CTkLabel(base, text="Источник данных:").pack(anchor="w")
        source_names = [v["name"] for v in pipeline.SOURCES.values()]
        self.source_dd = ctk.CTkOptionMenu(
            base, values=source_names, width=400,
            command=lambda _choice: self._on_source_changed(),
        )
        self.source_dd.set(source_names[0])
        self.source_dd.pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(base, text="Файл с данными:").pack(anchor="w")
        row1 = ctk.CTkFrame(base, fg_color="transparent")
        row1.pack(fill="x", pady=(0, 15))
        self.input_tf = ctk.CTkEntry(row1, placeholder_text="Выберите файл...", width=600)
        self.input_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(row1, text="Обзор", width=100,
                      command=lambda: self._pick_file(self.input_tf)).pack(side="right")

        # Сводка того, что обычный режим выбрал за пользователя. Пакуется
        # и распаковывается в _on_mode_changed() вместе с self.adv_box.
        self.simple_info_lbl = ctk.CTkLabel(
            base, text="", justify="left", text_color="gray60", wraplength=700,
        )
        self.simple_tmpl_lbl = ctk.CTkLabel(
            base, text="", justify="left", text_color="#4CAF50", wraplength=700,
        )

        # ==================================================================
        # Продвинутый блок — скрывается целиком в обычном режиме.
        # ==================================================================
        adv = ctk.CTkFrame(scroll, fg_color="transparent")
        adv.pack(fill="x")
        self.adv_box = adv

        # --- Шаг 1 промта "HRC / TopMed": выбор референсной панели -------
        # В обычном режиме панель не показывается и всегда берётся
        # pipeline.DEFAULT_PANEL (HRC): смена панели тянет за собой другую
        # сборку генома, отдельный референс, отдельный кэш доноров и
        # лифтовер координат — это осознанный выбор, а не рутинная настройка.
        ctk.CTkLabel(adv, text="Референсная панель импутации:").pack(anchor="w")
        panel_names = [v["display_name"] for v in pipeline.REFERENCE_PANELS.values()]
        self.panel_dd = ctk.CTkOptionMenu(
            adv, values=panel_names, width=400,
            command=lambda _choice: self._on_panel_changed(),
        )
        self.panel_dd.set(pipeline.REFERENCE_PANELS[pipeline.DEFAULT_PANEL]["display_name"])
        self.panel_dd.pack(anchor="w", pady=(0, 5))

        self.panel_warning_lbl = ctk.CTkLabel(
            adv, text="", justify="left", text_color="#F9A825", wraplength=700,
        )
        self.panel_warning_lbl.pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(
            adv,
            text=("ℹ Референсный геном для FTDNA/MyHeritage проверяется\n"
                  "автоматически при запуске под выбранную выше панель. Если\n"
                  "файла нет в reference/<панель>/ — он будет скачан и распакован\n"
                  "(размер зависит от сборки генома выбранной панели)."),
            justify="left", text_color="gray60",
        ).pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(
            adv,
            text=("Трафарет (по умолчанию берётся из samples/, здесь можно "
                  "указать свой файл):"),
        ).pack(anchor="w")
        row3 = ctk.CTkFrame(adv, fg_color="transparent")
        row3.pack(fill="x", pady=(0, 15))
        self.tmpl_tf = ctk.CTkEntry(row3, placeholder_text="Выберите файл...", width=600)
        self.tmpl_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(row3, text="Обзор", width=100,
                      command=lambda: self._pick_file(self.tmpl_tf)).pack(side="right")

        ctk.CTkLabel(adv, text="Папка с бинарниками htslib:").pack(anchor="w")
        row_bin = ctk.CTkFrame(adv, fg_color="transparent")
        row_bin.pack(fill="x", pady=(0, 5))
        self.bin_tf = ctk.CTkEntry(row_bin, width=600)
        self.bin_tf.insert(0, str(_detect_bin_dir()))
        self.bin_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(
            row_bin, text="🌐 Диагностика сети", width=170,
            command=self._on_diagnose_network,
        ).pack(side="right")
        ctk.CTkLabel(
            adv,
            text=("ℹ Проверяет CA-сертификаты для libcurl (bcftools) и наличие "
                  "конфликтующего curl.exe в папке бинарников — нужно только "
                  "для ускоренного удалённого скачивания доноров (Этап 3), "
                  "само скачивание работает и без этого через обычный "
                  "полный путь."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 15))

        ctk.CTkFrame(adv, height=2, fg_color="gray40").pack(fill="x", pady=15)

        ctk.CTkLabel(
            adv, text="Параметры вывода",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(adv, text="Формат вывода:").pack(anchor="w")
        self.format_dd = ctk.CTkOptionMenu(
            adv,
            values=["v3 (LF, ~97% call rate)", "v5 (CRLF, ~92% call rate)"],
            width=400,
        )
        self.format_dd.set("v3 (LF, ~97% call rate)")
        self.format_dd.pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(
            adv, text="Порог Rsq (качество импутации, от 0 до 1):",
            font=ctk.CTkFont(weight="bold"),
        ).pack(anchor="w")
        rsq_info = (
            "Чем выше Rsq — тем надёжнее генотип, но тем меньше позиций проходит фильтр.\n"
            "0.30 — стандартный порог MIS: максимум позиций, часть — с невысоким качеством.\n"
            "0.80 — баланс количества и качества.\n"
            "0.90 — только высокоточные варианты (меньше позиций, но надёжнее).\n"
            "0.95+ — максимальное качество для критичных задач."
        )
        ctk.CTkLabel(adv, text=rsq_info, justify="left", text_color="gray60").pack(
            anchor="w", pady=(4, 8)
        )

        row_rsq = ctk.CTkFrame(adv, fg_color="transparent")
        row_rsq.pack(fill="x", pady=(0, 5))
        ctk.CTkLabel(row_rsq, text="Значение:").pack(side="left", padx=(0, 10))
        self.rsq_entry = ctk.CTkEntry(row_rsq, width=100, placeholder_text="0.30")
        self.rsq_entry.insert(0, "0.30")
        self.rsq_entry.pack(side="left")
        self.rsq_entry.bind("<KeyRelease>", lambda e: self._validate_rsq_entry())

        self.rsq_status_lbl = ctk.CTkLabel(
            adv, text="✓ Порог принят: 0.30", text_color="#4CAF50",
        )
        self.rsq_status_lbl.pack(anchor="w", pady=(0, 15))

        # Задача C: опциональная нормализация multiallelic-сайтов
        # (bcftools norm -m-both). НЕ входит в критический путь фикса
        # Invalid alleles (это делают Задачи A/B) — отдельная оптимизация,
        # поэтому выключена по умолчанию в продвинутом режиме, как и в CLI
        # (--normalize); обычный режим включает её сам.
        self.normalize_var = ctk.BooleanVar(value=False)
        self.normalize_cb = ctk.CTkCheckBox(
            adv,
            text="Нормализовать multiallelic-сайты перед split (bcftools norm -m-both, опционально)",
            variable=self.normalize_var,
        )
        self.normalize_cb.pack(anchor="w", pady=(0, 15))

        # Задача D: опциональное переиспользование доноров между разными
        # людьми на одном чипе (широкая сигнатура вместо строгой).
        self.reuse_donors_var = ctk.BooleanVar(value=False)
        self.reuse_donors_cb = ctk.CTkCheckBox(
            adv,
            text=("Переиспользовать доноров между разными людьми на одном чипе "
                  "(экспериментально, Задача D)"),
            variable=self.reuse_donors_var,
        )
        self.reuse_donors_cb.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            adv,
            text=("ℹ Строит сигнатуру чипа по всем измеренным позициям вместо позиций, "
                  "прошедших личный QC. Позволяет не перекачивать ~22 файла доноров "
                  "для каждого нового человека на том же чипе. Требует пересборки "
                  "batch_merged без принудительной подстановки 0/0 на пропусках — "
                  "это уже обеспечено в этой версии."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 15))

        # Промт "Monomorphic sites / настраиваемое количество EUR-доноров":
        # раньше число доноров-образцов 1000 Genomes было жёстко зашито
        # как 20 (download_donors.create_eur20_list()). На такой маленькой
        # случайной выборке многие сайты, полиморфные в популяции в целом,
        # случайно оказываются мономорфными во всех 20 взятых образцах —
        # MIS отбрасывает такие сайты из QC ("Monomorphic sites"), снижая
        # итоговое покрытие. По умолчанию в продвинутом режиме используется
        # ВСЯ доступная EUR-подвыборка панели (обычно порядка 500 человек);
        # обычный режим сознательно берёт 20 — ради разумного трафика и
        # времени первого запуска.
        self.eur_all_var = ctk.BooleanVar(value=True)
        self.eur_all_cb = ctk.CTkCheckBox(
            adv,
            text=("Использовать всех доступных EUR-доноров 1000 Genomes "
                  "(уменьшает Monomorphic sites на QC MIS, но увеличивает "
                  "трафик/время скачивания доноров)"),
            variable=self.eur_all_var,
            command=self._on_eur_all_toggled,
        )
        self.eur_all_cb.pack(anchor="w", pady=(0, 5))

        row_eur_count = ctk.CTkFrame(adv, fg_color="transparent")
        row_eur_count.pack(fill="x", pady=(0, 5))
        ctk.CTkLabel(row_eur_count, text="Или конкретное число доноров:").pack(
            side="left", padx=(0, 10)
        )
        self.eur_count_entry = ctk.CTkEntry(
            row_eur_count, width=100, placeholder_text="напр. 20",
            state="disabled",
        )
        self.eur_count_entry.pack(side="left")
        self.eur_count_entry.bind("<KeyRelease>", lambda e: self._validate_eur_count_entry())
        attach_input_features(self.eur_count_entry)

        self.eur_count_status_lbl = ctk.CTkLabel(
            adv, text="✓ Будут использованы все доступные EUR-доноры (~500)",
            text_color="#4CAF50",
        )
        self.eur_count_status_lbl.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            adv,
            text=("ℹ По умолчанию (галочка включена) используются ВСЕ EUR-образцы "
                  "из панели 1000 Genomes — это снижает долю Monomorphic sites, "
                  "которые MIS исключает из QC на маленькой (20 образцов) случайной "
                  "выборке. Снимите галочку и укажите число (например, 20), если "
                  "важнее скорость/трафик, а не полнота покрытия. Изменение этого "
                  "значения между запусками требует перекачки доноров — старый кэш "
                  "с другим числом образцов не подходит и не будет использован "
                  "автоматически."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 15))

        # Промт "Доноры для VCF-источника: понятная отмена + общий кэш
        # сырых хромосом", Шаг 5: общий кэш ЕЩЁ НЕ отфильтрованных полных
        # хромосом 1000 Genomes, переиспользуемый МЕЖДУ ВСЕМИ источниками
        # (ftdna/myheritage/vcf) и чипами ОДНОЙ референсной сборки —
        # экономит трафик при повторных запусках/переключении источника
        # ценой постоянного места на диске. По умолчанию выключен в
        # продвинутом режиме; обычный режим включает его.
        self.raw_cache_var = ctk.BooleanVar(value=False)
        self.raw_cache_cb = ctk.CTkCheckBox(
            adv,
            text=("Хранить сырые (нефильтрованные) хромосомы 1000 Genomes для "
                  "повторного использования между разными источниками/чипами "
                  "(~десятки ГБ на диске, экономит трафик при последующих запусках)"),
            variable=self.raw_cache_var,
        )
        self.raw_cache_cb.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            adv,
            text=("ℹ Доноры каждого источника/чипа всё равно хранятся отдельно "
                  "(donors/<source>/<panel>/) — этот кэш касается только ПОЛНЫХ, "
                  "ещё не отфильтрованных хромосом 1000 Genomes "
                  "(donors/_raw_chromosomes/<сборка_генома>/), которые при "
                  "полном скачивании (без удалённой фильтрации) одинаковы для "
                  "любого источника. Занимает дополнительно ~13-20 ГБ на диске "
                  "(все 22 хромосомы 1000 Genomes phase3). Можно включить позже "
                  "и для повторного запуска — уже скачанные хромосомы других "
                  "источников/чипов будут переиспользованы сразу."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 15))

        for entry in (self.input_tf, self.tmpl_tf, self.bin_tf, self.rsq_entry):
            attach_input_features(entry)

        # Первичная синхронизация предупреждения под панель по умолчанию.
        self._on_panel_changed()
        # Применяем режим, восстановленный из ui_state.json: в обычном
        # режиме это скроет self.adv_box и проставит автоматические
        # значения (в том числе трафарет из samples/).
        self._on_mode_changed(self.mode_switch.get())

    # -----------------------------------------------------------------------
    # Вкладка "Запуск"
    # -----------------------------------------------------------------------
    def _build_run_tab(self):
        # ==================================================================
        # Вкладка "Запуск" — три подвкладки по числу шагов.
        #
        # Раньше всё содержимое трёх шагов лежало одним длинным скроллом, и
        # пользователь одновременно видел кнопку запуска подготовки, поля
        # для письма MIS и инструкцию для сайта — хотя в каждый момент
        # времени осмысленно ровно одно из трёх. Теперь каждый шаг живёт в
        # своей подвкладке, а программа сама переключает её по мере
        # выполнения (_set_wizard_step). Переключиться руками тоже можно —
        # посмотреть вперёд или вернуться никто не мешает.
        #
        # Общее для всех шагов (выбор запуска и шкала прогресса) остаётся
        # НАД подвкладками: прогресс должен быть виден, на какой бы шаг ни
        # переключился пользователь.
        # ==================================================================
        header = ctk.CTkFrame(self.tab_run, fg_color="transparent")
        header.pack(fill="x", padx=10, pady=(10, 0))

        run_row = ctk.CTkFrame(header, fg_color="transparent")
        run_row.pack(fill="x", pady=(0, 5))
        ctk.CTkLabel(run_row, text="Запуск:").pack(side="left", padx=(0, 10))
        self.run_history_dd = ctk.CTkOptionMenu(run_row, values=["(нет запусков)"], width=340)
        self.run_history_dd.pack(side="left", padx=(0, 10))
        ctk.CTkButton(
            run_row, text="＋ Новый", width=105, command=self._on_new_run,
        ).pack(side="left", padx=(0, 5))
        ctk.CTkButton(
            run_row, text="📂 Папка запуска", width=150, command=self._on_open_run_folder,
        ).pack(side="left", padx=(0, 5))
        self.run_details_btn = ctk.CTkButton(
            run_row, text="▾ Подробнее", width=125,
            fg_color="transparent", border_width=1,
            command=self._toggle_run_details,
        )
        self.run_details_btn.pack(side="left")

        self.active_run_lbl = ctk.CTkLabel(header, text="Активный запуск: нет", text_color="gray60")
        self.active_run_lbl.pack(anchor="w", pady=(0, 5))

        # --- Свёрнутые подробности запуска -------------------------------
        self.run_details_box = ctk.CTkFrame(header, fg_color="transparent")
        self._run_details_visible = False

        ctk.CTkLabel(self.run_details_box, text="Название нового запуска:").pack(anchor="w")
        run_name_row = ctk.CTkFrame(self.run_details_box, fg_color="transparent")
        run_name_row.pack(fill="x", pady=(0, 5))
        self.run_name_tf = ctk.CTkEntry(run_name_row, width=200)
        self.run_name_tf.pack(side="left", padx=(0, 10))
        ctk.CTkButton(
            run_name_row, text="▶ Продолжить (Шаг 3)", width=190,
            command=self._on_continue_run,
        ).pack(side="left", padx=(0, 5))
        ctk.CTkButton(
            run_name_row, text="✏ Переименовать", width=160,
            command=self._on_rename_run,
        ).pack(side="left", padx=(0, 5))
        ctk.CTkButton(
            run_name_row, text="⟳ Обновить историю", width=180,
            command=self._refresh_run_history,
        ).pack(side="left")
        attach_input_features(self.run_name_tf)

        ctk.CTkLabel(
            self.run_details_box,
            text=("ℹ Каждый запуск пишет рабочие файлы в свою папку "
                  "output/runs/<название>/ — донор-кэш (donors/) общий для "
                  "всех запусков и не дублируется. Итоговый файл кладётся "
                  f"отдельно, в папку {RESULTS_DIR_NAME}/ рядом с программой.\n"
                  "«Продолжить» нужен, если письмо от MIS пришло уже после "
                  "перезапуска программы: выберите свой запуск в списке "
                  "выше и сделайте его активным, не прогоняя Шаг 1 заново."),
            justify="left", text_color="gray60", wraplength=760,
        ).pack(anchor="w", pady=(0, 5))

        # --- Общая шкала прогресса ---------------------------------------
        self.stage_lbl = ctk.CTkLabel(header, text="Готов к запуску",
                                      font=ctk.CTkFont(size=16, weight="bold"))
        self.stage_lbl.pack(anchor="w", pady=(8, 3))

        self.progress = ctk.CTkProgressBar(header, height=15)
        self.progress.pack(fill="x", pady=(0, 3))
        self.progress.set(0)

        # Оценка времени: скачивание доноров идёт часами, и голая полоса
        # без времени читается как "зависло".
        self.eta_lbl = ctk.CTkLabel(header, text="", text_color="gray60")
        self.eta_lbl.pack(anchor="w")

        # Последняя строка лога прямо здесь — чтобы при ошибке не нужно
        # было догадываться переключиться на вкладку "Лог".
        log_peek_row = ctk.CTkFrame(header, fg_color="transparent")
        log_peek_row.pack(fill="x", pady=(0, 5))
        self.last_log_lbl = ctk.CTkLabel(
            log_peek_row, text="", text_color="gray60", anchor="w", justify="left",
        )
        self.last_log_lbl.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(
            log_peek_row, text="Показать лог", width=130,
            fg_color="transparent", border_width=1,
            command=lambda: self.tabview.set("Лог"),
        ).pack(side="right")

        # --- Подвкладки шагов --------------------------------------------
        self.run_tabs = ctk.CTkTabview(self.tab_run)
        self.run_tabs.pack(fill="both", expand=True, padx=10, pady=(5, 10))
        # Текущие подписи подвкладок: rename() меняет ключ, по которому
        # работает .set(), поэтому актуальные имена держим здесь.
        self._run_tab_names = list(RUN_TAB_BASE_NAMES)
        for name in self._run_tab_names:
            self.run_tabs.add(name)

        self._build_step1_tab(self.run_tabs.tab(self._run_tab_names[0]))
        self._build_step2_tab(self.run_tabs.tab(self._run_tab_names[1]))
        self._build_step3_tab(self.run_tabs.tab(self._run_tab_names[2]))

        self._refresh_run_instructions()
        self._refresh_run_name_suggestion()
        self._refresh_run_history()
        self._set_wizard_step(1)
        self._refresh_mis_btn_state()

    # --- Шаг 1: подготовка файлов (этапы 1-6) ----------------------------
    def _build_step1_tab(self, parent):
        scroll = ctk.CTkScrollableFrame(parent)
        scroll.pack(fill="both", expand=True)

        ctk.CTkLabel(
            scroll,
            text=("Программа прочитает ваш файл, скачает донорские хромосомы "
                  "1000 Genomes и подготовит 22 файла для загрузки на сервер "
                  "импутации. Самая долгая часть — скачивание доноров (этап 3), "
                  "оно может идти несколько часов."),
            justify="left", text_color="gray60", wraplength=760,
        ).pack(anchor="w", pady=(0, 10))

        self.start_btn = ctk.CTkButton(
            scroll, text="Запустить подготовку (этапы 1-6)",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color="#2E7D32", hover_color="#1B5E20",
            height=50, corner_radius=12,
            command=self._on_start,
        )
        self.start_btn.pack(fill="x", pady=(0, 10))

        # Кнопка остановки живёт в собственном контейнере: её пакуют и
        # распаковывают (_show_cancel_donor_btn/_hide_cancel_donor_btn), а
        # не просто гасят — неактивная красная кнопка, которая почти всё
        # время бесполезна, только сбивала с толку. Реально прерывает
        # текущий subprocess (curl), а не ждёт окончания хромосомы.
        self.stop_box = ctk.CTkFrame(scroll, fg_color="transparent")
        self.cancel_donor_btn = ctk.CTkButton(
            self.stop_box, text="⏹ Остановить скачивание доноров",
            fg_color="#B71C1C", hover_color="#7F0000",
            command=self._on_cancel_donor_download,
        )
        self.cancel_donor_btn.pack(fill="x")

        # --- Живая карта 22 хромосом (этап 3) ----------------------------
        # Скачивание идёт в несколько потоков и часами; общий процент
        # «Обработано хромосом: 7/22» не отвечает на вопрос «оно вообще
        # шевелится?». Здесь по каждой хромосоме видно, что именно с ней
        # происходит прямо сейчас: качается (с процентами), фильтруется,
        # готова или упала. Панель появляется, только когда пошли
        # сообщения про хромосомы, и прячется в начале нового запуска.
        self.donor_panel = ctk.CTkFrame(scroll, border_width=1, border_color="gray40")
        ctk.CTkLabel(
            self.donor_panel, text="Этап 3 · Донорские хромосомы 1000 Genomes",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).grid(row=0, column=0, columnspan=DONOR_GRID_COLUMNS, sticky="w",
               padx=10, pady=(8, 6))

        self.donor_chr_lbls: dict[int, ctk.CTkLabel] = {}
        for chrom in range(1, 23):
            idx = chrom - 1
            row = 1 + idx // DONOR_GRID_COLUMNS
            col = idx % DONOR_GRID_COLUMNS
            lbl = ctk.CTkLabel(
                self.donor_panel, text=f"chr{chrom} —", anchor="w",
                text_color="gray50", width=150,
            )
            lbl.grid(row=row, column=col, sticky="w", padx=10, pady=2)
            self.donor_chr_lbls[chrom] = lbl
        for col in range(DONOR_GRID_COLUMNS):
            self.donor_panel.grid_columnconfigure(col, weight=1)

        base_row = 1 + (21 // DONOR_GRID_COLUMNS) + 1
        # Что качается прямо сейчас: имя файла, сколько мегабайт уже на
        # диске и с какой скоростью растёт. Заполняется _poll_donor_files().
        self.donor_active_lbl = ctk.CTkLabel(
            self.donor_panel, text="", text_color="#42A5F5", justify="left",
        )
        self.donor_active_lbl.grid(
            row=base_row, column=0, columnspan=DONOR_GRID_COLUMNS,
            sticky="w", padx=10, pady=(8, 2),
        )

        self.donor_summary_lbl = ctk.CTkLabel(
            self.donor_panel, text="", text_color="gray60", justify="left",
        )
        self.donor_summary_lbl.grid(
            row=base_row + 1, column=0,
            columnspan=DONOR_GRID_COLUMNS, sticky="w", padx=10, pady=(2, 10),
        )

    # --- Шаг 2: ручная работа на сайте MIS -------------------------------
    def _build_step2_tab(self, parent):
        scroll = ctk.CTkScrollableFrame(parent)
        scroll.pack(fill="both", expand=True)

        # Текст держим коротким: кнопки под ним открывают сайт и папку, а
        # значения формы вынесены в отдельный заметный блок ниже — в
        # абзаце они терялись. Формируется динамически из выбранной
        # панели (_refresh_run_instructions).
        self.run_instructions_lbl = ctk.CTkLabel(scroll, text="", justify="left")
        self.run_instructions_lbl.pack(anchor="w", pady=(0, 10))

        mis_actions = ctk.CTkFrame(scroll, fg_color="transparent")
        mis_actions.pack(fill="x", pady=(0, 5))
        ctk.CTkButton(
            mis_actions, text="🌐 Открыть сервер импутации", width=240,
            command=self._on_open_mis_site,
        ).pack(side="left", padx=(0, 8))
        ctk.CTkButton(
            mis_actions, text="📂 Папка с 22 файлами", width=210,
            command=self._on_open_upload_folder,
        ).pack(side="left")

        self.mis_actions_status_lbl = ctk.CTkLabel(scroll, text="", text_color="gray60")
        self.mis_actions_status_lbl.pack(anchor="w", pady=(0, 8))

        # Параметры формы MIS отдельным заметным блоком: ошибиться здесь
        # дорого (неверный Array Build — провал QC и потерянные часы).
        # Кнопки "скопировать" тут нет намеренно: на сайте это выпадающие
        # списки, вставлять в них нечего — значения нужно ВИДЕТЬ.
        params_box = ctk.CTkFrame(scroll, border_width=1, border_color="#42A5F5")
        params_box.pack(fill="x", pady=(0, 15))
        ctk.CTkLabel(
            params_box, text="Выберите в форме на сайте именно эти значения:",
            text_color="gray70",
        ).pack(anchor="w", padx=12, pady=(10, 4))
        self.mis_params_lbl = ctk.CTkLabel(
            params_box, text="", justify="left",
            font=ctk.CTkFont(size=15, weight="bold"),
        )
        self.mis_params_lbl.pack(anchor="w", padx=12, pady=(0, 4))
        ctk.CTkLabel(
            params_box,
            text="Остальные поля формы оставьте со значениями по умолчанию.",
            text_color="gray60",
        ).pack(anchor="w", padx=12, pady=(0, 10))

        ctk.CTkLabel(
            scroll,
            text=("Когда сервер закончит импутацию, он пришлёт письмо со "
                  "ссылкой и паролем — с ними переходите на Шаг 3. Программу "
                  "можно закрыть и вернуться позже: выберите свой запуск в "
                  "списке сверху и нажмите «▶ Продолжить (Шаг 3)» в подробностях."),
            justify="left", text_color="gray60", wraplength=760,
        ).pack(anchor="w")

    # --- Шаг 3: скачивание результатов MIS и сборка (этап 7) -------------
    def _build_step3_tab(self, parent):
        scroll = ctk.CTkScrollableFrame(parent)
        scroll.pack(fill="both", expand=True)

        ctk.CTkLabel(scroll, text="curl-команда из письма MIS:").pack(anchor="w")
        self.curl_tf = ctk.CTkTextbox(scroll, height=80, width=700)
        self.curl_tf.pack(fill="x", pady=(0, 3))
        self.curl_tf.bind("<KeyRelease>", lambda e: self._validate_curl_entry())
        self.curl_tf.bind("<<Paste>>", lambda e: self.after(50, self._validate_curl_entry))
        self.curl_tf.bind("<FocusOut>", lambda e: self._validate_curl_entry())

        self.curl_status_lbl = ctk.CTkLabel(
            scroll, text="", text_color="gray60", justify="left", wraplength=760,
        )
        self.curl_status_lbl.pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="Пароль из письма:").pack(anchor="w")
        pwd_row = ctk.CTkFrame(scroll, fg_color="transparent")
        pwd_row.pack(fill="x", pady=(0, 15))
        self.pwd_tf = ctk.CTkEntry(pwd_row, show="*", width=520, placeholder_text="Пароль")
        self.pwd_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        self.pwd_tf.bind("<KeyRelease>", lambda e: self._refresh_mis_btn_state())
        self.pwd_tf.bind("<<Paste>>", lambda e: self.after(50, self._refresh_mis_btn_state))
        self.show_pwd_btn = ctk.CTkButton(
            pwd_row, text="👁", width=50,
            command=self._toggle_pwd,
        )
        self.show_pwd_btn.pack(side="right", padx=(5, 0))
        ctk.CTkButton(
            pwd_row, text="🔍 Диагностика", width=140,
            command=self._on_diagnose_password,
        ).pack(side="right")

        self.mis_btn = ctk.CTkButton(
            scroll, text="📥 Скачать результаты и собрать финальный файл",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color="#1565C0", hover_color="#0D47A1",
            height=50, corner_radius=12,
            state="disabled",
            command=self._on_mis,
        )
        self.mis_btn.pack(fill="x", pady=(0, 3))

        # Почему кнопка серая — раньше это было неочевидно, и пользователь
        # просто жал по ней впустую.
        self.mis_hint_lbl = ctk.CTkLabel(
            scroll, text="", text_color="gray60", justify="left", wraplength=760,
        )
        self.mis_hint_lbl.pack(anchor="w", pady=(0, 15))

        # Что скачивается прямо сейчас. Скачивание результатов MIS — это
        # десяток зашифрованных архивов на несколько гигабайт, и раньше на
        # всё это время шкала стояла на нуле с подписью "Скачивание
        # результатов MIS...", что неотличимо от зависания. Панель
        # появляется на время скачивания и убирается после.
        self.mis_files_box = ctk.CTkFrame(scroll, border_width=1, border_color="gray40")
        ctk.CTkLabel(
            self.mis_files_box, text="Скачивание результатов",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(anchor="w", padx=12, pady=(10, 4))
        self.mis_files_lbl = ctk.CTkLabel(
            self.mis_files_box, text="", justify="left",
            text_color="gray60", wraplength=740,
        )
        self.mis_files_lbl.pack(anchor="w", padx=12, pady=(0, 10))

        result_box = ctk.CTkFrame(scroll, border_width=1, border_color="gray40")
        result_box.pack(fill="x")
        self.result_box = result_box
        ctk.CTkLabel(
            result_box, text="Итоговый файл",
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(anchor="w", padx=12, pady=(10, 4))
        self.result_lbl = ctk.CTkLabel(
            result_box, text="", justify="left", text_color="gray60", wraplength=740,
        )
        self.result_lbl.pack(anchor="w", padx=12, pady=(0, 8))
        ctk.CTkButton(
            result_box, text="📁 Открыть папку с итоговым файлом", width=320,
            command=self._on_open_results_folder,
        ).pack(anchor="w", padx=12, pady=(0, 12))

        attach_input_features(self.curl_tf)
        attach_input_features(self.pwd_tf)
        self._refresh_result_label()

    # -----------------------------------------------------------------------
    # Вкладка "Лог"
    # -----------------------------------------------------------------------
    def _build_log_tab(self):
        toolbar = ctk.CTkFrame(self.tab_log, fg_color="transparent")
        toolbar.pack(fill="x", padx=10, pady=(10, 5))

        ctk.CTkButton(
            toolbar, text="Копировать весь лог", width=180,
            command=self._copy_all_log,
        ).pack(side="left")

        ctk.CTkButton(
            toolbar, text="Очистить лог", width=150,
            command=self._clear_log,
        ).pack(side="right")

        ctk.CTkLabel(
            toolbar, text="Текст в логе можно выделять мышью и копировать через Ctrl+C",
            text_color="gray60",
        ).pack(side="left", padx=20)

        self.log_text = ctk.CTkTextbox(self.tab_log, wrap="word", font=ctk.CTkFont(size=12))
        self.log_text.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        attach_input_features(self.log_text)

    # -----------------------------------------------------------------------
    # Вспомогательные методы
    # -----------------------------------------------------------------------
    def _pick_file(self, entry: ctk.CTkEntry):
        path = filedialog.askopenfilename()
        if path:
            entry.delete(0, "end")
            entry.insert(0, path)

    # -----------------------------------------------------------------------
    # Промт "обычные / продвинутые настройки"
    # -----------------------------------------------------------------------
    def _is_simple_mode(self) -> bool:
        return self.mode_switch.get() == MODE_SIMPLE

    def _on_mode_changed(self, choice: str | None = None):
        """
        Переключение между обычным и продвинутым режимом вкладки
        "Подготовка". Продвинутый блок (self.adv_box) не пересоздаётся, а
        снимается/возвращается через pack_forget()/pack() — виджеты в нём
        живы всегда, поэтому весь остальной код (валидация, чтение
        значений при запуске) работает одинаково в обоих режимах.

        pack() возвращает блок в конец self.adv_box'ового родителя
        (scroll), а базовый блок и сводка лежат внутри self.basic_box,
        который остаётся на месте — порядок элементов не нарушается.
        """
        simple = (choice or self.mode_switch.get()) == MODE_SIMPLE
        if simple:
            self.adv_box.pack_forget()
            self.simple_info_lbl.pack(anchor="w", pady=(0, 4))
            self.simple_tmpl_lbl.pack(anchor="w", pady=(0, 15))
            self._apply_simple_presets()
        else:
            self.simple_info_lbl.pack_forget()
            self.simple_tmpl_lbl.pack_forget()
            self.adv_box.pack(fill="x")
        _save_ui_mode(MODE_SIMPLE if simple else MODE_ADVANCED)

    def _on_source_changed(self):
        """
        В обычном режиме источник данных — единственная настройка, от
        которой зависят остальные (формат вывода и трафарет), поэтому при
        его смене пресеты пересчитываются. В продвинутом режиме смена
        источника ничего не трогает: там всё выбирает пользователь.
        """
        if self._is_simple_mode():
            self._apply_simple_presets()

    def _apply_simple_presets(self):
        """
        Проставляет в виджеты продвинутого блока значения, выбранные за
        пользователя в обычном режиме:

          * формат вывода   — по источнику (ftdna -> v3 LF, myheritage -> v5 CRLF);
          * трафарет        — samples/template_v3.txt / template_v5.txt;
          * порог Rsq       — 0.30 (стандартный порог MIS);
          * нормализация multiallelic-сайтов перед split — включена;
          * хранение сырых хромосом 1000 Genomes         — включено;
          * "использовать всех доступных EUR-доноров"     — выключено,
            вместо этого конкретное число доноров = 20.

        Значения пишутся именно в виджеты (а не подставляются в момент
        запуска), поэтому, переключившись в продвинутый режим, пользователь
        видит ровно то, что будет использовано.
        """
        source = self._get_source_key()
        fmt = SIMPLE_FORMAT_BY_SOURCE.get(source, "v3")

        # Референсная панель в обычном режиме всегда HRC (DEFAULT_PANEL):
        # это GRCh37, то есть та же сборка, в которой уже лежат координаты
        # чипа — не нужен ни лифтовер, ни второй комплект референса/доноров.
        panel_display = pipeline.REFERENCE_PANELS[pipeline.DEFAULT_PANEL]["display_name"]
        if self.panel_dd.get() != panel_display:
            self.panel_dd.set(panel_display)
            self._on_panel_changed()

        for value in self.format_dd.cget("values"):
            if value.startswith(fmt):
                self.format_dd.set(value)
                break

        self.rsq_entry.delete(0, "end")
        self.rsq_entry.insert(0, SIMPLE_RSQ)
        self._validate_rsq_entry()

        self.normalize_var.set(SIMPLE_NORMALIZE)
        self.raw_cache_var.set(SIMPLE_RAW_CACHE)

        self.eur_all_var.set(False)
        self.eur_count_entry.configure(state="normal")
        self.eur_count_entry.delete(0, "end")
        self.eur_count_entry.insert(0, str(SIMPLE_EUR_COUNT))
        self._on_eur_all_toggled()

        template = _find_sample_template(fmt)
        if template is not None:
            self.tmpl_tf.delete(0, "end")
            self.tmpl_tf.insert(0, str(template))
            self.simple_tmpl_lbl.configure(
                text=f"✓ Трафарет подставлен автоматически: {template}",
                text_color="#4CAF50",
            )
        else:
            expected = _samples_dir() / SAMPLE_TEMPLATE_NAMES.get(fmt, "template_v3.txt")
            self.simple_tmpl_lbl.configure(
                text=(f"⚠ Трафарет не найден: положите файл в {expected} "
                      f"— или переключитесь в продвинутый режим и укажите "
                      f"путь вручную."),
                text_color="#F9A825",
            )

        source_name = pipeline.SOURCES.get(source, {}).get("name", source)
        self.simple_info_lbl.configure(
            text=(
                f"Настройки подобраны автоматически под источник «{source_name}»:\n"
                f"    • референсная панель: {panel_display}\n"
                f"    • формат вывода: {fmt} "
                f"({'CRLF' if fmt == 'v5' else 'LF'})\n"
                f"    • порог Rsq: {SIMPLE_RSQ}\n"
                f"    • нормализация multiallelic-сайтов перед split: включена\n"
                f"    • хранение сырых хромосом 1000 Genomes: включено\n"
                f"    • число EUR-доноров: {SIMPLE_EUR_COUNT} "
                f"(не «все доступные»)\n"
                f"Чтобы изменить их вручную, выберите «{MODE_ADVANCED}» выше."
            )
        )

    def _validate_rsq_entry(self) -> bool:
        text = self.rsq_entry.get().strip()
        try:
            value = float(text)
            if not (0.30 <= value <= 0.99):
                raise ValueError
        except ValueError:
            self.rsq_entry.configure(border_color="#F44336")
            self.rsq_status_lbl.configure(
                text="⚠ Введите число от 0.30 до 0.99", text_color="#F44336",
            )
            return False
        self.rsq_entry.configure(border_color=("gray70", "gray30"))
        self.rsq_status_lbl.configure(
            text=f"✓ Порог принят: {value:.2f}", text_color="#4CAF50",
        )
        return True

    def _get_rsq_threshold(self) -> float:
        return float(self.rsq_entry.get().strip())

    def _on_eur_all_toggled(self):
        """
        Промт "настраиваемое количество EUR-доноров": включает/выключает
        поле ручного ввода числа доноров в зависимости от галочки
        "Использовать всех доступных EUR-доноров". Оба варианта взаимно
        исключающие — поле ввода имеет смысл только когда галочка снята.
        """
        if self.eur_all_var.get():
            self.eur_count_entry.configure(state="disabled")
            self.eur_count_status_lbl.configure(
                text="✓ Будут использованы все доступные EUR-доноры (~500)",
                text_color="#4CAF50",
            )
        else:
            self.eur_count_entry.configure(state="normal")
            self._validate_eur_count_entry()

    def _validate_eur_count_entry(self) -> bool:
        """Валидация ручного ввода числа EUR-доноров — тот же паттерн, что
        и _validate_rsq_entry(). Актуально, только пока галочка "все
        доступные" снята; при включённой галочке поле неактивно и его
        содержимое не используется."""
        if self.eur_all_var.get():
            return True
        text = self.eur_count_entry.get().strip()
        try:
            value = int(text)
            if not (1 <= value <= download_donors.MAX_EUR_SAMPLE_COUNT):
                raise ValueError
        except ValueError:
            self.eur_count_entry.configure(border_color="#F44336")
            self.eur_count_status_lbl.configure(
                text=f"⚠ Введите целое число от 1 до {download_donors.MAX_EUR_SAMPLE_COUNT}",
                text_color="#F44336",
            )
            return False
        self.eur_count_entry.configure(border_color=("gray70", "gray30"))
        self.eur_count_status_lbl.configure(
            text=f"✓ Будет использовано {value} EUR-доноров", text_color="#4CAF50",
        )
        return True

    def _get_eur_sample_count(self) -> int | None:
        """
        None — использовать всю доступную EUR-подвыборку панели (галочка
        включена, поведение по умолчанию, см.
        download_donors.EUR_SAMPLE_COUNT_ALL). Иначе — явное число из
        self.eur_count_entry.
        """
        if self.eur_all_var.get():
            return None
        return int(self.eur_count_entry.get().strip())

    def _toggle_pwd(self):
        self._pwd_visible = not self._pwd_visible
        if self._pwd_visible:
            self.pwd_tf.configure(show="")
            self.show_pwd_btn.configure(text="🔒")
        else:
            self.pwd_tf.configure(show="*")
            self.show_pwd_btn.configure(text="👁")

    def _on_diagnose_password(self):
        """
        Диагностика пароля (Задача 3/7) — теперь целиком через
        core.archive_utils.diagnose_password(), без дублирования логики
        поиска 7z.exe и санитайзинга пароля.

        ⚠ Фикс: раньше и test_archive, и поиск 7z.exe были жёстко
        завязаны на PROJECT_ROOT/"bin", независимо от того, что
        пользователь ввёл в self.bin_tf ("Папка с бинарниками htslib") —
        та же папка, куда обычно кладут 7z.exe вместе с bcftools.exe/
        tabix.exe. Из-за этого diagnose_password() мог ложно сообщать
        "7z.exe НЕ найден", если пользователь указал другую папку
        бинарников, и тестовый архив _password_test.7z (если он лежит
        именно в выбранной пользователем папке) не находился вовсе.
        Теперь оба пути строятся от self.bin_tf, с откатом на
        PROJECT_ROOT/"bin" только если поле не заполнено.
        """
        pwd = self.pwd_tf.get()
        bd = Path(self.bin_tf.get()) if self.bin_tf.get() else _detect_bin_dir()
        test_archive = bd / "_password_test.7z"
        sevenzip_candidate = str(bd / "7z.exe") if os.name == "nt" else str(bd / "7z")
        lines = archive_utils.diagnose_password(pwd, test_archive, sevenzip_path=sevenzip_candidate)
        messagebox.showinfo("Диагностика пароля", "\n".join(lines))

    def _on_diagnose_network(self):
        """
        v13 (промт "Диагностика реальной удалённой фильтрации + устойчивая
        настройка CA-сертификатов"): кнопка "🌐 Диагностика сети" на
        вкладке "Подготовка". Проверяет ОБА фикса из network_utils.py
        (CA-сертификаты + конфликт bin/curl.exe) и затем реальную
        удалённую фильтрацию (download_donors.diagnose_remote_filter(),
        не только "-h" — см. докстринг функции), а не просто пробует
        открыть URL.

        Запускается в фоновом потоке — реальный тест фильтрации может
        занять до REMOTE_CHROM_TIMEOUT секунд при плохой сети/битом
        зеркале, блокировать интерфейс на это время не годится. Результат
        показывается через self.after(0, ...), как и остальные
        thread-safe диалоги в этом классе.
        """
        bd = Path(self.bin_tf.get()) if self.bin_tf.get() else None

        def _worker():
            htslib = download_donors.HtslibTools(bd)
            lines: list[str] = []
            if not htslib.has_bcftools:
                lines.append("bcftools не найден — укажите папку бинарников")
                self.after(0, messagebox.showinfo, "Диагностика сети", "\n".join(lines))
                return

            network_utils.ensure_network_ready(bd)
            ca_ok = bool(os.environ.get("CURL_CA_BUNDLE"))
            lines.append(
                f"{'✓' if ca_ok else '⚠'} CA-сертификаты (CURL_CA_BUNDLE): "
                f"{os.environ.get('CURL_CA_BUNDLE', 'не установлены')}"
            )
            conflicting = network_utils.find_conflicting_bin_curl(bd)
            if conflicting:
                lines.append(
                    f"⚠ В папке бинарников найден собственный {conflicting.name} — "
                    f"может конфликтовать с системным (см. лог для подробностей). "
                    f"Приложение игнорирует его при поиске curl для скачивания."
                )
            else:
                lines.append("✓ Конфликтующий curl.exe в папке бинарников не найден")

            lines.append("")
            lines.append("Проверка реальной удалённой фильтрации bcftools...")
            report = download_donors.diagnose_remote_filter(
                htslib, PROJECT_ROOT / "output",
            )
            if report["ok"]:
                lines.append(
                    f"✓ Удалённая фильтрация РАБОТАЕТ: {report['records']} "
                    f"записей за {report['duration_sec']:.1f}с"
                )
                lines.append(f"  URL: {report['url']}")
            else:
                lines.append("Удалённая фильтрация НЕ работает — будет использовано полное скачивание (это нормально, просто медленнее)")
                if report["stderr_tail"]:
                    lines.append(f"  Причина: {report['stderr_tail']}")

            self.after(0, messagebox.showinfo, "Диагностика сети", "\n".join(lines))

        threading.Thread(target=_worker, daemon=True).start()

    def _copy_all_log(self):
        content = self.log_text.get("1.0", "end").strip()
        if not content:
            messagebox.showinfo("Лог пуст", "В логе пока нет сообщений.")
            return
        self.clipboard_clear()
        self.clipboard_append(content)
        messagebox.showinfo("Скопировано", "Весь лог скопирован в буфер обмена.")

    def _clear_log(self):
        self.log_text.delete("1.0", "end")
        # Метки прогресса ссылались на позиции внутри только что стёртого
        # текста — без сброса следующий _upsert_progress_line() решил бы,
        # что метка ещё жива, и мог бы обновить не ту строку.
        self._progress_keys.clear()

    def _get_source_key(self) -> str:
        source_name = self.source_dd.get()
        return next((k for k, v in pipeline.SOURCES.items()
                     if v["name"] == source_name), "ftdna")

    def _get_panel_key(self) -> str:
        """
        Шаг 1 промта "HRC / TopMed": читает выбранную панель из
        self.panel_dd и возвращает её ключ ('hrc'/'topmed') для
        REFERENCE_PANELS. Откатывается на pipeline.DEFAULT_PANEL, если
        по какой-то причине название не нашлось (не должно происходить
        в норме — защита от рассинхронизации списка панелей).
        """
        panel_display = self.panel_dd.get()
        return next(
            (k for k, v in pipeline.REFERENCE_PANELS.items()
             if v["display_name"] == panel_display),
            pipeline.DEFAULT_PANEL,
        )

    def _get_format_key(self) -> str:
        fmt = self.format_dd.get()
        return "v5" if "v5" in fmt else "v3"

    def _on_panel_changed(self):
        """
        Callback при смене self.panel_dd. Обновляет предупреждение под
        выпадающим списком и динамический текст инструкции на вкладке
        "Запуск".

        Промт "встроить лифтовер HRC/TopMed в gui/app.py": раньше здесь
        было предупреждение "лифтовер координат ещё не выполняется" — это
        устарело, лифтовер (core/liftover.py::ChainLiftover,
        pipeline._build_liftover()) реализован и подключается в
        _run_stages_1_6()/_run_stage_7(). Текст ниже — нейтральное
        информационное сообщение о том, что для не-HRC панели скачиваются
        дополнительные файлы (референс другой сборки, доноры, chain-файлы
        лифтовера), а не предупреждение о некорректности результата.

        Единственное реальное ограничение, которое всё ещё может повлиять
        на корректность результата — source='vcf' не поддерживает
        лифтовер (pipeline._supports_liftover()). Это предупреждение
        сознательно НЕ дублируется здесь: оно показывается один раз,
        непосредственно перед запуском, в _run_stages_1_6() — там уже
        известен финальный source (после возможного переключения по
        автодетекту), тогда как здесь в момент смены панели source мог бы
        быть неактуален.
        """
        panel = self._get_panel_key()
        cfg = pipeline.REFERENCE_PANELS[panel]
        if panel == pipeline.DEFAULT_PANEL:
            self.panel_warning_lbl.configure(text="")
        else:
            self.panel_warning_lbl.configure(
                text=(
                    f"ℹ Панель «{cfg['display_name']}» использует сборку генома "
                    f"{cfg['genome_build'].upper()}, отличную от HRC (GRCh37). "
                    f"Референс, доноры и chain-файл лифтовера координат для "
                    f"этой панели хранятся отдельно от HRC и будут скачаны при "
                    f"первом запуске — это может занять дополнительное время "
                    f"(десятки МБ для chain-файла лифтовера, гигабайты для "
                    f"референса/доноров)."
                )
            )
        # Пустое предупреждение всё равно занимало вертикальный отступ и
        # оставляло дыру в макете — убираем метку с экрана, когда текста нет.
        if self.panel_warning_lbl.cget("text"):
            if not self.panel_warning_lbl.winfo_manager():
                self.panel_warning_lbl.pack(
                    anchor="w", pady=(0, 15), after=self.panel_dd,
                )
        else:
            self.panel_warning_lbl.pack_forget()

        # Метод может вызываться до построения вкладки "Запуск" (первичная
        # синхронизация в конце _build_settings_tab) — тогда просто пропускаем.
        if hasattr(self, "run_instructions_lbl"):
            self._refresh_run_instructions()

    def _refresh_run_instructions(self):
        panel = self._get_panel_key()
        cfg = pipeline.REFERENCE_PANELS[panel]
        build = "GRCh38/hg38" if cfg["genome_build"] == "grch38" else "GRCh37/hg19"
        # Текст держим коротким: кнопки под ним открывают сайт и папку, а
        # значения формы вынесены в отдельный заметный блок ниже — в
        # абзаце они терялись.
        self.run_instructions_lbl.configure(text=(
            "Этот шаг делается руками на сайте импутации:\n"
            "1. Откройте сайт и загрузите на него 22 файла из папки запуска "
            "(обе кнопки ниже).\n"
            "2. В форме выберите параметры из синей рамки.\n"
            "3. Дождитесь письма со ссылкой и паролем и вставьте их в Шаге 3.\n"
            "💡 Ctrl+V для вставки, правая кнопка мыши — контекстное меню"
        ))

        # Блок параметров создаётся позже самого первого вызова этого
        # метода (панель синхронизируется ещё на вкладке "Подготовка") —
        # поэтому проверяем наличие, как и с run_instructions_lbl выше.
        if hasattr(self, "mis_params_lbl"):
            self.mis_params_lbl.configure(text=(
                f"Reference Panel:  {cfg['mis_panel_value']}\n"
                f"Array Build:  {build}\n"
                f"Population:  EUR"
            ))

    # -----------------------------------------------------------------------
    # Промт "Именованные папки запуска": имя/история запусков
    # -----------------------------------------------------------------------
    def _refresh_run_name_suggestion(self):
        """Предзаполняет self.run_name_tf следующим свободным номером
        запуска (output/runs/<N>/) — редактируемо пользователем до
        нажатия «Запустить этапы 1-6»."""
        try:
            runs_root = PROJECT_ROOT / "output" / pipeline.RUNS_SUBDIR_NAME
            suggestion = pipeline._next_run_name(runs_root)
        except Exception:
            suggestion = "1"
        self.run_name_tf.delete(0, "end")
        self.run_name_tf.insert(0, suggestion)

    def _refresh_run_history(self):
        """
        Обновляет self.run_history_dd списком существующих запусков
        (output/runs/*, новые сверху) с человекочитаемыми подписями из
        run_info.json (pipeline.format_run_label()). self._run_history_map
        хранит соответствие "подпись -> Path папки запуска", потому что
        CTkOptionMenu показывает пользователю только текст подписи.
        """
        runs = pipeline.list_runs(PROJECT_ROOT / "output")
        self._run_history_map = {}
        if not runs:
            self.run_history_dd.configure(values=["(нет запусков)"])
            self.run_history_dd.set("(нет запусков)")
            return

        labels: list[str] = []
        for run_dir in runs:
            base_label = pipeline.format_run_label(run_dir)
            label = base_label
            n = 2
            while label in self._run_history_map:
                label = f"{base_label} ({n})"
                n += 1
            self._run_history_map[label] = run_dir
            labels.append(label)

        self.run_history_dd.configure(values=labels)
        # Не перетираем выбор пользователя, если он всё ещё существует.
        if self.run_history_dd.get() not in labels:
            self.run_history_dd.set(labels[0])

    def _selected_history_run(self) -> Path | None:
        label = self.run_history_dd.get()
        return self._run_history_map.get(label)

    def _detach_run_log_handler(self):
        """Снимает и закрывает FileHandler текущего активного запуска (если
        есть) — обязательно перед переименованием его папки на Windows
        (открытый файловый дескриптор мешает rename) и перед переключением
        на другой активный запуск."""
        if self._run_log_handler is not None:
            logging.getLogger().removeHandler(self._run_log_handler)
            try:
                self._run_log_handler.close()
            except Exception:
                pass
            self._run_log_handler = None

    def _set_active_run(self, run_dir: Path, run_name: str):
        """
        Делает run_dir активным запуском (self.current_run_dir/
        self.current_run_name) для всего последующего кода
        _run_stages_1_6()/_run_stage_7() — единая точка входа что для
        только что созданного нового запуска (_on_start), что для
        выбранного из истории (_on_continue_run/_on_rename_run).
        """
        self._detach_run_log_handler()
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)
        self.current_run_dir = run_dir
        self.current_run_name = run_name
        self._run_log_handler = pipeline.attach_run_log_handler(run_dir)
        self.active_run_lbl.configure(text=f"Активный запуск: {run_name} ({run_dir})")
        self.run_name_tf.delete(0, "end")
        self.run_name_tf.insert(0, run_name)

    def _on_continue_run(self):
        """
        Кнопка «▶ Продолжить (Шаг 3)» — делает выбранный из истории
        запуск активным без повторного прогона Этапов 1-6. Основной
        сценарий: пользователь получил письмо MIS уже после перезапуска
        GUI, и self.current_run_dir из предыдущей сессии потерян.
        """
        if self.running:
            return
        run_dir = self._selected_history_run()
        if run_dir is None:
            messagebox.showwarning("Предупреждение", "Выберите запуск из списка истории")
            return
        if not (run_dir / "parse_result.pkl").exists():
            if not messagebox.askyesno(
                "Внимание",
                f"В папке запуска «{run_dir.name}» не найден parse_result.pkl "
                f"— похоже, Этап 1-6 для него не был завершён. Этап 7 "
                f"(скачивание результатов MIS) не сможет собрать финальный "
                f"файл без него.\n\nВсё равно сделать этот запуск активным?",
            ):
                return
        self._set_active_run(run_dir, run_dir.name)
        self._set_wizard_step(2)
        self._refresh_mis_btn_state()
        messagebox.showinfo(
            "Запуск выбран",
            f"Активный запуск: «{run_dir.name}».\n"
            f"Вставьте curl-команду и пароль из письма MIS и нажмите "
            f"«📥 Скачать результаты и собрать финальный файл».",
        )

    def _on_rename_run(self):
        """Кнопка «✏ Переименовать» — переименовывает папку запуска на
        диске (output/runs/<старое> -> output/runs/<новое>), с валидацией
        имени (pipeline.validate_run_name) и защитой от коллизии."""
        if self.running:
            messagebox.showwarning(
                "Предупреждение", "Дождитесь завершения текущего запуска перед переименованием",
            )
            return
        run_dir = self._selected_history_run()
        if run_dir is None:
            messagebox.showwarning("Предупреждение", "Выберите запуск из списка истории")
            return

        dialog = ctk.CTkInputDialog(
            text=f"Новое имя для запуска «{run_dir.name}»:",
            title="Переименовать запуск",
        )
        new_name = dialog.get_input()
        if not new_name:
            return
        try:
            new_name = pipeline.validate_run_name(new_name)
        except ValueError as e:
            messagebox.showerror("Ошибка", str(e))
            return

        new_dir = run_dir.parent / new_name
        if new_dir.exists():
            messagebox.showerror("Ошибка", f"Папка запуска «{new_name}» уже существует")
            return

        was_active = (self.current_run_dir is not None and self.current_run_dir == run_dir)
        if was_active:
            # Открытый FileHandler на run.log мешает переименованию папки
            # на Windows — снимаем его перед rename() и переустанавливаем
            # на новый путь сразу после.
            self._detach_run_log_handler()
        try:
            run_dir.rename(new_dir)
        except OSError as e:
            messagebox.showerror("Ошибка", f"Не удалось переименовать папку: {e}")
            if was_active:
                # Переименование не удалось — восстанавливаем handler на
                # старом пути, чтобы логирование текущего запуска не оборвалось.
                self._run_log_handler = pipeline.attach_run_log_handler(run_dir)
            return

        if was_active:
            self._set_active_run(new_dir, new_name)
        self._refresh_run_history()
        messagebox.showinfo("Готово", f"Запуск переименован в «{new_name}»")

    def _on_open_run_folder(self):
        """Кнопка «📂 Открыть папку» — открывает папку выбранного в списке
        запуска, а если в списке ничего нет — папку активного запуска."""
        run_dir = self._selected_history_run() or self.current_run_dir
        if run_dir is None:
            messagebox.showwarning(
                "Предупреждение",
                "Пока нет ни одного запуска. Запустите Шаг 1 — папка "
                "создастся автоматически.",
            )
            return
        self._open_in_file_manager(Path(run_dir))

    # -----------------------------------------------------------------------
    # Промт "сделать вкладку Запуск юзерфрендли"
    # -----------------------------------------------------------------------
    def _set_wizard_step(self, step: int):
        """
        Переключает подвкладку "Запуска" на нужный шаг и помечает
        пройденные галочкой в самой подписи вкладки. step: 1 (подготовка),
        2 (импутация на MIS), 3 (сборка файла).

        CTkTabview.rename() меняет ключ, по которому работает .set(),
        поэтому актуальные подписи держатся в self._run_tab_names — без
        этого второй вызов .set() ушёл бы к несуществующему имени.
        """
        self._wizard_step = step
        for i, base in enumerate(RUN_TAB_BASE_NAMES, start=1):
            desired = f"✓ {base}" if i < step else base
            current = self._run_tab_names[i - 1]
            if current != desired:
                self.run_tabs.rename(current, desired)
                self._run_tab_names[i - 1] = desired
        self.run_tabs.set(self._run_tab_names[step - 1])

    def _toggle_run_details(self):
        """Разворачивает/сворачивает блок подробностей о папке запуска."""
        if self._run_details_visible:
            self.run_details_box.pack_forget()
            self.run_details_btn.configure(text="▾ Подробнее")
        else:
            # Возвращаем блок на его место — сразу после active_run_lbl, а
            # не в конец контейнера (pack по умолчанию добавляет в хвост).
            self.run_details_box.pack(fill="x", after=self.active_run_lbl)
            self.run_details_btn.configure(text="▴ Свернуть")
        self._run_details_visible = not self._run_details_visible

    def _show_run_details(self):
        if not self._run_details_visible:
            self._toggle_run_details()

    def _on_new_run(self):
        """
        «＋ Новый»: предлагает следующее свободное имя запуска и
        разворачивает подробности, чтобы имя можно было сразу поправить.
        Саму папку не создаёт — это делает _on_start() при запуске Шага 1.
        """
        if self.running:
            messagebox.showwarning(
                "Предупреждение", "Дождитесь завершения текущего запуска",
            )
            return
        self._refresh_run_name_suggestion()
        self._show_run_details()
        self._set_wizard_step(1)

    def _current_upload_dir(self) -> Path | None:
        """Папка с 22 файлами для загрузки на MIS у активного запуска."""
        if self.current_run_dir is None:
            return None
        return self.current_run_dir / "upload"

    def _on_open_upload_folder(self):
        """
        «📂 Папка с 22 файлами»: открывает output/runs/<запуск>/upload —
        иначе пользователю приходится искать её в проводнике руками.
        """
        upload_dir = self._current_upload_dir()
        if upload_dir is None:
            messagebox.showinfo(
                "Папка ещё не создана",
                "Сначала выполните Шаг 1 (подготовку файлов) или выберите "
                "готовый запуск в списке и нажмите «▶ Продолжить (Шаг 3)».",
            )
            return
        if not upload_dir.is_dir():
            messagebox.showinfo(
                "Папка ещё не создана",
                f"У запуска «{self.current_run_name}» ещё нет папки с файлами "
                f"для загрузки:\n{upload_dir}\n\nОна появляется в конце Шага 1.",
            )
            return
        self._open_in_file_manager(upload_dir)

    def _on_open_results_folder(self):
        """
        «📁 Открыть папку с итоговым файлом» — открывает results/ рядом с
        программой. Туда складываются ИТОГОВЫЕ файлы всех запусков, отдельно
        от рабочих папок output/runs/<...>, где лежат десятки промежуточных
        VCF и логов и где готовый файл терялся из виду.
        """
        self._open_in_file_manager(_results_dir())

    def _refresh_result_label(self):
        """Строка о том, где лежит итоговый файл (или где он появится)."""
        if not hasattr(self, "result_lbl"):
            return
        if self._last_result_path is not None:
            self.result_lbl.configure(
                text=f"✓ Готов: {self._last_result_path}", text_color="#4CAF50",
            )
        else:
            self.result_lbl.configure(
                text=(f"Появится после сборки в папке {_results_dir()} — "
                      f"туда складываются итоговые файлы всех запусков."),
                text_color="gray60",
            )

    def _open_in_file_manager(self, path: Path):
        """Открывает папку в проводнике; на не-Windows — xdg-open, а если и
        его нет, просто показывает путь, чтобы его можно было скопировать."""
        path = Path(path)
        if os.name == "nt":
            os.startfile(str(path))
            return
        try:
            subprocess.run(["xdg-open", str(path)], check=False)
        except Exception:
            messagebox.showinfo("Папка", str(path))

    def _on_open_mis_site(self):
        """
        «🌐 Открыть сервер импутации» — адрес берётся из конфигурации той
        панели, которая выбрана на вкладке "Подготовка": TOPMed r3 живёт
        не на Michigan, а на BioData Catalyst (см. REFERENCE_PANELS в
        main.py), и открывать всегда Michigan было бы ошибкой.
        """
        cfg = pipeline.REFERENCE_PANELS[self._get_panel_key()]
        webbrowser.open(cfg["mis_upload_url"])
        self.mis_actions_status_lbl.configure(
            text=f"Открыт {cfg['mis_upload_url']}", text_color="gray60",
        )

    # --- Наблюдение за файлами доноров на диске -------------------------
    def _start_file_watch(self, dirs, target: str = "donors"):
        """
        Запускает опрос папок с качающимися файлами (раз в
        _DONOR_WATCH_INTERVAL_MS). Вызывается дважды за прогон: перед
        этапом 3 Шага 1 (донорские хромосомы, target="donors") и перед
        скачиванием результатов MIS на Шаге 3 (target="mis").

        Почему опрос диска, а не разбор вывода загрузчика: сколько именно
        мегабайт скачано, надёжно знает только сам файл на диске. Вывод
        зависит от того, каким инструментом идёт закачка (aria2c печатает
        «1.0GiB/1.2GiB(83%)», curl --progress-bar — только полоску из
        решёток с процентом, bcftools при удалённой фильтрации не печатает
        прогресса вообще), от того, установлен ли aria2c у пользователя, и
        от прореживания вывода по времени. Размер файла не зависит ни от
        чего из этого и растёт всегда, когда закачка жива, — именно на
        этот вопрос («оно шевелится или зависло?») пользователь и смотрит.
        """
        self._donor_watch_dirs = [Path(d) for d in dirs if d]
        self._donor_file_sizes = {}
        self._watch_target = target
        if not self._donor_watch_dirs:
            return
        if target == "donors":
            self._show_donor_panel()
        else:
            self.mis_files_box.pack(fill="x", pady=(0, 15), before=self.result_box)
            self.mis_files_lbl.configure(
                text="Жду начала скачивания...", text_color="gray60",
            )
            # На этом отрезке общий объём заранее неизвестен (сколько
            # архивов пришлёт MIS и какого размера — видно только по факту),
            # поэтому вместо застывшей на нуле шкалы честнее бегущая
            # полоса: «работаю, но сколько осталось — не знаю». Точные
            # цифры при этом идут строкой ниже, по файлам.
            self.progress.configure(mode="indeterminate")
            self.progress.start()
        if self._donor_watch_id is None:
            self._poll_watch_files()

    def _stop_file_watch(self):
        if self._watch_target == "mis":
            try:
                self.progress.stop()
                self.progress.configure(mode="determinate")
                self.progress.set(0)
            except Exception:
                pass
        self._donor_watch_dirs = []
        if self._donor_watch_id is not None:
            try:
                self.after_cancel(self._donor_watch_id)
            except Exception:
                pass
            self._donor_watch_id = None
        if hasattr(self, "donor_active_lbl"):
            self.donor_active_lbl.configure(text="")

    def _poll_watch_files(self):
        """Один цикл опроса: обновляет ячейки растущих файлов и список
        активных закачек. Перепланирует сам себя, пока идёт скачивание."""
        self._donor_watch_id = None
        if not self._donor_watch_dirs:
            return

        now = time.monotonic()
        donors = self._watch_target == "donors"
        active: list[tuple[str, int, float]] = []
        total_files = 0
        total_bytes = 0
        for directory in self._donor_watch_dirs:
            try:
                entries = list(directory.iterdir())
            except OSError:
                continue
            for path in entries:
                # На Шаге 1 растут только VCF-файлы доноров; на Шаге 3
                # приходят zip-архивы MIS и распакованное из них, поэтому
                # там смотрим на всё подряд.
                if donors and ".vcf" not in path.name:
                    continue
                try:
                    if not path.is_file():
                        continue
                    size = path.stat().st_size
                except OSError:
                    continue
                total_files += 1
                total_bytes += size
                previous = self._donor_file_sizes.get(path)
                self._donor_file_sizes[path] = (now, size)
                if previous is None:
                    continue
                prev_time, prev_size = previous
                if size <= prev_size:
                    continue
                speed = (size - prev_size) / max(0.001, now - prev_time)
                active.append((path.name, size, speed))

                match = _DONOR_CHR_RE.search(path.name) if donors else None
                if match:
                    chrom = int(match.group(1))
                    if 1 <= chrom <= 22:
                        kind = "индекс" if path.name.endswith(".tbi") else "⬇"
                        self._set_donor_cell(
                            chrom,
                            f"{kind} {_fmt_size(size)} ({_fmt_size(speed)}/с)",
                            "#42A5F5", _DONOR_RANK_DOWNLOAD,
                        )

        target_lbl = self.donor_active_lbl if donors else self.mis_files_lbl
        if active:
            active.sort(key=lambda item: -item[1])
            lines = [
                f"⬇ {name} — скачано {_fmt_size(size)}, {_fmt_size(speed)}/с"
                for name, size, speed in active[:4]
            ]
            if len(active) > 4:
                lines.append(f"… и ещё {len(active) - 4} файл(ов)")
            if not donors:
                lines.append(
                    f"Всего в папке результатов: {total_files} файл(ов), "
                    f"{_fmt_size(total_bytes)}"
                )
            target_lbl.configure(text="\n".join(lines), text_color="#42A5F5")
        elif donors:
            # Пусто — не значит "зависло": между хромосомами идёт
            # фильтрация bcftools, она диск почти не растит.
            target_lbl.configure(
                text="Сейчас ничего не скачивается — идёт фильтрация/индексация.",
                text_color="gray60",
            )
        else:
            target_lbl.configure(
                text=(f"Сейчас ничего не скачивается — идёт распаковка или "
                      f"проверка архивов.\nВсего в папке результатов: "
                      f"{total_files} файл(ов), {_fmt_size(total_bytes)}"),
                text_color="gray60",
            )

        self._donor_watch_id = self.after(
            _DONOR_WATCH_INTERVAL_MS, self._poll_watch_files,
        )

    # --- Живая карта донорских хромосом (этап 3) -------------------------
    def _reset_donor_panel(self):
        """Сбрасывает карту хромосом и убирает её с экрана — вызывается в
        начале нового запуска, чтобы не смешивать состояния разных прогонов."""
        if not hasattr(self, "donor_chr_lbls"):
            return
        for chrom, lbl in self.donor_chr_lbls.items():
            lbl.configure(text=f"chr{chrom} —", text_color="gray50")
        self.donor_summary_lbl.configure(text="")
        self.donor_active_lbl.configure(text="")
        self.donor_panel.pack_forget()
        self._donor_states = {}
        self._donor_file_sizes = {}

    def _update_donor_panel(self, msg: str):
        """
        Обновляет карту хромосом по строке лога скачивания доноров.
        Разбор — в _parse_donor_state() (модульная функция, чтобы её было
        видно и тестировать отдельно от виджетов).
        """
        parsed = _parse_donor_state(msg)
        if parsed is None:
            return
        self._set_donor_cell(*parsed)

    def _show_donor_panel(self):
        if hasattr(self, "donor_panel") and not self.donor_panel.winfo_manager():
            self.donor_panel.pack(fill="x", pady=(10, 0))

    def _set_donor_cell(self, chrom: int, text: str, color: str, rank: int):
        """
        Ставит состояние одной хромосоме. Единая точка входа и для разбора
        строк лога (_update_donor_panel), и для опроса файлов на диске
        (_poll_donor_files) — счётчик готовых и показ панели считаются
        в одном месте, а не дублируются.
        """
        lbl = self.donor_chr_lbls.get(chrom)
        if lbl is None:
            return

        # Сообщения от параллельных потоков приходят вперемешку, и строка
        # прогресса скачивания может прийти уже ПОСЛЕ "готово" (её напечатал
        # другой поток чуть раньше). Ранг не даёт состоянию откатиться назад.
        # Равные ранги разрешены: так обновляются мегабайты внутри закачки.
        if self._donor_states.get(chrom, (0,))[0] > rank:
            return
        self._donor_states[chrom] = (rank, text)
        lbl.configure(text=f"chr{chrom} {text}", text_color=color)
        self._show_donor_panel()

        # Строго == DONE: у ранга FAILED число больше (ошибку не должно
        # перекрывать запоздалое "готово" от другого потока), но в счётчик
        # готовых упавшая хромосома, разумеется, не входит.
        done = sum(1 for r, _ in self._donor_states.values() if r == _DONOR_RANK_DONE)
        failed = sum(1 for r, _ in self._donor_states.values() if r == _DONOR_RANK_FAILED)
        summary = f"Готово {done} из 22"
        if failed:
            summary += f", с ошибкой {failed}"
        self.donor_summary_lbl.configure(text=summary)

    def _validate_curl_entry(self) -> bool:
        """
        Проверяет вставленную curl-команду сразу при вставке, а не при
        нажатии кнопки — тот же приём, что и у порога Rsq на вкладке
        "Подготовка". Возвращает True, если поле непустое (кнопку Шага 3
        не блокируем по одной лишь эвристике: письма MIS со временем
        меняют формат, и ложное срабатывание не должно останавливать
        работу — непохожий текст только помечается предупреждением).
        """
        text = self.curl_tf.get("1.0", "end").strip()
        if not text:
            self.curl_status_lbl.configure(text="", text_color="gray60")
            self._refresh_mis_btn_state()
            return False

        looks_like_curl = "curl" in text.lower() and "http" in text.lower()
        job = _MIS_JOB_RE.search(text)
        if looks_like_curl and job:
            self.curl_status_lbl.configure(
                text=f"✓ Распознано задание {job.group(0)}", text_color="#4CAF50",
            )
        elif looks_like_curl:
            self.curl_status_lbl.configure(
                text="✓ Похоже на curl-команду (идентификатор задания не "
                     "распознан — это нормально, формат письма мог измениться)",
                text_color="#4CAF50",
            )
        else:
            self.curl_status_lbl.configure(
                text="⚠ Не похоже на curl-команду из письма MIS. Нужна строка, "
                     "которая начинается с «curl» и содержит ссылку — скопируйте "
                     "её из письма целиком. Кнопку ниже это не блокирует.",
                text_color="#F9A825",
            )
        self._refresh_mis_btn_state()
        return True

    def _refresh_mis_btn_state(self):
        """
        Единая точка правды о доступности кнопки Шага 3 и о том, ПОЧЕМУ она
        серая. Раньше кнопка просто включалась после Шага 1, а недостающие
        curl/пароль обнаруживались только по нажатию — во всплывающем окне.
        """
        if self.running:
            self.mis_btn.configure(state="disabled")
            self.mis_hint_lbl.configure(
                text="Идёт выполнение — дождитесь завершения.", text_color="gray60",
            )
            return

        missing: list[str] = []
        if self.current_run_dir is None:
            missing.append("выполнить Шаг 1 (или выбрать готовый запуск и нажать "
                           "«▶ Продолжить (Шаг 3)» в подробностях)")
        if not self.curl_tf.get("1.0", "end").strip():
            missing.append("вставить curl-команду из письма MIS")
        if not self.pwd_tf.get().strip():
            missing.append("вставить пароль из письма")

        if missing:
            self.mis_btn.configure(state="disabled")
            self.mis_hint_lbl.configure(
                text="Кнопка станет активной, когда: " + "; ".join(missing) + ".",
                text_color="gray60",
            )
        else:
            self.mis_btn.configure(state="normal")
            self.mis_hint_lbl.configure(
                text="✓ Всё готово — можно скачивать результаты и собирать файл.",
                text_color="#4CAF50",
            )

    def _update_eta(self, frac: float):
        """
        Оценка оставшегося времени по доле выполненного. Оценка грубая
        (этапы очень разные по длительности: скачивание доноров занимает
        часы, разбивка по хромосомам — минуты), поэтому в тексте стоит
        "примерно" и точное время не обещается. До 3% прогресса вообще
        ничего не показываем — там оценка бессмысленна.
        """
        if self._run_started_at is None or frac < 0.03:
            self.eta_lbl.configure(text="")
            return
        elapsed = time.monotonic() - self._run_started_at
        remaining = elapsed * (1.0 - frac) / frac
        self.eta_lbl.configure(
            text=f"Прошло {_fmt_duration(elapsed)} · осталось примерно "
                 f"{_fmt_duration(remaining)}"
        )

    def _notify_done(self, success: bool):
        """
        Звук + мигание кнопки в панели задач по окончании длинной операции:
        Шаг 1 идёт часами, и к экрану в этот момент обычно никто не сидит.
        """
        try:
            self.bell()
        except Exception:
            pass
        self._flash_taskbar()
        if success:
            self.eta_lbl.configure(text="Готово.", text_color="#4CAF50")
        else:
            self.eta_lbl.configure(text="Завершено с ошибкой — см. лог.",
                                   text_color="#F44336")

    def _flash_taskbar(self):
        """
        Мигание кнопки приложения в панели задач Windows (FlashWindowEx).
        Только Windows; любые ошибки проглатываются — уведомление не та
        вещь, из-за которой приложение имеет право упасть.
        """
        if os.name != "nt":
            return
        try:
            import ctypes
            from ctypes import wintypes

            class _FLASHWINFO(ctypes.Structure):
                _fields_ = [
                    ("cbSize", wintypes.UINT),
                    ("hwnd", wintypes.HWND),
                    ("dwFlags", wintypes.DWORD),
                    ("uCount", wintypes.UINT),
                    ("dwTimeout", wintypes.DWORD),
                ]

            # winfo_id() возвращает дочернее окно Tk, а мигать должно
            # окно верхнего уровня — берём его через GetParent().
            hwnd = ctypes.windll.user32.GetParent(self.winfo_id()) or self.winfo_id()
            info = _FLASHWINFO(
                ctypes.sizeof(_FLASHWINFO), hwnd,
                0x00000003 | 0x0000000C,  # FLASHW_ALL | FLASHW_TIMERNOFG
                5, 0,
            )
            ctypes.windll.user32.FlashWindowEx(ctypes.byref(info))
        except Exception:
            pass

    # -----------------------------------------------------------------------
    # Промт "обратная связь автору"
    # -----------------------------------------------------------------------
    def _collect_diagnostics(self) -> str:
        """
        Технические данные, которые подставляются в письмо: версия, система
        и настройки запуска. Сюда сознательно НЕ попадают пути к файлам —
        в них видно имя пользователя Windows и имя генетического файла.
        Лог, где такое встречается, прикладывается только по явной галочке.
        """
        try:
            import platform
            os_line = f"{platform.system()} {platform.release()} ({platform.version()})"
            arch = platform.machine()
            py = platform.python_version()
        except Exception:
            os_line, arch, py = "неизвестно", "неизвестно", "неизвестно"

        try:
            eur = self._get_eur_sample_count()
            eur_text = "все доступные" if eur is None else str(eur)
        except Exception:
            eur_text = "неизвестно"

        lines = [
            f"Версия: {__version__}",
            f"ОС: {os_line}, {arch}",
            f"Python: {py}",
            f"Режим настроек: {self.mode_switch.get()}",
            f"Источник данных: {self.source_dd.get()}",
            f"Референсная панель: {self.panel_dd.get()}",
            f"Формат вывода: {self._get_format_key()}",
            f"Порог Rsq: {self.rsq_entry.get().strip()}",
            f"EUR-доноров: {eur_text}",
            f"Нормализация multiallelic: {'да' if self.normalize_var.get() else 'нет'}",
            f"Кэш сырых хромосом: {'да' if self.raw_cache_var.get() else 'нет'}",
            f"Переиспользование доноров: {'да' if self.reuse_donors_var.get() else 'нет'}",
            f"Активный запуск: {self.current_run_name or 'нет'}",
            f"Шаг мастера: {self._wizard_step}",
        ]
        return "\n".join(lines)

    def _collect_log_tail(self, max_lines: int = FEEDBACK_LOG_LINES) -> str:
        try:
            text = self.log_text.get("1.0", "end")
        except Exception:
            return ""
        lines = [line for line in text.splitlines() if line.strip()]
        return "\n".join(lines[-max_lines:])

    def _build_feedback_text(self, kind: str, subject: str, description: str,
                             include_log: bool) -> str:
        parts = [
            f"Тип обращения: {kind}",
            f"Тема: {subject}",
            "",
            "Описание:",
            description.strip() or "(не заполнено)",
            "",
            "--- Технические данные (заполнено программой) ---",
            self._collect_diagnostics(),
        ]
        if include_log:
            tail = self._collect_log_tail()
            parts += ["", f"--- Последние строки лога ---", tail or "(лог пуст)"]
        return "\n".join(parts)

    def _open_feedback_dialog(self):
        """
        Окно обратной связи. Письмо НЕ отправляется программой само:
        открывается почтовый клиент с уже заполненным письмом, и последнее
        слово — за пользователем. Так в дистрибутиве не появляется ни
        SMTP-пароля (его извлёк бы из exe любой желающий и разослал бы с
        этого ящика спам), ни отправки чего-либо за спиной пользователя —
        он видит текст письма целиком и может его отредактировать.
        """
        dialog = ctk.CTkToplevel(self)
        dialog.title("Обратная связь")
        dialog.geometry("820x660")
        dialog.minsize(700, 520)
        dialog.transient(self)
        # grab_set() до появления окна на экране на Windows иногда падает —
        # ставим модальность следующим тиком.
        dialog.after(200, lambda: dialog.grab_set())

        # Кнопки живут в ЗАКРЕПЛЁННОЙ нижней панели, а не внутри
        # прокручиваемой области: раньше они прокручивались вместе с
        # содержимым и не помещались по ширине — «Закрыть» упиралась в
        # край окна и обрезалась. Панель пакуется первой (side="bottom"),
        # иначе прокрутка с expand=True забрала бы всё место себе.
        buttons = ctk.CTkFrame(dialog, fg_color="transparent")
        buttons.pack(side="bottom", fill="x", padx=15, pady=(0, 15))

        frame = ctk.CTkScrollableFrame(dialog)
        frame.pack(fill="both", expand=True, padx=15, pady=(15, 10))

        ctk.CTkLabel(
            frame, text="Сообщить об ошибке или предложить улучшение",
            font=ctk.CTkFont(size=18, weight="bold"),
        ).pack(anchor="w", pady=(0, 8))

        ctk.CTkLabel(
            frame,
            text=(f"Письмо уйдёт на {FEEDBACK_EMAIL}. Программа сама ничего не "
                  f"отправляет: она откроет вашу почтовую программу с готовым "
                  f"письмом — перед отправкой его можно прочитать и поправить."),
            justify="left", text_color="gray60", wraplength=620,
        ).pack(anchor="w", pady=(0, 12))

        kind_var = ctk.StringVar(value=FEEDBACK_KINDS[0])
        ctk.CTkSegmentedButton(
            frame, values=list(FEEDBACK_KINDS), variable=kind_var, width=340,
        ).pack(anchor="w", pady=(0, 12))

        ctk.CTkLabel(frame, text="Коротко о чём (попадёт в тему письма):").pack(anchor="w")
        subject_tf = ctk.CTkEntry(frame, width=640,
                                  placeholder_text="например: не скачиваются доноры chr14")
        subject_tf.pack(fill="x", pady=(0, 12))
        attach_input_features(subject_tf)

        ctk.CTkLabel(frame, text="Подробности — что делали и что произошло:").pack(anchor="w")
        body_tf = ctk.CTkTextbox(frame, height=170)
        body_tf.pack(fill="x", pady=(0, 12))
        attach_input_features(body_tf)

        log_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            frame, text=f"Приложить последние {FEEDBACK_LOG_LINES} строк лога",
            variable=log_var,
        ).pack(anchor="w")
        ctk.CTkLabel(
            frame,
            text=("ℹ С логом разбираться в проблеме заметно проще. Учтите: в "
                  "нём встречаются пути к файлам, а значит имя пользователя "
                  "Windows и имя вашего файла с ДНК-данными. Письмо перед "
                  "отправкой открывается в почтовой программе — лишнее можно "
                  "удалить прямо там."),
            justify="left", text_color="gray60", wraplength=620,
        ).pack(anchor="w", pady=(2, 12))

        status_lbl = ctk.CTkLabel(frame, text="", justify="left",
                                  text_color="gray60", wraplength=620)
        status_lbl.pack(anchor="w", pady=(0, 10))

        def _texts():
            subject = subject_tf.get().strip() or "без темы"
            full_subject = f"{FEEDBACK_SUBJECT_PREFIX} v{__version__} — {kind_var.get()}: {subject}"
            body = self._build_feedback_text(
                kind_var.get(), subject, body_tf.get("1.0", "end"), log_var.get(),
            )
            return full_subject, body

        def _mailto_url(subject: str, body: str) -> str:
            return (f"mailto:{FEEDBACK_EMAIL}"
                    f"?subject={quote(subject)}&body={quote(body)}")

        def _on_mail():
            full_subject, body = _texts()
            saved_note = ""
            url = _mailto_url(full_subject, body)

            # Слишком длинное письмо почтовый клиент может обрезать на
            # полуслове. Практически это случается только с приложенным
            # логом: он в разы длиннее всего остального. Тогда полный текст
            # сохраняем файлом рядом с программой, а в письме оставляем всё
            # то же самое, кроме лога, и просим приложить файл вложением.
            if len(url) > FEEDBACK_MAILTO_URL_LIMIT:
                path = self._save_feedback_file(body)
                if path is not None:
                    saved_note = (f"\n\nЛог не поместился в письмо, поэтому "
                                  f"полный текст обращения вместе с ним сохранён "
                                  f"отдельным файлом:\n{path}\n"
                                  f"Приложите его к письму вложением.")
                short_body = self._build_feedback_text(
                    kind_var.get(), subject_tf.get().strip() or "без темы",
                    body_tf.get("1.0", "end"), include_log=False,
                ) + saved_note
                url = _mailto_url(full_subject, short_body)
                # Даже без лога не влезло — значит очень длинное описание.
                # Тогда в письме остаётся суть и путь к файлу с полным
                # текстом. Никаких циклов подгонки: приписка про файл сама
                # по себе длинная, и цикл «отрезать половину» на ней
                # никогда бы не сошёлся (на этом приложение уже висло).
                if len(url) > FEEDBACK_MAILTO_URL_LIMIT:
                    short_body = (
                        f"Тип обращения: {kind_var.get()}\n"
                        f"Тема: {subject_tf.get().strip() or 'без темы'}\n\n"
                        f"Описание оказалось слишком длинным для письма."
                        f"{saved_note}"
                    )
                    url = _mailto_url(full_subject, short_body)

            try:
                webbrowser.open(url)
            except Exception as e:
                status_lbl.configure(
                    text=f"⚠ Не удалось открыть почтовую программу: {e}\n"
                         f"Воспользуйтесь кнопкой «Скопировать текст».",
                    text_color="#F9A825",
                )
                return
            status_lbl.configure(
                text=("✓ Почтовая программа открыта — проверьте письмо и нажмите "
                      "«Отправить». Если окно не появилось, почтовый клиент не "
                      "настроен: скопируйте текст кнопкой ниже." + saved_note),
                text_color="#4CAF50" if not saved_note else "#F9A825",
            )

        def _on_copy():
            full_subject, body = _texts()
            self.clipboard_clear()
            self.clipboard_append(f"Кому: {FEEDBACK_EMAIL}\nТема: {full_subject}\n\n{body}")
            self.update()
            status_lbl.configure(
                text="✓ Текст письма скопирован — вставьте его в почту вручную.",
                text_color="#4CAF50",
            )

        def _on_save():
            _, body = _texts()
            path = self._save_feedback_file(body, ask=True)
            if path is not None:
                status_lbl.configure(text=f"✓ Сохранено: {path}", text_color="#4CAF50")

        # «Закрыть» пакуется ПЕРВОЙ и прижимается вправо: при нехватке
        # ширины pack обделяет тех, кто добавлен позже, — так ужмутся
        # вспомогательные кнопки, а выход из окна останется на месте.
        ctk.CTkButton(buttons, text="Закрыть", width=110,
                      fg_color="transparent", border_width=1,
                      command=dialog.destroy).pack(side="right", padx=(10, 0))
        ctk.CTkButton(buttons, text="✉ Открыть в почте", width=190,
                      command=_on_mail).pack(side="left", padx=(0, 8))
        ctk.CTkButton(buttons, text="📋 Скопировать", width=160,
                      fg_color="transparent", border_width=1,
                      command=_on_copy).pack(side="left", padx=(0, 8))
        ctk.CTkButton(buttons, text="💾 В файл", width=140,
                      fg_color="transparent", border_width=1,
                      command=_on_save).pack(side="left")

        # Esc закрывает окно — привычнее, чем искать кнопку.
        dialog.bind("<Escape>", lambda e: dialog.destroy())
        subject_tf.focus_set()

    def _save_feedback_file(self, body: str, ask: bool = False) -> Path | None:
        """Сохраняет текст обращения в файл. ask=True — со стандартным
        диалогом сохранения, иначе молча рядом с программой."""
        default_name = f"helixfilldna_feedback_{datetime.now():%Y%m%d_%H%M%S}.txt"
        if ask:
            chosen = filedialog.asksaveasfilename(
                defaultextension=".txt", initialfile=default_name,
                filetypes=[("Текстовый файл", "*.txt")],
            )
            if not chosen:
                return None
            target = Path(chosen)
        else:
            target = PROJECT_ROOT / default_name
        try:
            target.write_text(body, encoding="utf-8")
        except OSError:
            return None
        return target

    def _peek_log_line(self, msg: str):
        """
        Дублирует последнюю содержательную строку лога под шкалой на
        вкладке "Запуск". Разделители («====») и пустые строки
        пропускаются — они бы просто гасили полезное сообщение.
        """
        text = msg.strip()
        if not text or set(text) <= set("=-— "):
            return
        if len(text) > 140:
            text = text[:137] + "..."
        color = "gray60"
        if text.startswith(("✓", "✅")):
            color = "#4CAF50"
        elif text.startswith(("✗", "❌", "ОШИБКА", "⚠")):
            color = "#F44336"
        self.last_log_lbl.configure(text=text, text_color=color)

    def _upsert_progress_line(self, key: str, text: str):
        """
        Показывает text как ОДНУ строку в логе для данного key — повторный
        вызов с тем же key заменяет её содержимое на месте (через именованные
        метки Tkinter Text), а не добавляет новую строку. Первый вызов для
        нового key создаёт строку в конце лога и запоминает её границы.
        """
        start_mark, end_mark = f"prog_start_{key}", f"prog_end_{key}"
        if key in self._progress_keys:
            try:
                start = self.log_text.index(start_mark)
                end = self.log_text.index(end_mark)
                self.log_text.delete(start, end)
                self.log_text.insert(start, text, "progress")
                self.log_text.mark_set(end_mark, f"{start}+{len(text)}c")
                return
            except tk.TclError:
                # Метки потерялись (например, лог был очищен кнопкой
                # "Очистить лог") — создаём строку заново, как для нового key.
                self._progress_keys.discard(key)

        self.log_text.insert("end", text + "\n", "progress")
        end_index = self.log_text.index("end-1c")  # сразу перед добавленным \n
        start_index = f"{end_index}-{len(text)}c"
        self.log_text.mark_set(start_mark, start_index)
        self.log_text.mark_set(end_mark, end_index)
        self.log_text.mark_gravity(start_mark, "left")
        self._progress_keys.add(key)

    def _poll_logs(self):
        try:
            while True:
                msg = self.log_q.get_nowait()
                progress = _parse_progress_line(msg)
                if progress is not None:
                    key, text = progress
                    self._upsert_progress_line(key, text)
                elif any(msg.startswith(p) for p in ("✓", "✅")):
                    self.log_text.insert("end", msg + "\n", "success")
                elif any(msg.startswith(p) for p in ("✗", "❌", "ОШИБКА")):
                    self.log_text.insert("end", msg + "\n", "error")
                elif msg.startswith("["):
                    self.log_text.insert("end", msg + "\n", "stage")
                else:
                    self.log_text.insert("end", msg + "\n")
                self.log_text.see("end")
                self._peek_log_line(msg)
                self._update_donor_panel(msg)
        except queue.Empty:
            pass
        self.after(100, self._poll_logs)

    def _validate_settings(self) -> str | None:
        if not self.input_tf.get():
            return "Выберите файл с данными"
        if not self.tmpl_tf.get():
            if self._is_simple_mode():
                fmt = SIMPLE_FORMAT_BY_SOURCE.get(self._get_source_key(), "v3")
                expected = _samples_dir() / SAMPLE_TEMPLATE_NAMES.get(fmt, "template_v3.txt")
                return (
                    f"Трафарет для формата {fmt} не найден. Положите файл в\n"
                    f"{expected}\n"
                    f"или переключитесь в «{MODE_ADVANCED}» и укажите путь вручную."
                )
            return "Выберите трафарет"
        if not self._validate_eur_count_entry():
            return "Укажите корректное число EUR-доноров (или включите галочку 'все доступные')"
        return None

    def _check_donors(
        self, chip_signature: str, source: str | None = None, panel: str | None = None,
    ) -> list[Path]:
        """
        Тонкая обёртка над pipeline.check_donor_cache() (Задача A/B) —
        единая точка правды для CLI (main.main()) и GUI, без дублирования
        логики. Проверяет donors/<source>/<panel>/, сравнивает
        chip_signature.txt с сигнатурой текущего запуска и НЕ имеет права
        её перезаписывать. При отсутствии/несовпадении кэша бросает
        RuntimeError с понятной инструкцией (--source, --panel,
        --donors-subdir).

        Задача 2: source необязательный параметр — раньше всегда читался
        из выпадающего списка через self._get_source_key(), что создавало
        гонку после появления переключения источника в середине
        _run_stages_1_6() (обновление виджета уходит через self.after(0, ...)
        асинхронно). Явный параметр устраняет гонку.

        Шаг 1 промта "HRC / TopMed": аналогичный необязательный параметр
        `panel` — по умолчанию читается из self.panel_dd через
        self._get_panel_key(), но вызывающий код (_ensure_donors) передаёт
        его явно по тем же причинам, что и source.
        """
        if source is None:
            source = self._get_source_key()
        if panel is None:
            panel = self._get_panel_key()
        donors_root = PROJECT_ROOT / "donors"
        return pipeline.check_donor_cache(chip_signature, source, donors_root, panel=panel)

    # -----------------------------------------------------------------------
    # Задача 1: автопредложение скачать доноров через GUI
    # -----------------------------------------------------------------------
    def _prompt_yes_no(self, title: str, message: str) -> bool:
        """
        Thread-safe запрос messagebox.askyesno() из ФОНОВОГО потока.

        messagebox нельзя вызывать напрямую не из главного потока Tkinter —
        поэтому диалог показывается через self.after(0, ...), а этот метод
        (вызываемый из фонового потока) блокируется на threading.Event до
        тех пор, пока пользователь не ответит. Сам GUI при этом остаётся
        отзывчивым, т.к. блокируется только фоновый поток, а не mainloop.
        """
        response_event = threading.Event()
        response_holder = {"value": False}

        def ask_dialog():
            response_holder["value"] = messagebox.askyesno(title, message)
            response_event.set()

        self.after(0, ask_dialog)
        response_event.wait()
        return response_holder["value"]

    def _prompt_file_download_retry(self, filename: str, error: str) -> bool:
        """
        Thread-safe запрос "Повторить скачивание файла?" — вызывается из
        core.mis_adapter.MISAdapter.download_results() (через
        pipeline.download_mis_results_smart(..., on_file_error=...)) при
        неудаче скачивания КОНКРЕТНОГО ZIP-архива результатов MIS на
        Этапе 7. Тот же паттерн синхронизации, что и в _prompt_yes_no().

        При ответе "Да" download_results() повторяет попытку именно для
        этого файла (остальные уже скачанные файлы не трогаются). При
        "Нет" файл добавляется в список неудавшихся, но скачивание
        остальных файлов продолжается — MISAdapterError со списком всех
        проблемных файлов будет показана только в конце, после попытки
        скачать всё остальное.
        """
        return self._prompt_yes_no(
            "Ошибка скачивания файла",
            f"Не удалось скачать файл «{filename}»:\n\n{error}\n\n"
            f"Повторить попытку скачивания именно этого файла?",
        )

    def _prompt_info(self, title: str, message: str) -> None:
        """
        Thread-safe messagebox.showinfo() из ФОНОВОГО потока — тот же
        паттерн синхронизации, что и в _prompt_yes_no(), но без выбора:
        используется там, где нужно просто ПОКАЗАТЬ и дождаться, чтобы
        пользователь прочитал разъяснение (Шаг 2 промта "Доноры для
        VCF-источника: понятная отмена + общий кэш сырых хромосом"),
        прежде чем фоновый поток продолжит (обычно — бросит исключение,
        прерывающее запуск).
        """
        response_event = threading.Event()

        def show_dialog():
            messagebox.showinfo(title, message)
            response_event.set()

        self.after(0, show_dialog)
        response_event.wait()

    # -----------------------------------------------------------------------
    # Задача 2: диалог при несоответствии автодетекта и выбранного источника
    # -----------------------------------------------------------------------
    def _prompt_source_mismatch(self, detected: str, selected: str, confidence: float) -> str:
        """
        Thread-safe диалог с ТРЕМЯ вариантами при несовпадении автодетекта
        источника (pipeline.detect_source_from_file()) с тем, что выбрано
        в выпадающем списке GUI:
          "switch"   — сменить источник на detected
          "continue" — продолжить с уже выбранным источником
          "cancel"   — прервать запуск

        messagebox.askyesno() (как в _prompt_yes_no()) поддерживает только
        2 осмысленных исхода — здесь нужен третий ("Отмена"), поэтому
        диалог собран как модальное окно CTkToplevel с тремя кнопками.
        Синхронизация с фоновым потоком — тот же паттерн, что и в
        _prompt_yes_no(): показ окна уходит через self.after(0, ...),
        вызывающий (фоновый) поток блокируется на threading.Event.wait()
        до нажатия одной из кнопок или закрытия окна крестиком (это тоже
        трактуется как "cancel").
        """
        response_event = threading.Event()
        response_holder = {"value": "cancel"}

        def show_dialog():
            detected_name = pipeline.SOURCES.get(detected, {}).get("name", detected)
            selected_name = pipeline.SOURCES.get(selected, {}).get("name", selected)

            dialog = ctk.CTkToplevel(self)
            dialog.title("Возможное несоответствие источника")
            dialog.geometry("480x260")
            dialog.resizable(False, False)
            dialog.transient(self)
            dialog.grab_set()

            msg = (
                f"Похоже, выбранный файл на самом деле в формате «{detected_name}» "
                f"(уверенность {confidence:.0%}), а в настройках выбран источник "
                f"«{selected_name}».\n\n"
                f"Что сделать?"
            )
            ctk.CTkLabel(dialog, text=msg, wraplength=440, justify="left").pack(
                padx=20, pady=(20, 15), fill="x"
            )

            btn_col = ctk.CTkFrame(dialog, fg_color="transparent")
            btn_col.pack(fill="x", padx=20, pady=(0, 20))

            def choose(value: str):
                response_holder["value"] = value
                response_event.set()
                dialog.destroy()

            ctk.CTkButton(
                btn_col, text=f"Сменить на «{detected_name}»",
                fg_color="#2E7D32", hover_color="#1B5E20",
                command=lambda: choose("switch"),
            ).pack(fill="x", pady=(0, 8))
            ctk.CTkButton(
                btn_col, text=f"Продолжить с «{selected_name}»",
                command=lambda: choose("continue"),
            ).pack(fill="x", pady=(0, 8))
            ctk.CTkButton(
                btn_col, text="Отмена", fg_color="#B71C1C", hover_color="#7F0000",
                command=lambda: choose("cancel"),
            ).pack(fill="x")

            dialog.protocol("WM_DELETE_WINDOW", lambda: choose("cancel"))

        self.after(0, show_dialog)
        response_event.wait()
        return response_holder["value"]

    def _show_cancel_donor_btn(self):
        """Показывает кнопку остановки — только на том этапе, где она
        реально работает (скачивание доноров). В остальное время её на
        экране нет вовсе, а не «есть, но серая»."""
        self.cancel_donor_btn.configure(
            state="normal", text="⏹ Остановить скачивание доноров",
        )
        self.stop_box.pack(fill="x", pady=(0, 20), after=self.start_btn)

    def _hide_cancel_donor_btn(self):
        self.stop_box.pack_forget()
        # Единственная точка, которую вызывают ВСЕ пути выхода из этапа 3
        # (успех, ошибка, отмена, finally) — заодно снимаем опрос файлов,
        # чтобы after() не тикал вхолостую после конца скачивания.
        self._stop_file_watch()

    def _on_cancel_donor_download(self):
        """
        Обработчик кнопки "Отменить скачивание доноров". Устанавливает
        threading.Event, который download_donors.download_donors_for_chip()
        проверяет перед каждой хромосомой И внутри самих сетевых операций
        (через _run_cancelable/Popen) — отмена прерывает уже идущую
        загрузку файла (terminate()/kill() дочернего curl), а не только
        ждёт завершения текущей хромосомы.
        """
        self._cancel_donor_download.set()
        self.cancel_donor_btn.configure(state="disabled", text="⏳ Останавливаю...")
        print("⏳ Запрошена отмена скачивания доноров — завершаю текущую операцию...")

    def _ensure_donors(
        self, source: str, signature: str, positions_json: Path, panel: str,
    ) -> list[Path]:
        """
        Этап 3 с автопредложением скачать доноров (Задача 1).

        1. Пытается pipeline.check_donor_cache() как раньше.
        2. Если он бросает RuntimeError (доноров нет или сигнатура не
           совпадает) — спрашивает пользователя через messagebox.askyesno
           (thread-safe, см. _prompt_yes_no).
        3. При согласии — качает доноров ПРЯМО в этом же фоновом потоке
           через download_donors.download_donors_for_chip(), с прогрессом
           в _set_subprogress(3, ...) и поддержкой отмены/повтора.
        4. После успешного скачивания повторно вызывает
           pipeline.check_donor_cache() — единая точка правды, гарантирует,
           что chip_signature.txt реально записан и список путей финальный
           (так и задумано в докстринге download_donors_for_chip()).
        5. При отказе пользователя или отказе от повтора после
           отмены/ошибки — бросает RuntimeError, что штатно ловится общим
           except в _run_stages_1_6() (кнопки разблокируются как обычно).

        CLI (main.py::main()) эту функцию не использует и продолжает
        падать с RuntimeError как раньше — автопредложение только в GUI.

        Задача D: signature передаётся явно (строгая или широкая — в
        зависимости от чекбокса "Переиспользовать доноров..."), вместо
        того чтобы читаться из result.chip_signature внутри — иначе
        широкий режим здесь тихо игнорировался бы.

        Шаг 1 промта "HRC / TopMed": panel также передаётся явно (по той
        же причине, что и source, см. докстринг _check_donors) — папка
        для скачивания доноров теперь donors/<source>/<panel>/, а не
        donors/<source>/.
        """
        try:
            return self._check_donors(signature, source, panel)
        except RuntimeError as e:
            print(f"⚠ {e}")

        donors_root = PROJECT_ROOT / "donors"
        output_dir = pipeline._donor_source_dir(source, donors_root, panel=panel)

        # Промт "Доноры для VCF-источника: понятная отмена + общий кэш
        # сырых хромосом", Шаг 1: явно объясняем ДО показа диалога
        # "скачать?", почему доноры 1000 Genomes нужны ЛЮБОМУ источнику
        # (включая 'vcf', которому НЕ нужен референсный геном — это
        # единственное, о чём пользователю где-либо явно сообщалось, и
        # легко перепутать "не нужен референс" с "вообще ничего не
        # нужно"). Печатается в лог один раз перед диалогом, а не только
        # в самом диалоге — так объяснение остаётся видимым в run.log,
        # даже если пользователь не читал диалог внимательно.
        print(
            "ℹ Доноры 1000 Genomes нужны для ЛЮБОГО источника данных (FTDNA, "
            "MyHeritage, готовый VCF) — они не связаны с референсным геномом "
            "(которого 'vcf' действительно не требует) и служат для самого "
            "шага импутации на Michigan Imputation Server. Каждый новый "
            "источник данных или новый чип/набор позиций требует своего "
            "кэша доноров, отфильтрованного именно под эти позиции — "
            "поэтому для текущего источника/чипа кэш ещё не готов."
        )

        while True:
            proceed = self._prompt_yes_no(
                "Доноры отсутствуют или устарели",
                f"Донорские файлы для источника «{source}» "
                f"(панель «{pipeline.REFERENCE_PANELS[panel]['display_name']}») отсутствуют "
                f"или устарели (не подходят под текущий чип).\n\n"
                f"Доноры 1000 Genomes нужны для ЛЮБОГО источника (в том числе "
                f"«Готовый VCF») — они не связаны с референсным геномом "
                f"(который для VCF действительно не нужен), а нужны для "
                f"самого шага импутации на Michigan Imputation Server. Каждый "
                f"новый источник/чип требует своего кэша доноров.\n\n"
                f"Скачать ~22 файла доноров 1000 Genomes автоматически?\n"
                f"Трафик скачивания может составить несколько десятков ГБ "
                f"(итоговый размер отфильтрованных файлов на диске "
                f"меньше).",
            )
            if not proceed:
                # Шаг 2: прежде чем прерывать запуск, явно объясняем
                # последствия отказа — иначе человек, ответивший "Нет" не
                # разобравшись, видит только внезапное прерывание запуска.
                self._prompt_info(
                    "Запуск будет прерван",
                    "Без доноров сборка не может продолжиться — Этапы 4-6 "
                    "(объединение с донорами, разбивка по хромосомам) "
                    f"требуют donors/{source}/{panel}/.\n\n"
                    "Запуск будет прерван. Файлы, подготовленные до этого "
                    "шага (sample.vcf.gz), останутся в папке запуска, но "
                    "собрать 22 файла для загрузки на MIS не получится, пока "
                    "доноры не будут скачаны.",
                )
                raise UserCancelledRun(
                    "Скачивание доноров отменено пользователем — запуск прерван."
                )

            self._cancel_donor_download.clear()
            self.after(0, self._show_cancel_donor_btn)

            def donor_progress(frac: float, text: str) -> None:
                self.after(0, self._set_subprogress, 3, frac, text)

            # Промт "...общий кэш сырых хромосом", Шаг 5: общий кэш
            # ЕЩЁ НЕ отфильтрованных полных хромосом 1000 Genomes,
            # переиспользуемый между всеми источниками/чипами этой
            # референсной сборки (donors/_raw_chromosomes/<genome_build>/)
            # — включается только по явному желанию пользователя
            # (self.raw_cache_var), по умолчанию выключен.
            raw_cache_dir = (
                pipeline.raw_chromosome_cache_dir(donors_root, panel)
                if self.raw_cache_var.get() else None
            )

            try:
                # Промт "HRC / TopMed" (v4-фикс): genome_build обязателен —
                # без него download_donors_for_chip() молча использует
                # дефолт DEFAULT_GENOME_BUILD="grch37" независимо от
                # выбранной панели, и для panel="topmed" доноры качались бы
                # с GRCh37-зеркал 1000 Genomes вместо GRCh38 (GRCH38_MIRRORS).
                # Промт "не видно, какой файл качается и сколько мегабайт":
                # следим за размерами файлов в папке доноров (и в общем
                # кэше сырых хромосом, если он включён) — это единственный
                # источник, не зависящий от того, чем идёт закачка.
                self.after(0, self._start_file_watch,
                           [output_dir, raw_cache_dir], "donors")
                download_donors.download_donors_for_chip(
                    positions_json, source, output_dir, pipeline.HTSLIB,
                    progress_cb=donor_progress,
                    cancel_check=self._cancel_donor_download.is_set,
                    raw_cache_dir=raw_cache_dir,
                    genome_build=pipeline.REFERENCE_PANELS[panel]["genome_build"],
                    # Промт "Monomorphic sites / настраиваемое количество
                    # EUR-доноров": None (галочка "все доступные", по
                    # умолчанию) — вся EUR-подвыборка панели, иначе —
                    # явное число из self.eur_count_entry.
                    eur_sample_count=self._get_eur_sample_count(),
                )
            except download_donors.DownloadCancelled:
                self.after(0, self._hide_cancel_donor_btn)
                print("⚠ Скачивание доноров отменено пользователем.")
                if self._prompt_yes_no(
                    "Скачивание отменено",
                    "Скачивание доноров было отменено.\n\nПопробовать снова?",
                ):
                    continue
                raise UserCancelledRun(
                    "Скачивание доноров отменено пользователем — запуск прерван."
                )
            except RuntimeError as e:
                self.after(0, self._hide_cancel_donor_btn)
                print(f"❌ Ошибка скачивания доноров: {e}")
                if self._prompt_yes_no(
                    "Ошибка скачивания",
                    f"Не удалось скачать доноров:\n\n{e}\n\n"
                    f"Повторить попытку? (уже скачанные хромосомы не "
                    f"будут перекачаны заново — докачка продолжится с "
                    f"места обрыва).",
                ):
                    continue
                raise
            else:
                self.after(0, self._hide_cancel_donor_btn)
                # Единая точка правды (см. докстринг
                # download_donors_for_chip): убеждаемся, что сигнатура
                # реально записана, и получаем финальный список путей.
                return self._check_donors(signature, source, panel)

    # --- Задача 6: плавный прогресс этапов 1-6 ---------------------------
    def _set_stage(self, n: int, text: str):
        """Оставлен для обратной совместимости — эквивалент sub_progress=0."""
        self._set_subprogress(n, 0.0, text)

    def _set_subprogress(self, stage_n: int, sub_progress: float, text: str):
        """
        Плавный прогресс внутри этапа: общий прогресс = (stage_n - 1 + sub_progress) / STAGES_TOTAL.
        sub_progress — доля выполнения текущего этапа (0.0 .. 1.0).
        Вызывается только через self.after(0, ...) из фонового потока.

        Делитель — STAGES_TOTAL (6), а не 7: эта шкала показывает прогресс
        кнопки "Запустить этапы 1-6 (до MIS)", то есть ровно шести этапов.
        Раньше делитель был 7 (с учётом Этапа 7 — сборки после MIS), из-за
        чего по завершении подготовки шкала замирала на 6/7 (~86%) и
        выглядела недоделанной, хотя работа была закончена. Этап 7 живёт в
        отдельной секции вкладки и рисует эту же шкалу своим собственным
        _set_stage7_progress() от 0 до 1.
        """
        sub_progress = max(0.0, min(1.0, sub_progress))
        overall = (stage_n - 1 + sub_progress) / STAGES_TOTAL
        self.progress.set(overall)
        self._update_eta(overall)
        self.stage_lbl.configure(text=f"[{stage_n}/{STAGES_TOTAL}] {text}")
        print(f"[{stage_n}/{STAGES_TOTAL}] {text}")

    # --- Задача 7: прогресс + краткие сообщения этапа 7 -------------------
    def _set_stage7_progress(self, frac: float, text: str):
        """frac в диапазоне 0.0..1.0: 0-0.5 скачивание, 0.5-1.0 сборка."""
        frac = max(0.0, min(1.0, frac))
        self.progress.set(frac)
        self._update_eta(frac)
        self.stage_lbl.configure(text=text)
        print(text)

    # -----------------------------------------------------------------------
    # Обработчики кнопок
    # -----------------------------------------------------------------------
    def _on_start(self):
        if self.running:
            return
        err = self._validate_settings()
        if err:
            messagebox.showwarning("Предупреждение", err)
            return

        # Промт "Именованные папки запуска": резолвим папку ЭТОГО запуска
        # ДО старта фонового потока — коллизия имён (пользователь вручную
        # ввёл имя уже существующего запуска) должна быть явной ошибкой
        # здесь и сейчас, а не тихой перезаписью чужих файлов где-то в
        # середине _run_stages_1_6().
        run_name_input = self.run_name_tf.get().strip() or None
        try:
            run_dir, run_name = pipeline.resolve_run_dir(
                PROJECT_ROOT / "output", run_name_input, must_exist=False,
            )
        except RuntimeError as e:
            messagebox.showwarning("Предупреждение", str(e))
            return
        self._set_active_run(run_dir, run_name)

        self.running = True
        self._run_started_at = time.monotonic()
        self.start_btn.configure(state="disabled")
        self.progress.set(0)
        self.stage_lbl.configure(text="Запуск...")
        self.eta_lbl.configure(text="", text_color="gray60")
        self._reset_donor_panel()
        self._set_wizard_step(1)
        self._refresh_mis_btn_state()
        threading.Thread(target=self._run_stages_1_6, daemon=True).start()

    def _on_mis(self):
        if self.running:
            return
        if self.current_run_dir is None:
            messagebox.showwarning(
                "Предупреждение",
                "Не выбран активный запуск. Выполните Шаг 1, либо "
                "разверните «Подробнее», выберите завершённый запуск в "
                "списке и нажмите «▶ Продолжить (Шаг 3)».",
            )
            return
        curl = self.curl_tf.get("1.0", "end").strip()
        pwd = self.pwd_tf.get().strip()
        if not curl or not pwd:
            messagebox.showwarning("Предупреждение", "Укажите curl-команду и пароль")
            return
        if not self._validate_rsq_entry():
            messagebox.showwarning(
                "Предупреждение", "Порог Rsq должен быть числом от 0.30 до 0.99",
            )
            return
        # Баг-фикс: раньше pipeline.HTSLIB (пути к bcftools.exe/tabix.exe)
        # переустанавливался ТОЛЬКО внутри Этапа 1-6 (_run_stages_1_6) —
        # если пользователь закрывал и заново открывал GUI, а затем сразу
        # выбирал «▶ Продолжить (Этап 2)» и переходил к Этапу 7, минуя
        # Этап 1-6 в этой сессии, pipeline.HTSLIB оставался дефолтным
        # (bin_dir=None, из первоначального импорта main.py), и код искал
        # bcftools/tabix только в системном PATH — где их обычно нет, они
        # лежат в папке --bin-dir. Результат — "[WinError 2] Не удаётся
        # найти указанный файл" при попытке subprocess.run(["tabix", ...])
        # с голым именем без пути. Теперь HTSLIB инициализируется из
        # self.bin_tf здесь же, перед стартом Этапа 7 — тем же способом,
        # каким это уже делает _run_stages_1_6() перед Этапом 1-6.
        bd = Path(self.bin_tf.get()) if self.bin_tf.get() else None
        pipeline.HTSLIB = pipeline.HtslibTools(bd)
        if bd and str(bd) not in os.environ.get("PATH", ""):
            os.environ["PATH"] = str(bd) + os.pathsep + os.environ.get("PATH", "")
        self.running = True
        self._run_started_at = time.monotonic()
        self.progress.set(0)
        self.eta_lbl.configure(text="", text_color="gray60")
        self._set_wizard_step(3)
        self._refresh_mis_btn_state()
        threading.Thread(target=self._run_stage_7, daemon=True).start()

    # -----------------------------------------------------------------------
    # Этапы 1-6 (Задача 6: плавный прогресс внутри каждого этапа)
    # -----------------------------------------------------------------------
    def _run_stages_1_6(self):
        # Промт "Именованные папки запуска": self.current_run_dir уже
        # выставлен в _on_start() ДО запуска этого потока — печатаем
        # print()-сообщения и в очередь GUI, и в <run_dir>/run.log.
        # logger.*-сообщения дублируются отдельно, через FileHandler,
        # уже подключённый в _set_active_run()/pipeline.attach_run_log_handler().
        run_log_path = (self.current_run_dir / "run.log") if self.current_run_dir else None
        stdout_redirector = LogRedirector(self.log_q, run_log_path)
        stderr_redirector = LogRedirector(self.log_q, run_log_path)
        try:
            with contextlib.redirect_stdout(stdout_redirector), \
                 contextlib.redirect_stderr(stderr_redirector):

                if self.current_run_dir is None:
                    raise RuntimeError(
                        "Внутренняя ошибка: папка запуска не была выбрана "
                        "перед стартом Этапа 1-6."
                    )

                bd = Path(self.bin_tf.get()) if self.bin_tf.get() else None
                pipeline.HTSLIB = pipeline.HtslibTools(bd)
                if bd and str(bd) not in os.environ.get("PATH", ""):
                    os.environ["PATH"] = str(bd) + os.pathsep + os.environ.get("PATH", "")

                # v13 (промт "Диагностика + устойчивая настройка CA-сертификатов"):
                # ставим bin_dir в НАЧАЛО PATH (строка выше) — это удобно
                # для поиска bgzip/tabix/bcftools по имени, но на Windows
                # именно так в PATH может попасть собственный curl.exe из
                # бандла htslib, конфликтующий с системным (см. докстринг
                # core/network_utils.py). ensure_network_ready() ниже:
                # (1) гарантирует CA-сертификаты для libcurl в bcftools
                # программно (os.environ), без требования к пользователю
                # вручную выставлять $env:CURL_CA_BUNDLE в каждой сессии;
                # (2) предупреждает в лог, если конфликтующий curl.exe
                # обнаружен — сам обход конфликта (игнорирование bin_dir
                # при поиске curl) сделан внутри download_donors.py через
                # which_curl_ignoring_dir(), не здесь.
                network_utils.ensure_network_ready(bd)

                if not pipeline.HTSLIB.has_bcftools:
                    raise RuntimeError("bcftools не найден в указанной папке бинарников")

                source = self._get_source_key()
                panel = self._get_panel_key()
                panel_cfg = pipeline.REFERENCE_PANELS[panel]
                print(f"ℹ Референсная панель: {panel_cfg['display_name']}")

                csv_path = Path(self.input_tf.get())
                if not csv_path.exists():
                    raise RuntimeError(f"Файл с данными не найден: {csv_path}")
                # Промт "Именованные папки запуска": все артефакты этого
                # запуска (sample.vcf.gz, upload/, parse_result.pkl и т.д.)
                # пишутся в СВОЮ папку output/runs/<run_name>/, а не в общую
                # output/ — donors/<source>/<panel>/ остаётся общим кэшем,
                # см. pipeline._donor_source_dir() ниже, без изменений.
                output_dir = self.current_run_dir
                output_dir.mkdir(parents=True, exist_ok=True)
                print(f"ℹ Папка запуска: {output_dir} (имя запуска: {self.current_run_name!r})")

                # --- Задача 2: автодетект источника vs выбранный в GUI --
                # Вызывается СРАЗУ после проверки, что csv_path существует,
                # и ДО "[0/7] Проверка референсного генома" — та может
                # качать/проверять несколько ГБ, а источнику 'vcf' референс
                # вообще не нужен, так что не тратим на это время до
                # проверки соответствия источника файлу.
                detected_source, confidence = pipeline.detect_source_from_file(csv_path)
                if detected_source:
                    print(
                        f"ℹ Автодетект: файл похож на {detected_source} "
                        f"(уверенность {confidence:.2f}), выбрано {source}"
                    )
                    if confidence >= 0.8 and detected_source != source:
                        choice = self._prompt_source_mismatch(detected_source, source, confidence)
                        if choice == "cancel":
                            raise UserCancelledRun(
                                "Запуск отменён пользователем: обнаружено "
                                "несоответствие источника и формата файла."
                            )
                        if choice == "switch":
                            new_name = pipeline.SOURCES[detected_source]["name"]
                            self.after(0, self.source_dd.set, new_name)
                            # CTkOptionMenu.set() не вызывает command=,
                            # поэтому в обычном режиме пресеты (формат
                            # вывода/трафарет) пересчитываем явно.
                            self.after(0, self._on_source_changed)
                            source = detected_source
                            print(f"✓ Источник переключён на: {new_name}")
                        # choice == "continue" — ничего не меняем, source
                        # остаётся тем, что выбрал пользователь изначально.
                else:
                    print(f"ℹ Автодетект: не удалось определить формат файла {csv_path}")

                # Промт "Именованные папки запуска" (run_info.json,
                # Шаг 4): фиксируем метаданные запуска, как только source/
                # panel окончательно определены (после возможного
                # переключения из-за автодетекта выше).
                pipeline.save_run_info(
                    output_dir,
                    run_name=self.current_run_name,
                    started_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    source=source,
                    panel=panel,
                    csv_filename=csv_path.name,
                    format=self._get_format_key(),
                    rsq_threshold=self._get_rsq_threshold(),
                    normalize=self.normalize_var.get(),
                    reuse_donors_across_people=self.reuse_donors_var.get(),
                    raw_chromosome_cache=self.raw_cache_var.get(),
                    # save_run_info() игнорирует поля со значением None —
                    # передаём "all" вместо None, иначе выбор "все доступные
                    # EUR-доноры" вообще не попал бы в run_info.json/историю
                    # запусков.
                    eur_sample_count=self._get_eur_sample_count() or "all",
                )

                # --- Этап 0: приведение файла к оформлению 23andMe v3 ---
                # Только для источников из
                # pipeline._SOURCES_NEEDING_CONVERSION (сейчас это
                # 'ancestry'); для остальных prepare_source_file() —
                # no-op, возвращающий тот же путь, поэтому FTDNA/
                # MyHeritage/VCF этот блок никак не задевает.
                #
                # Стоит ДО "[0/7] Проверка референсного генома" по той же
                # причине, что и автодетект источника выше: проверка
                # референса может качать и хешировать гигабайты, и если
                # файл окажется не того формата, узнать об этом лучше до
                # неё, а не после.
                #
                # Конвертированный файл кладём в папку запуска — он
                # самостоятельный результат: его можно проверить глазами
                # и загрузить в Генотек как есть, не дожидаясь импутации.
                csv_for_parsing = csv_path
                if source in pipeline._SOURCES_NEEDING_CONVERSION:
                    self.after(0, self._set_subprogress, 1, 0.0,
                               "Приведение к формату 23andMe v3...")
                    csv_for_parsing, conversion_stats = pipeline.prepare_source_file(
                        source, csv_path, output_dir, Path(self.tmpl_tf.get()),
                    )
                    print(f"✓ Этап 0: {conversion_stats.summary()}")
                    if not conversion_stats.skipped:
                        pipeline.save_run_info(
                            output_dir,
                            converted_file=Path(conversion_stats.out_path).name,
                        )

                print("[0/7] Проверка референсного генома")
                reference = None
                if source != "vcf":
                    ref_path = pipeline.ensure_reference_genome(PROJECT_ROOT, panel=panel)
                    reference = ReferenceGenome(ref_path)

                # --- Лифтовер вперёд (промт "встроить лифтовер HRC/TopMed в
                # gui/app.py", аналог Этапа 2.5) ---------------------------
                # В этой архитектуре лифтовер точечный — происходит ВНУТРИ
                # parser_fn (adapters/ftdna_v3.py::parse_ftdna_v3()/
                # adapters/myheritage_v5.py::parse_myheritage_v5()), а не
                # отдельным проходом над уже готовым VCF, поэтому chain-файл
                # нужно подготовить ДО вызова parser_fn, но ПОСЛЕ референса
                # (использует тот же project_root/панель). Для panel="hrc"
                # pipeline._build_liftover() всегда возвращает None (в
                # REFERENCE_PANELS["hrc"] нет liftover_chain_url) — поведение
                # для HRC не меняется ни на йоту, liftover=None прокидывается
                # в parser_fn точно как раньше.
                liftover = None
                if pipeline._supports_liftover(source):
                    self.after(0, self._set_subprogress, 1, 0.0, "Проверка chain-файла лифтовера...")
                    liftover = pipeline._build_liftover(
                        panel, project_root=PROJECT_ROOT,
                        progress_cb=lambda p, t: self.after(0, self._set_subprogress, 1, 0.05 * p, t),
                    )
                    if panel != pipeline.DEFAULT_PANEL and liftover is None:
                        print(
                            f"⚠ Не удалось построить лифтовер для панели "
                            f"'{panel}' (chain-файл отсутствует в конфигурации) "
                            f"— координаты останутся в GRCh37, результат "
                            f"импутации, скорее всего, будет некорректным."
                        )
                elif panel != pipeline.DEFAULT_PANEL:
                    # Единственное реальное ограничение (source='vcf' не
                    # поддерживает лифтовер, см. pipeline._supports_liftover())
                    # — предупреждаем здесь один раз, где уже известен
                    # финальный source (после возможного переключения по
                    # автодетекту выше), а не в _on_panel_changed().
                    print(
                        f"⚠ Источник '{source}' пока не поддерживает лифтовер "
                        f"координат (см. pipeline._supports_liftover()) — если "
                        f"исходные координаты чипа не в сборке "
                        f"{panel_cfg['genome_build'].upper()}, результат "
                        f"импутации будет некорректным."
                    )

                # --- Этап 1: Парсинг --------------------------------------
                self.after(0, self._set_subprogress, 1, 0.05, "Чтение файла...")
                parser_fn = pipeline.SOURCES[source]["parser"]
                # csv_for_parsing — результат Этапа 0 (для 'ancestry')
                # либо сам csv_path (для всех остальных источников).
                if pipeline._supports_liftover(source):
                    result = parser_fn(csv_for_parsing, reference, liftover=liftover)
                else:
                    result = parser_fn(csv_for_parsing, reference)
                self.after(0, self._set_subprogress, 1, 0.5, "Разрешение ориентации...")
                print(f"Годных вариантов: {len(result.variants)}, сигнатура: {result.chip_signature}")
                if getattr(result, "lift_failed", 0):
                    print(
                        f"⚠ Не перенесено лифтовером на целевую сборку "
                        f"(lift_failed): {result.lift_failed}"
                    )

                # Задача B (доп. пункт): сохраняем позиции чипа тем же
                # парсером, который уже отработал — download_donors.py
                # потом читает этот JSON через --positions-json вместо
                # повторного (и потенциально неверного для MyHeritage/VCF)
                # парсинга CSV.
                # Задача D: строгая/широкая сигнатура — единая логика с CLI.
                signature, save_pos_fn = pipeline._resolve_chip_signature_mode(
                    result, source,
                    reuse_donors_across_people=self.reuse_donors_var.get(),
                )
                # Шаг 1 промта "HRC / TopMed": позиции чипа теперь тоже
                # хранятся раздельно по панели (donors/<source>/<panel>/) —
                # см. pipeline._donor_source_dir().
                positions_cache_dir = pipeline._donor_source_dir(
                    source, PROJECT_ROOT / "donors", panel=panel,
                )
                positions_json = save_pos_fn(positions_cache_dir, result)
                print(f"Позиции чипа сохранены: {positions_json}")

                pipeline.save_run_info(
                    output_dir,
                    chip_signature=result.chip_signature,
                    chip_signature_broad=getattr(result, "chip_signature_broad", None),
                )

                # ВАЖНО (Задача A): здесь больше НЕТ вызова
                # pipeline._save_chip_signature() — сигнатура пишется
                # только в download_donors.py после свежего скачивания
                # доноров. Запись здесь, ДО проверки кэша на Этапе 3,
                # делала сравнение сигнатур бессмысленным (кэш от другого
                # чипа тихо принимался как валидный).
                with (output_dir / "parse_result.pkl").open("wb") as f:
                    pickle.dump(result, f)
                self.after(0, self._set_subprogress, 1, 1.0, "Парсинг завершён")

                # --- Этап 2: VCF -------------------------------------------
                self.after(0, self._set_subprogress, 2, 0.0, "Сборка VCF...")
                sample_vcf = output_dir / "sample.vcf.gz"
                # Промт "HRC / TopMed" (v4-фикс): chrom_prefix обязателен —
                # без него CHROM в sample.vcf.gz для panel="topmed" писался
                # бы без "chr", и bcftools merge с GRCh38-донорами (Этап 5)
                # молча не находил бы пересечений по CHROM.
                pipeline.build_vcf(result, sample_vcf, sample_name="genotek",
                                   bgzip_path=pipeline.HTSLIB.bgzip_path,
                                   chrom_prefix=panel_cfg["chrom_prefix"])
                self.after(0, self._set_subprogress, 2, 0.7, "Индексация...")
                pipeline._index_vcf(sample_vcf)
                self.after(0, self._set_subprogress, 2, 1.0, "VCF готов")

                # --- Этап 3: Кэш доноров (Задача 1: автопредложение) ------
                self.after(0, self._set_subprogress, 3, 0.0, "Проверка кэша доноров...")
                donor_vcfs = self._ensure_donors(source, signature, positions_json, panel)
                self.after(0, self._set_subprogress, 3, 1.0, "Доноры готовы")

                # --- Этап 4: Concat -----------------------------------------
                self.after(0, self._set_subprogress, 4, 0.0, "Объединение 22 доноров...")
                kgp_all = output_dir / "kgp_all.vcf.gz"
                pipeline._concat_donors(donor_vcfs, kgp_all)
                self.after(0, self._set_subprogress, 4, 1.0, "Готово")

                # --- Этап 5: Merge -------------------------------------------
                self.after(0, self._set_subprogress, 5, 0.0, "Подготовка к merge...")
                merged = output_dir / "batch_merged.vcf.gz"
                self.after(0, self._set_subprogress, 5, 0.3, "Проверка индексов...")
                self.after(0, self._set_subprogress, 5, 0.7, "Merge sample + доноры...")
                pipeline._merge_with_donors_bcftools(sample_vcf, kgp_all, merged)
                self.after(0, self._set_subprogress, 5, 1.0, "Готово")

                # --- Этап 6: (опц.) нормализация + post-merge intersect + Split ----
                self.after(0, self._set_subprogress, 6, 0.0, "Проверка позиций...")

                # Задача C, необязательный чекбокс: bcftools norm -m-both.
                # НЕ входит в критический путь фикса Invalid alleles
                # (это делают Задачи A/B) — отдельная оптимизация под
                # multiallelic-сайты, выключена по умолчанию.
                if self.normalize_var.get():
                    if reference is None:
                        print(
                            f"⚠ Нормализация запрошена, но у источника '{source}' "
                            f"нет референса — пропускаю"
                        )
                    else:
                        self.after(0, self._set_subprogress, 6, 0.15,
                                   "Нормализация (bcftools norm -m-both)...")
                        normalized = output_dir / "batch_merged.norm.vcf.gz"
                        merged = pipeline._normalize_vcf(
                            merged, reference.fasta_path, normalized,
                            pipeline.HTSLIB.bcftools_path,
                        )

                # Задача C: post-merge intersect как диагностический/
                # защитный слой. После исправлений Задач A/B должен быть
                # no-op (0 удалено). Если удаляется больше 0 — это сигнал
                # регрессии (кэш доноров не соответствует текущему чипу),
                # выводим явное предупреждение в лог, но не прерываем сборку —
                # это диагностика, а не основной фикс.
                self.after(0, self._set_subprogress, 6, 0.3,
                           "Post-merge intersect (обязательная фильтрация)...")
                checked = output_dir / "batch_merged.checked.vcf.gz"
                try:
                    merged, before_n, after_n = pipeline._post_merge_intersect(
                        merged, donor_vcfs, checked, pipeline.HTSLIB.bcftools_path,
                        kgp_all_vcf=kgp_all,
                    )
                    removed_n = before_n - after_n
                    if removed_n > 0:
                        print(
                            f"✓ Post-merge intersect: {before_n} → {after_n} позиций "
                            f"({removed_n} удалено — позиции, отсутствующие в донорской "
                            f"подвыборке, корректно отфильтрованы)."
                        )
                    else:
                        print(
                            f"⚠ Post-merge intersect: 0 удалено ({before_n} → {after_n}). "
                            f"Проверьте итоговый call rate после Этапа 7 — для "
                            f"panel != 'hrc' 0 удалённых позиций нетипично."
                        )
                except (pipeline.PureCoreError, subprocess.CalledProcessError) as e:
                    # ⚠ Больше НЕ диагностика — сбой этого шага фатален (см.
                    # докстринг pipeline._post_merge_intersect()): без
                    # фильтрации по донорским позициям в отправку на сервер
                    # импутации попадут десятки/сотни тысяч позиций,
                    # физически отсутствующих в донорской подвыборке 1000
                    # Genomes — сервер гарантированно провалит QC.
                    # Пробрасываем как обычное исключение — подхватывается
                    # уже существующим `except Exception as e:` ниже в этом
                    # методе, которое покажет ошибку в логе GUI, остановит
                    # выполнение ДО Этапа 6/загрузки на MIS и разблокирует
                    # кнопку «Запустить».
                    raise RuntimeError(
                        f"Post-merge intersect не удался: {e}\n\n"
                        f"Это НЕ безобидная диагностика: без фильтрации по донорским "
                        f"позициям в отправку на сервер импутации попадут десятки/сотни "
                        f"тысяч позиций, которых физически нет в донорской подвыборке "
                        f"1000 Genomes — сервер (Michigan Imputation Server / BioData "
                        f"Catalyst) гарантированно провалит QC ('Invalid Alleles', "
                        f"'SNPs call rate < 90%', 'No chunks passed the QC step'). "
                        f"Запуск остановлен до создания 22 файлов для загрузки, чтобы "
                        f"не тратить время на заведомо провальное задание."
                    ) from e

                self.after(0, self._set_subprogress, 6, 0.5, "Разбивка по хромосомам...")
                upload_dir = output_dir / "upload"
                outputs = pipeline.split_autosomes(merged, upload_dir,
                                                   bgzip_path=pipeline.HTSLIB.bgzip_path,
                                                   chrom_prefix=panel_cfg["chrom_prefix"])
                print(f"Создано {len(outputs)} файлов в {upload_dir}")
                self.after(0, self._set_subprogress, 6, 1.0, "22 файла готовы")

                self.after(0, self._log_success, "=" * 60)
                self.after(0, self._log_success, "✅ Этапы 1-6 завершены!")
                self.after(0, self._log_success, f"📂 Загрузите 22 файла из: {upload_dir}")
                self.after(0, self._log_success, f"📁 Папка запуска: {output_dir}")
                self.after(0, self._log_success, "=" * 60)

                self.after(0, self._enable_mis_btn)
                self.after(0, self._refresh_run_history)
                # Файлы готовы — пользователю дальше на сайт MIS.
                self.after(0, self._set_wizard_step, 2)
                self.after(0, self._notify_done, True)

                if os.name == "nt":
                    os.startfile(str(upload_dir))
                # Промт "поправить ссылку для TopMed": Michigan Imputation
                # Server не предоставляет панель TOPMed r3 — открываем URL
                # именно ТОЙ панели, которая использовалась в этом запуске
                # (panel_cfg уже вычислен выше по коду метода), а не
                # захардкоженный Michigan для обеих панелей.
                webbrowser.open(panel_cfg["mis_upload_url"])

        except UserCancelledRun as e:
            # Промт "Доноры для VCF-источника: понятная отмена + общий
            # кэш сырых хромосом", Шаг 3: это ШТАТНОЕ, ожидаемое
            # прерывание запуска по решению пользователя (отказ качать
            # доноров, отмена при несоответствии источника и т.п.) —
            # никакой Python-traceback здесь не несёт диагностической
            # ценности и только пугает, выглядя как программный сбой.
            full_msg = str(e)
            short_msg = full_msg if len(full_msg) <= 120 else full_msg[:117] + "..."
            self.after(0, self._log_error, f"⚠ Запуск прерван: {full_msg}")
            self.after(0, lambda m=short_msg: self.stage_lbl.configure(text=f"⚠ Прервано: {m}"))
            self.after(0, self._notify_done, False)
        except Exception as e:
            # Формируем текст СРАЗУ (а не внутри lambda) — переменная
            # исключения 'e' удаляется Python'ом при выходе из блока except,
            # так что ленивое обращение к ней в lambda из self.after
            # приводило к "NameError: cannot access free variable 'e'"
            # вместо показа настоящей ошибки.
            full_msg = str(e)
            short_msg = full_msg if len(full_msg) <= 120 else full_msg[:117] + "..."
            tb_text = traceback.format_exc()
            # Печатаем полный traceback и в консоль (мы уже вне блока
            # redirect_stdout/redirect_stderr, поэтому LogRedirector его не
            # перехватит) — это упрощает диагностику будущих ошибок.
            print(tb_text, file=sys.stderr)
            self.after(0, self._log_error, f"❌ Ошибка: {full_msg}")
            self.after(0, self._log_error, tb_text)
            self.after(0, lambda m=short_msg: self.stage_lbl.configure(text=f"❌ Ошибка: {m}"))
            self.after(0, self._notify_done, False)
        finally:
            stdout_redirector.close()
            stderr_redirector.close()
            self._cancel_donor_download.clear()
            self.after(0, self._hide_cancel_donor_btn)
            self.after(0, self._enable_start_btn)

    # -----------------------------------------------------------------------
    # Этап 7 (Задача 7: прогресс 0-50%/50-100% + краткие сообщения)
    # -----------------------------------------------------------------------
    def _run_stage_7(self):
        rsq_threshold = self._get_rsq_threshold()
        fmt = self._get_format_key()
        output_path = None
        # Промт "Именованные папки запуска": та же папка запуска, что
        # использовалась на Этапе 1-6 (self.current_run_dir) — выставляется
        # либо в _on_start()/_run_stages_1_6(), либо в _on_continue_run()
        # (после выбора завершённого запуска из истории, в т.ч. после
        # перезапуска GUI). _on_mis() уже проверил, что она не None.
        run_log_path = (self.current_run_dir / "run.log") if self.current_run_dir else None
        stdout_redirector = LogRedirector(self.log_q, run_log_path)
        stderr_redirector = LogRedirector(self.log_q, run_log_path)
        try:
            with contextlib.redirect_stdout(stdout_redirector), \
                 contextlib.redirect_stderr(stderr_redirector):

                if self.current_run_dir is None:
                    raise RuntimeError(
                        "Внутренняя ошибка: папка запуска не была выбрана "
                        "перед стартом Этапа 7."
                    )

                curl = self.curl_tf.get("1.0", "end").strip()
                pwd = self.pwd_tf.get().strip()

                output_dir = self.current_run_dir
                # Промт "встроить лифтовер HRC/TopMed в gui/app.py": панель
                # ЭТОГО конкретного запуска читается из run_info.json
                # (записан в _run_stages_1_6() через pipeline.save_run_info()),
                # а НЕ из текущего состояния self.panel_dd — пользователь мог
                # сменить выбор в выпадающем списке между Этапом 1-6 и
                # Этапом 7 в рамках одной сессии, либо продолжить Этап 7 для
                # запуска из истории (_on_continue_run()), где self.panel_dd
                # вообще не выставлялся под этот конкретный запуск.
                run_info = pipeline.load_run_info(output_dir)
                panel = run_info.get("panel", pipeline.DEFAULT_PANEL)

                results_dir = output_dir / "rerun_results"
                results_dir.mkdir(exist_ok=True)

                self.after(0, self._set_stage7_progress, 0.0, "Скачивание результатов MIS...")
                # Промт "не видно прогресса на Шаге 3": следим за папкой
                # rerun_results теми же средствами, что и за донорами —
                # download_mis_results_smart() своего прогресса не отдаёт,
                # а размеры файлов на диске честны независимо от того,
                # чем именно качается и распаковывается архив.
                self.after(0, self._start_file_watch, [results_dir], "mis")
                # Промт "проверять уже скачанные файлы + предлагать повтор
                # при ошибке": уже присутствующие в results_dir непустые
                # файлы пропускаются автоматически (см.
                # MISAdapter.download_results()), а при сбое скачивания
                # конкретного файла показывается диалог с предложением
                # повторить именно его, не прерывая скачивание остальных.
                pipeline.download_mis_results_smart(
                    curl, results_dir, pwd,
                    on_file_error=self._prompt_file_download_retry,
                )
                self.after(0, self._stop_file_watch)
                self.after(0, self.mis_files_box.pack_forget)
                self.after(0, self._set_stage7_progress, 0.5, "Скачивание завершено, начинаю сборку...")

                tmpl_path = Path(self.tmpl_tf.get())
                skeleton = pipeline.extract_skeleton(tmpl_path, autosomes_only=False)
                # panel_pos — ВРЕМЕННАЯ копия координат скелета, используемая
                # только как позиционный фильтр (-R) для load_imputed_genotypes().
                # Сам skeleton (список SkeletonRow) ниже, в assemble_final(),
                # передаётся БЕЗ изменений и остаётся в GRCh37, как и раньше —
                # template/skeleton.py и template/assembler.py по замыслу
                # промта не должны знать о существовании лифтовера вообще.
                panel_pos = [(r.chrom, r.pos) for r in skeleton]

                # --- Лифтовер результата MIS (промт "встроить лифтовер
                # HRC/TopMed в gui/app.py", Этап 7) -------------------------
                # Результат Michigan Imputation Server для panel="topmed"
                # приходит в координатах GRCh38 (панель загружалась на MIS
                # именно в этой сборке — см. build_vcf(chrom_prefix=...) на
                # Этапе 2), тогда как skeleton/panel_pos всегда в GRCh37.
                # Без форвард-лифтовера panel_pos фильтр -R в
                # load_imputed_genotypes() почти ничего не найдёт (координаты
                # одного и того же локуса в разных сборках, как правило,
                # разные числа) и МОЛЧА вернёт пустой/почти пустой результат.
                forward_liftover = None
                reverse_liftover = None
                if panel != pipeline.DEFAULT_PANEL:
                    self.after(0, self._set_stage7_progress, 0.55, "Подготовка лифтовера результата...")
                    forward_liftover = pipeline._build_liftover(panel, direction="forward")
                    reverse_liftover = pipeline._build_liftover(panel, direction="reverse")
                    if forward_liftover is not None:
                        panel_pos = pipeline.liftover_positions_forward(panel_pos, forward_liftover)
                    else:
                        print(
                            f"⚠ Не удалось построить форвард-лифтовер для "
                            f"панели '{panel}' — фильтрация результата MIS "
                            f"будет выполнена по нелифтованным (GRCh37) "
                            f"координатам, что для этой панели, вероятно, "
                            f"приведёт к почти пустому результату."
                        )
                    if reverse_liftover is None:
                        print(
                            f"⚠ Не удалось построить реверс-лифтовер для "
                            f"панели '{panel}' — импутированные генотипы "
                            f"останутся в координатах панели и, скорее всего, "
                            f"не совпадут со скелетом трафарета (GRCh37)."
                        )

                self.after(0, self._set_stage7_progress, 0.6, "Загрузка импутированных генотипов...")
                imputed = pipeline.load_imputed_genotypes(
                    results_dir, "genotek", panel_pos,
                    rsq_threshold=rsq_threshold,
                    bcftools_path=pipeline.HTSLIB.bcftools_path,
                    tabix_path=pipeline.HTSLIB.tabix_path,
                )

                if reverse_liftover is not None:
                    self.after(0, self._set_stage7_progress, 0.65, "Обратный лифтовер результата в GRCh37...")
                    imputed, dropped = pipeline.liftback_imputed_genotypes(imputed, reverse_liftover)
                    print(
                        f"ℹ Обратный лифтовер: {len(imputed)} позиций "
                        f"перенесено в GRCh37, {dropped} отброшено"
                    )

                with (output_dir / "parse_result.pkl").open("rb") as f:
                    result = pickle.load(f)
                measured = pipeline.load_measured_genotypes(result.variants)

                # ⚠ Находка (промт "встроить лифтовер HRC/TopMed в
                # gui/app.py", проверка перед сдачей): result.variants для
                # panel="topmed" УЖЕ в координатах GRCh38 — форвард-лифтовер
                # применяется ВНУТРИ parser_fn на Этапе 1 (см. docstring
                # adapters/base.py::ParsedVariant, adapters/ftdna_v3.py/
                # adapters/myheritage_v5.py::parse_*(liftover=...)), ДО того
                # как result сохраняется в parse_result.pkl. Раньше здесь
                # лифтовался обратно только `imputed` (результат MIS), а
                # `measured` (напрямую измеренные чипом позиции) — нет.
                # Без этого исправления genotypes.get(f"{chrom}_{pos}") в
                # assemble_final() (координаты скелета — всегда GRCh37) не
                # находил бы НИ ОДНОЙ измеренной позиции для panel="topmed" —
                # все они молча уходили бы в "--", как будто человек вообще
                # не тестировался на этих позициях чипа. Тот же
                # reverse_liftover, что уже используется для imputed выше.
                if reverse_liftover is not None:
                    measured, measured_dropped = pipeline.liftback_imputed_genotypes(
                        measured, reverse_liftover,
                    )
                    print(
                        f"ℹ Обратный лифтовер измеренных позиций: "
                        f"{len(measured)} перенесено в GRCh37, "
                        f"{measured_dropped} отброшено"
                    )

                genotypes = pipeline.merge_dictionaries(imputed, measured)

                self.after(0, self._set_stage7_progress, 0.85, "Сборка финального файла...")
                # Промт "итоговый файл в отдельной папке": результат больше
                # не ложится в рабочую папку запуска вперемешку с
                # промежуточными VCF и логами — все итоговые файлы всех
                # запусков собираются в results/ рядом с программой, с
                # именем запуска в начале, чтобы не путались между собой.
                # _unique_result_path() не даёт повторной сборке того же
                # запуска (например с другим Rsq) молча затереть прошлый файл.
                results_root = _results_dir()
                output_path = _unique_result_path(
                    results_root / f"{self.current_run_name}_genotek_23andme_{fmt}.txt"
                )
                pipeline.assemble_final(
                    skeleton, genotypes, output_path,
                    format_version=fmt,
                    template_path=tmpl_path,
                )

                self.after(0, self._set_stage7_progress, 0.95, "Проверка результата...")
                validation = pipeline.validate_output(output_path, tmpl_path, fmt)
                self.after(0, self._set_stage7_progress, 1.0, "Готово")

                if validation.is_valid:
                    short_msg = f"✅ Готово! {output_path.name} (call rate: {validation.call_rate:.2f}%)"
                    pipeline.save_run_info(
                        output_dir,
                        call_rate=round(validation.call_rate, 2),
                        finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        format=fmt,
                        rsq_threshold=rsq_threshold,
                        result_file=str(output_path),
                    )
                    self._last_result_path = output_path
                    self.after(0, self._refresh_result_label)
                    self.after(0, self._refresh_run_history)
                    self.after(0, self._log_success, "=" * 60)
                    self.after(0, self._log_success, f"✅ ГОТОВО! Файл: {output_path}")
                    self.after(0, self._log_success, f"   Call rate: {validation.call_rate:.2f}%")
                    self.after(0, self._log_success, f"   Формат: {fmt}")
                    self.after(0, self._log_success, f"   Порог Rsq: {rsq_threshold:.2f}")
                    if validation.template_duplicate_positions:
                        # Не ошибка сборки (см. докстринг
                        # ValidationResult.template_duplicate_positions в
                        # assembler.py) — позиции унаследованы из самого
                        # трафарета (template.txt), просто информируем.
                        self.after(
                            0, self._log_success,
                            f"   ℹ В трафарете есть {len(validation.template_duplicate_positions)} "
                            f"дублирующихся позиций (напр. {validation.template_duplicate_positions[0]}) "
                            f"— это особенность template.txt, не ошибка сборки.",
                        )
                    self.after(0, self._log_success, "=" * 60)
                    self.after(0, lambda: self.stage_lbl.configure(text=short_msg))
                    self.after(0, self._notify_done, True)
                    # Открываем именно папку с результатом, а не рабочую
                    # папку запуска: пользователю нужен собранный файл.
                    self.after(0, self._open_in_file_manager, results_root)
                else:
                    first_err = validation.errors[0] if validation.errors else "неизвестная ошибка валидации"
                    short_msg = f"❌ Ошибка: {first_err}"
                    for err in validation.errors:
                        self.after(0, self._log_error, f"❌ {err}")
                    self.after(0, lambda m=short_msg: self.stage_lbl.configure(text=m))
                    self.after(0, self._notify_done, False)

        except Exception as e:
            full = str(e)
            short = full if len(full) <= 120 else full[:117] + "..."
            self.after(0, self._log_error, f"❌ Ошибка: {full}")
            self.after(0, lambda m=short: self.stage_lbl.configure(text=f"❌ {m}"))
            self.after(0, self._notify_done, False)
        finally:
            stdout_redirector.close()
            stderr_redirector.close()
            self.after(0, self._stop_file_watch)
            self.after(0, self._enable_mis_btn)

    # -----------------------------------------------------------------------
    # Thread-safe обновления UI
    # -----------------------------------------------------------------------
    def _log_success(self, msg: str):
        self.log_text.insert("end", msg + "\n", "success")
        self.log_text.see("end")

    def _log_error(self, msg: str):
        self.log_text.insert("end", msg + "\n", "error")
        self.log_text.see("end")

    def _enable_start_btn(self):
        self.start_btn.configure(state="normal")
        self.running = False
        self._run_started_at = None
        self._refresh_mis_btn_state()

    def _enable_mis_btn(self):
        """
        Раньше просто включал кнопку Шага 3. Теперь решение о её
        доступности принимает _refresh_mis_btn_state(): одной готовности
        файлов мало — нужны ещё curl-команда и пароль из письма, и когда
        их нет, пользователь видит под кнопкой, чего именно не хватает.
        """
        self.running = False
        self._run_started_at = None
        self._refresh_mis_btn_state()


# ---------------------------------------------------------------------------
# Точка входа
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = App()

    app.log_text.tag_config("success", foreground="#4CAF50")
    app.log_text.tag_config("error", foreground="#F44336")
    app.log_text.tag_config("stage", foreground="#2196F3")
    app.log_text.tag_config("progress", foreground="#64B5F6")

    app.mainloop()