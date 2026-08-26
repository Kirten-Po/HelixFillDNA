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
import logging
import os
import pickle
import queue
import shutil
import subprocess
import sys
import threading
import tkinter as tk
import traceback
import webbrowser
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox

import customtkinter as ctk

# === Корень проекта ===
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

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


# === Импорты пайплайна ===
import main as pipeline
import download_donors
from adapters.ftdna_v3 import ReferenceGenome
from core import archive_utils
from core import network_utils


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


def attach_hotkeys(widget):
    """
    Привязывает Ctrl+C/V/X/A к виджету.

    Задача 1: каждый обработчик возвращает "break". Без этого Tk/CustomTkinter
    сначала выполняет наш обработчик (который сам вставляет/копирует текст),
    а затем ЕЩЁ РАЗ прогоняет встроенный биндинг того же события — в
    результате текст вставлялся дважды. "break" останавливает дальнейшую
    обработку события этим виджетом, поэтому срабатывает только наш код.
    """
    def _bind(seq, action_fn):
        def handler(event, w=widget, fn=action_fn):
            fn(w)
            return "break"
        widget.bind(seq, handler, add="+")

    _bind("<Control-c>", _do_copy)
    _bind("<Control-C>", _do_copy)
    _bind("<Control-v>", _do_paste)
    _bind("<Control-V>", _do_paste)
    _bind("<Control-x>", _do_cut)
    _bind("<Control-X>", _do_cut)
    _bind("<Control-a>", _do_select_all)
    _bind("<Control-A>", _do_select_all)


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
        super().__init__()

        self.title("HelixFillDNA")
        self.geometry("1000x750")
        self.minsize(850, 650)

        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

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
    def _build_settings_tab(self):
        scroll = ctk.CTkScrollableFrame(self.tab_settings)
        scroll.pack(fill="both", expand=True, padx=10, pady=10)

        ctk.CTkLabel(
            scroll, text="Основные настройки",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="Источник данных:").pack(anchor="w")
        source_names = [v["name"] for v in pipeline.SOURCES.values()]
        self.source_dd = ctk.CTkOptionMenu(scroll, values=source_names, width=400)
        self.source_dd.set(source_names[0])
        self.source_dd.pack(anchor="w", pady=(0, 15))

        # --- Шаг 1 промта "HRC / TopMed": выбор референсной панели -------
        ctk.CTkLabel(scroll, text="Референсная панель импутации:").pack(anchor="w")
        panel_names = [v["display_name"] for v in pipeline.REFERENCE_PANELS.values()]
        self.panel_dd = ctk.CTkOptionMenu(
            scroll, values=panel_names, width=400,
            command=lambda _choice: self._on_panel_changed(),
        )
        self.panel_dd.set(pipeline.REFERENCE_PANELS[pipeline.DEFAULT_PANEL]["display_name"])
        self.panel_dd.pack(anchor="w", pady=(0, 5))

        self.panel_warning_lbl = ctk.CTkLabel(
            scroll, text="", justify="left", text_color="#F9A825", wraplength=700,
        )
        self.panel_warning_lbl.pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(scroll, text="Файл с данными:").pack(anchor="w")
        row1 = ctk.CTkFrame(scroll, fg_color="transparent")
        row1.pack(fill="x", pady=(0, 15))
        self.input_tf = ctk.CTkEntry(row1, placeholder_text="Выберите файл...", width=600)
        self.input_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(row1, text="Обзор", width=100,
                      command=lambda: self._pick_file(self.input_tf)).pack(side="right")

        ctk.CTkLabel(
            scroll,
            text=("ℹ Референсный геном для FTDNA/MyHeritage проверяется\n"
                  "автоматически при запуске под выбранную выше панель. Если\n"
                  "файла нет в reference/<панель>/ — он будет скачан и распакован\n"
                  "(размер зависит от сборки генома выбранной панели)."),
            justify="left", text_color="gray60",
        ).pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(scroll, text="Трафарет (template_v3.txt):").pack(anchor="w")
        row3 = ctk.CTkFrame(scroll, fg_color="transparent")
        row3.pack(fill="x", pady=(0, 15))
        self.tmpl_tf = ctk.CTkEntry(row3, placeholder_text="Выберите файл...", width=600)
        self.tmpl_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(row3, text="Обзор", width=100,
                      command=lambda: self._pick_file(self.tmpl_tf)).pack(side="right")

        ctk.CTkLabel(scroll, text="Папка с бинарниками htslib:").pack(anchor="w")
        row_bin = ctk.CTkFrame(scroll, fg_color="transparent")
        row_bin.pack(fill="x", pady=(0, 5))
        self.bin_tf = ctk.CTkEntry(row_bin, width=600)
        self.bin_tf.insert(0, str(_detect_bin_dir()))
        self.bin_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
        ctk.CTkButton(
            row_bin, text="🌐 Диагностика сети", width=170,
            command=self._on_diagnose_network,
        ).pack(side="right")
        ctk.CTkLabel(
            scroll,
            text=("ℹ Проверяет CA-сертификаты для libcurl (bcftools) и наличие "
                  "конфликтующего curl.exe в папке бинарников — нужно только "
                  "для ускоренного удалённого скачивания доноров (Этап 3), "
                  "само скачивание работает и без этого через обычный "
                  "полный путь."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 15))

        ctk.CTkFrame(scroll, height=2, fg_color="gray40").pack(fill="x", pady=15)

        ctk.CTkLabel(
            scroll, text="Параметры вывода",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="Формат вывода:").pack(anchor="w")
        self.format_dd = ctk.CTkOptionMenu(
            scroll,
            values=["v3 (LF, ~97% call rate)", "v5 (CRLF, ~92% call rate)"],
            width=400,
        )
        self.format_dd.set("v3 (LF, ~97% call rate)")
        self.format_dd.pack(anchor="w", pady=(0, 15))

        ctk.CTkLabel(
            scroll, text="Порог Rsq (качество импутации, от 0 до 1):",
            font=ctk.CTkFont(weight="bold"),
        ).pack(anchor="w")
        rsq_info = (
            "Чем выше Rsq — тем надёжнее генотип, но тем меньше позиций проходит фильтр.\n"
            "0.30 — стандартный порог MIS: максимум позиций, часть — с невысоким качеством.\n"
            "0.80 — баланс количества и качества.\n"
            "0.90 — только высокоточные варианты (меньше позиций, но надёжнее).\n"
            "0.95+ — максимальное качество для критичных задач."
        )
        ctk.CTkLabel(scroll, text=rsq_info, justify="left", text_color="gray60").pack(
            anchor="w", pady=(4, 8)
        )

        row_rsq = ctk.CTkFrame(scroll, fg_color="transparent")
        row_rsq.pack(fill="x", pady=(0, 5))
        ctk.CTkLabel(row_rsq, text="Значение:").pack(side="left", padx=(0, 10))
        self.rsq_entry = ctk.CTkEntry(row_rsq, width=100, placeholder_text="0.30")
        self.rsq_entry.insert(0, "0.30")
        self.rsq_entry.pack(side="left")
        self.rsq_entry.bind("<KeyRelease>", lambda e: self._validate_rsq_entry())

        self.rsq_status_lbl = ctk.CTkLabel(
            scroll, text="✓ Порог принят: 0.30", text_color="#4CAF50",
        )
        self.rsq_status_lbl.pack(anchor="w", pady=(0, 15))

        # Задача C: опциональная нормализация multiallelic-сайтов
        # (bcftools norm -m-both). НЕ входит в критический путь фикса
        # Invalid alleles (это делают Задачи A/B) — отдельная оптимизация,
        # поэтому выключена по умолчанию, как и в CLI (--normalize).
        self.normalize_var = ctk.BooleanVar(value=False)
        self.normalize_cb = ctk.CTkCheckBox(
            scroll,
            text="Нормализовать multiallelic-сайты перед split (bcftools norm -m-both, опционально)",
            variable=self.normalize_var,
        )
        self.normalize_cb.pack(anchor="w", pady=(0, 15))

        # Задача D: опциональное переиспользование доноров между разными
        # людьми на одном чипе (широкая сигнатура вместо строгой).
        self.reuse_donors_var = ctk.BooleanVar(value=False)
        self.reuse_donors_cb = ctk.CTkCheckBox(
            scroll,
            text=("Переиспользовать доноров между разными людьми на одном чипе "
                  "(экспериментально, Задача D)"),
            variable=self.reuse_donors_var,
        )
        self.reuse_donors_cb.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            scroll,
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
        # итоговое покрытие. По умолчанию теперь используется ВСЯ доступная
        # EUR-подвыборка панели (обычно порядка 500 человек) — это честное
        # увеличение размера случайной выборки, а не подбор конкретных
        # образцов под конкретные позиции чипа (такой алгоритмический
        # "умный" подбор сознательно не реализован — см. пояснения в чате).
        self.eur_all_var = ctk.BooleanVar(value=True)
        self.eur_all_cb = ctk.CTkCheckBox(
            scroll,
            text=("Использовать всех доступных EUR-доноров 1000 Genomes "
                  "(уменьшает Monomorphic sites на QC MIS, но увеличивает "
                  "трафик/время скачивания доноров)"),
            variable=self.eur_all_var,
            command=self._on_eur_all_toggled,
        )
        self.eur_all_cb.pack(anchor="w", pady=(0, 5))

        row_eur_count = ctk.CTkFrame(scroll, fg_color="transparent")
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
            scroll, text="✓ Будут использованы все доступные EUR-доноры (~500)",
            text_color="#4CAF50",
        )
        self.eur_count_status_lbl.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            scroll,
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
        # ценой постоянного места на диске. По умолчанию выключен.
        self.raw_cache_var = ctk.BooleanVar(value=False)
        self.raw_cache_cb = ctk.CTkCheckBox(
            scroll,
            text=("Хранить сырые (нефильтрованные) хромосомы 1000 Genomes для "
                  "повторного использования между разными источниками/чипами "
                  "(~десятки ГБ на диске, экономит трафик при последующих запусках)"),
            variable=self.raw_cache_var,
        )
        self.raw_cache_cb.pack(anchor="w", pady=(0, 5))
        ctk.CTkLabel(
            scroll,
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

    # -----------------------------------------------------------------------
    # Вкладка "Запуск"
    # -----------------------------------------------------------------------
    def _build_run_tab(self):
        scroll = ctk.CTkScrollableFrame(self.tab_run)
        scroll.pack(fill="both", expand=True, padx=10, pady=10)

        # --- Промт "Именованные папки запуска": имя/история запусков ------
        ctk.CTkLabel(
            scroll, text="Запуск (папка результатов)",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="Название запуска:").pack(anchor="w")
        run_name_row = ctk.CTkFrame(scroll, fg_color="transparent")
        run_name_row.pack(fill="x", pady=(0, 5))
        self.run_name_tf = ctk.CTkEntry(run_name_row, width=200)
        self.run_name_tf.pack(side="left", padx=(0, 10))
        ctk.CTkButton(
            run_name_row, text="Обновить историю", width=170,
            command=self._refresh_run_history,
        ).pack(side="left")
        attach_input_features(self.run_name_tf)

        ctk.CTkLabel(
            scroll,
            text=("ℹ Каждый запуск пишет файлы в свою папку "
                  "output/runs/<название>/ — донор-кэш (donors/) общий "
                  "для всех запусков и не дублируется. По умолчанию "
                  "название — следующий свободный номер; можно ввести "
                  "своё (например имя человека)."),
            justify="left", text_color="gray60", wraplength=700,
        ).pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="История запусков:").pack(anchor="w")
        history_row = ctk.CTkFrame(scroll, fg_color="transparent")
        history_row.pack(fill="x", pady=(0, 5))
        self.run_history_dd = ctk.CTkOptionMenu(history_row, values=["(нет запусков)"], width=420)
        self.run_history_dd.pack(side="left", padx=(0, 10))
        ctk.CTkButton(
            history_row, text="Продолжить (Этап 2)", width=180,
            command=self._on_continue_run,
        ).pack(side="left", padx=(0, 5))
        ctk.CTkButton(
            history_row, text="Переименовать", width=140,
            command=self._on_rename_run,
        ).pack(side="left", padx=(0, 5))
        ctk.CTkButton(
            history_row, text="Открыть папку", width=140,
            command=self._on_open_run_folder,
        ).pack(side="left")

        self.active_run_lbl = ctk.CTkLabel(scroll, text="Активный запуск: нет", text_color="gray60")
        self.active_run_lbl.pack(anchor="w", pady=(0, 15))

        ctk.CTkFrame(scroll, height=2, fg_color="gray40").pack(fill="x", pady=(0, 15))

        ctk.CTkLabel(
            scroll, text="Этап 1: Подготовка файлов для MIS",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        self.stage_lbl = ctk.CTkLabel(scroll, text="Готов к запуску",
                                      font=ctk.CTkFont(size=16, weight="bold"))
        self.stage_lbl.pack(anchor="w", pady=(0, 5))

        self.progress = ctk.CTkProgressBar(scroll, height=15)
        self.progress.pack(fill="x", pady=(0, 15))
        self.progress.set(0)

        self.start_btn = ctk.CTkButton(
            scroll, text="Запустить этапы 1-6 (до MIS)",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color="#2E7D32", hover_color="#1B5E20",
            height=50, corner_radius=12,
            command=self._on_start,
        )
        self.start_btn.pack(fill="x", pady=(0, 10))

        # Задача 1: активна только пока идёт скачивание доноров (Этап 3).
        # Отмена реально прерывает текущий subprocess (curl) внутри
        # download_donors.py, а не ждёт окончания текущей хромосомы.
        self.cancel_donor_btn = ctk.CTkButton(
            scroll, text="Отменить скачивание доноров", width=280,
            fg_color="#B71C1C", hover_color="#7F0000",
            state="disabled",
            command=self._on_cancel_donor_download,
        )
        self.cancel_donor_btn.pack(fill="x", pady=(0, 20))

        ctk.CTkFrame(scroll, height=2, fg_color="gray40").pack(fill="x", pady=15)

        ctk.CTkLabel(
            scroll, text="Этап 2: Imputation на Michigan Server",
            font=ctk.CTkFont(size=20, weight="bold"),
        ).pack(anchor="w", pady=(0, 10))

        # Шаг 1 промта TopMed: текст инструкции про Reference Panel/Array
        # Build теперь формируется динамически из выбранной на вкладке
        # "Подготовка" панели (self.panel_dd), а не захардкожен под HRC —
        # обновляется через _refresh_run_instructions(), вызываемую и при
        # построении вкладки, и при смене self.panel_dd.
        self.run_instructions_lbl = ctk.CTkLabel(scroll, text="", justify="left")
        self.run_instructions_lbl.pack(anchor="w", pady=(0, 10))

        ctk.CTkLabel(scroll, text="curl-команда из письма MIS:").pack(anchor="w")
        self.curl_tf = ctk.CTkTextbox(scroll, height=80, width=700)
        self.curl_tf.pack(fill="x", pady=(0, 10))

        ctk.CTkLabel(scroll, text="Пароль из письма:").pack(anchor="w")
        pwd_row = ctk.CTkFrame(scroll, fg_color="transparent")
        pwd_row.pack(fill="x", pady=(0, 15))
        self.pwd_tf = ctk.CTkEntry(pwd_row, show="*", width=520, placeholder_text="Пароль")
        self.pwd_tf.pack(side="left", fill="x", expand=True, padx=(0, 10))
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
            scroll, text="Скачать результаты и собрать финальный файл",
            font=ctk.CTkFont(size=16, weight="bold"),
            fg_color="#1565C0", hover_color="#0D47A1",
            height=50, corner_radius=12,
            state="disabled",
            command=self._on_mis,
        )
        self.mis_btn.pack(fill="x", pady=(0, 10))

        attach_input_features(self.curl_tf)
        attach_input_features(self.pwd_tf)

        self._refresh_run_instructions()
        self._refresh_run_name_suggestion()
        self._refresh_run_history()

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
        # Метод может вызываться до построения вкладки "Запуск" (первичная
        # синхронизация в конце _build_settings_tab) — тогда просто пропускаем.
        if hasattr(self, "run_instructions_lbl"):
            self._refresh_run_instructions()

    def _refresh_run_instructions(self):
        panel = self._get_panel_key()
        cfg = pipeline.REFERENCE_PANELS[panel]
        text = (
            "1. После завершения этапов 1-6 откроется папка с 22 файлами.\n"
            f"2. Загрузите их на {cfg['mis_upload_url']}\n"
            f"3. Reference Panel: {cfg['mis_panel_value']}, Population: EUR\n"
            "4. Дождитесь письма со ссылкой и паролем, вставьте их ниже.\n"
            "💡 Совет: Ctrl+V для вставки, правая кнопка мыши — контекстное меню"
        )
        self.run_instructions_lbl.configure(text=text)

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
        Кнопка «▶ Продолжить (Этап 2)» — делает выбранный из истории
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
        self.mis_btn.configure(state="normal")
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
        """Кнопка «📂 Открыть папку» — открывает папку выбранного из
        истории запуска в системном файловом менеджере."""
        run_dir = self._selected_history_run()
        if run_dir is None:
            messagebox.showwarning("Предупреждение", "Выберите запуск из списка истории")
            return
        if os.name == "nt":
            os.startfile(str(run_dir))
        else:
            try:
                subprocess.run(["xdg-open", str(run_dir)], check=False)
            except Exception:
                messagebox.showinfo("Папка запуска", str(run_dir))

    def _poll_logs(self):
        try:
            while True:
                msg = self.log_q.get_nowait()
                if any(msg.startswith(p) for p in ("✓", "✅")):
                    self.log_text.insert("end", msg + "\n", "success")
                elif any(msg.startswith(p) for p in ("✗", "❌", "ОШИБКА")):
                    self.log_text.insert("end", msg + "\n", "error")
                elif msg.startswith("["):
                    self.log_text.insert("end", msg + "\n", "stage")
                else:
                    self.log_text.insert("end", msg + "\n")
                self.log_text.see("end")
        except queue.Empty:
            pass
        self.after(100, self._poll_logs)

    def _validate_settings(self) -> str | None:
        if not self.input_tf.get():
            return "Выберите файл с данными"
        if not self.tmpl_tf.get():
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
        self.cancel_donor_btn.configure(state="normal")

    def _hide_cancel_donor_btn(self):
        self.cancel_donor_btn.configure(state="disabled")

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
        self.cancel_donor_btn.configure(state="disabled")
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
        Плавный прогресс внутри этапа: общий прогресс = (stage_n - 1 + sub_progress) / 7.
        sub_progress — доля выполнения текущего этапа (0.0 .. 1.0).
        Вызывается только через self.after(0, ...) из фонового потока.
        """
        sub_progress = max(0.0, min(1.0, sub_progress))
        overall = (stage_n - 1 + sub_progress) / 7
        self.progress.set(overall)
        self.stage_lbl.configure(text=f"[{stage_n}/7] {text}")
        print(f"[{stage_n}/7] {text}")

    # --- Задача 7: прогресс + краткие сообщения этапа 7 -------------------
    def _set_stage7_progress(self, frac: float, text: str):
        """frac в диапазоне 0.0..1.0: 0-0.5 скачивание, 0.5-1.0 сборка."""
        frac = max(0.0, min(1.0, frac))
        self.progress.set(frac)
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
        self.start_btn.configure(state="disabled")
        self.mis_btn.configure(state="disabled")
        self.progress.set(0)
        self.stage_lbl.configure(text="Запуск...")
        threading.Thread(target=self._run_stages_1_6, daemon=True).start()

    def _on_mis(self):
        if self.running:
            return
        if self.current_run_dir is None:
            messagebox.showwarning(
                "Предупреждение",
                "Не выбран активный запуск. Выполните Этап 1-6, либо "
                "выберите завершённый запуск в «Истории запусков» и "
                "нажмите «▶ Продолжить (Этап 2)».",
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
        self.mis_btn.configure(state="disabled")
        self.progress.set(0)
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
                if pipeline._supports_liftover(source):
                    result = parser_fn(csv_path, reference, liftover=liftover)
                else:
                    result = parser_fn(csv_path, reference)
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
                output_path = output_dir / f"genotek_23andme_{fmt}.txt"
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
                    )
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
                    if os.name == "nt":
                        os.startfile(str(output_dir))
                else:
                    first_err = validation.errors[0] if validation.errors else "неизвестная ошибка валидации"
                    short_msg = f"❌ Ошибка: {first_err}"
                    for err in validation.errors:
                        self.after(0, self._log_error, f"❌ {err}")
                    self.after(0, lambda m=short_msg: self.stage_lbl.configure(text=m))

        except Exception as e:
            full = str(e)
            short = full if len(full) <= 120 else full[:117] + "..."
            self.after(0, self._log_error, f"❌ Ошибка: {full}")
            self.after(0, lambda m=short: self.stage_lbl.configure(text=f"❌ {m}"))
        finally:
            stdout_redirector.close()
            stderr_redirector.close()
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

    def _enable_mis_btn(self):
        self.mis_btn.configure(state="normal")
        self.running = False


# ---------------------------------------------------------------------------
# Точка входа
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app = App()

    app.log_text.tag_config("success", foreground="#4CAF50")
    app.log_text.tag_config("error", foreground="#F44336")
    app.log_text.tag_config("stage", foreground="#2196F3")

    app.mainloop()