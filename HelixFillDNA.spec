# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec для HelixFillDNA.
Сборка:
    pyinstaller HelixFillDNA.spec

Точка входа: gui/app.py (GUI на customtkinter).
Если нужен ещё и CLI-вариант (main.py) как отдельный exe — см. закомментированный
второй Analysis/EXE блок внизу файла.
"""

from PyInstaller.utils.hooks import collect_all
import os

block_cipher = None

PROJECT_ROOT = os.path.abspath(".")

# ---------------------------------------------------------------------------
# customtkinter и pyfaidx часто не подхватываются автоматически PyInstaller'ом
# (свои ассеты — темы/json/шрифты у customtkinter, компилируемые part у
# pyfaidx) — collect_all тянет datas + binaries + hiddenimports разом.
# ---------------------------------------------------------------------------
ctk_datas, ctk_binaries, ctk_hidden = collect_all("customtkinter")
pyfaidx_datas, pyfaidx_binaries, pyfaidx_hidden = collect_all("pyfaidx")

# ---------------------------------------------------------------------------
# Данные приложения:
#   - bin/  -> все htslib-бинарники и DLL (bcftools.exe/tabix.exe/bgzip.exe +
#     их зависимости + cacert.pem) — код ищет их по --bin-dir/PROJECT_ROOT/"bin",
#     поэтому кладём в подпапку bin/ рядом с exe.
#   - template/  -> только __init__.py (это python-пакет сборщика, сами
#     трафареты *.txt лежат не здесь, а в samples/, см. ниже).
#   - samples/   -> трафареты template_v3.txt (FTDNA) и template_v5.txt
#     (MyHeritage) СЮДА НЕ ВКЛЮЧАЮТСЯ СПЕЦИАЛЬНО. Промт "обычные /
#     продвинутые настройки": в обычном режиме GUI подставляет трафарет из
#     папки samples/ сам, и папка эта должна оказаться ровно в {app}\samples
#     у установленного приложения. PyInstaller, в зависимости от версии,
#     кладёт datas либо плоско рядом с exe, либо в _internal/ (та же
#     неопределённость, из-за которой существует _detect_bin_dir() в
#     gui/app.py) — поэтому трафареты кладёт НЕ PyInstaller, а Inno Setup:
#     см. секцию [Files] в HelixFillDNA.iss, там путь к {app}\samples задан
#     явно и от раскладки PyInstaller не зависит.
# ---------------------------------------------------------------------------
app_datas = [
    (os.path.join(PROJECT_ROOT, "bin"), "bin"),
]

a = Analysis(
    ["gui/app.py"],
    pathex=[PROJECT_ROOT],
    binaries=ctk_binaries + pyfaidx_binaries,
    datas=app_datas + ctk_datas + pyfaidx_datas,
    hiddenimports=[
        *ctk_hidden,
        *pyfaidx_hidden,
        # Модули проекта, импортируемые динамически/косвенно (main.py
        # делает `import main as pipeline` из gui/app.py, но сами
        # subpackage-и adapters/core/template тоже стоит перечислить явно
        # на случай, если PyInstaller не увидит их через анализ импортов
        # main.py/download_donors.py целиком).
        "adapters",
        "adapters.base",
        "adapters.ftdna_v3",
        "adapters.myheritage_v5",
        "adapters.vcf_source",
        "core",
        "core.archive_utils",
        "core.liftover",
        "core.network_utils",
        "core.pure_python_core",
        "template",
        "template.assembler",
        "template.skeleton",
        "main",
        "version",
        "mis_adapter",
        "download_donors",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # tests/ в exe не нужны — PyInstaller и так не тянет их сам, но
        # явный excludes страхует, если что-то импортирует их косвенно.
        "tests",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="HelixFillDNA",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,       # GUI-приложение — без консольного окна
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="app_icon.ico",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="HelixFillDNA",
    # ИСПРАВЛЕНИЕ БАГА "ModuleNotFoundError: No module named 'customtkinter'"
    # (v1.0.2): раньше здесь стоял contents_directory="." — идея была в том,
    # чтобы вернуть старый плоский макет (bin/ рядом с exe, без _internal),
    # так как main.py::HtslibTools._find() и gui/app.py::_detect_bin_dir()
    # исторически ждали bin/ именно там.
    #
    # Но выяснилось (см. известный баг PyInstaller 6.x, issue #8075 в их
    # трекере, закрыт мейнтейнерами как "won't fix"): файлы, добавленные
    # через collect_all() — а именно так подключены customtkinter и
    # pyfaidx выше в Analysis() — ФИЗИЧЕСКИ всё равно попадают в _internal,
    # даже когда contents_directory="." просит плоский макет. При этом
    # бутлоадер, ДОВЕРЯЯ contents_directory=".", настраивает sys.path БЕЗ
    # папки _internal вообще (он не знает, что она нужна). В итоге
    # customtkinter физически лежал в dist/HelixFillDNA/_internal/, но
    # Python туда не заглядывал при импорте — отсюда и
    # "ModuleNotFoundError: No module named 'customtkinter'" у конечных
    # пользователей, никак не связанный с антивирусом/UPX (тот фикс в
    # v1.0.1 был не по адресу, хоть и не вредный).
    #
    # Возвращаемся к ДЕФОЛТНОМУ макету с _internal — PyInstaller сам
    # корректно прописывает эту папку в sys.path, поэтому импорт снова
    # работает. gui/app.py::_detect_bin_dir() уже и так проверяет ОБА
    # варианта расположения bin/ (плоский и _internal/bin) — переход
    # обратно на _internal ничего не ломает в поиске bcftools/tabix/bgzip.
)

# ---------------------------------------------------------------------------
# Опционально: отдельный CLI-экзешник для main.py (например, для батч-запуска
# без GUI). Раскомментируйте блок ниже, если нужен второй exe.
# ---------------------------------------------------------------------------
# a_cli = Analysis(
#     ["main.py"],
#     pathex=[PROJECT_ROOT],
#     binaries=[],
#     datas=app_datas,
#     hiddenimports=[
#         "adapters", "adapters.base", "adapters.ftdna_v3",
#         "adapters.myheritage_v5", "adapters.vcf_source",
#         "core", "core.archive_utils", "core.liftover",
#         "core.network_utils", "core.pure_python_core",
#         "template", "template.assembler", "template.skeleton",
#         "mis_adapter", "download_donors",
#     ],
#     hookspath=[],
#     cipher=block_cipher,
# )
# pyz_cli = PYZ(a_cli.pure, a_cli.zipped_data, cipher=block_cipher)
# exe_cli = EXE(
#     pyz_cli, a_cli.scripts, [],
#     exclude_binaries=True,
#     name="HelixFillDNA-CLI",
#     console=True,
# )
