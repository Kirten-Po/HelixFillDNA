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
#   - template/  -> только __init__.py (шаблоны *.txt пользователь выбирает
#     сам через диалог "Обзор" — не включаем 38 МБ template_v3/v5.txt).
#     Если хотите дать пример трафарета — раскомментируйте отдельной строкой.
# ---------------------------------------------------------------------------
app_datas = [
    (os.path.join(PROJECT_ROOT, "bin"), "bin"),
    # Пример шаблона (раскомментировать при необходимости):
    # (os.path.join(PROJECT_ROOT, "template", "template_v3.txt"), "template"),
    # (os.path.join(PROJECT_ROOT, "template", "template_v5.txt"), "template"),
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
    # ИСПРАВЛЕНИЕ БАГА "bcftools.exe не найден":
    # начиная с PyInstaller 6.0 по умолчанию contents_directory="_internal" —
    # все данные и бинарники (в т.ч. bin/bcftools.exe) кладутся в
    # dist\HelixFillDNA\_internal\, а не плоско рядом с HelixFillDNA.exe.
    # main.py::HtslibTools._find() и gui/app.py вычисляют bin_dir как папку
    # РЯДОМ с exe (pre-6.0 layout), поэтому bcftools.exe не находился —
    # bgzip/tabix "работали" только если случайно были в системном PATH.
    # contents_directory="." возвращает старый плоский макет: bin/ снова
    # лежит прямо в dist\HelixFillDNA\bin\, без правки Python-кода.
    contents_directory=".",
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
