; ============================================================================
; HelixFillDNA.iss — скрипт Inno Setup для сборки Setup.exe.
;
; Как использовать (без командной строки):
;   1. Соберите папку dist\HelixFillDNA\ через PyInstaller:
;        python -m PyInstaller HelixFillDNA.spec
;   2. Установите Inno Setup: https://jrsoftware.org/isdl.php
;   3. Откройте этот файл (HelixFillDNA.iss) двойным кликом — откроется
;      Inno Setup Compiler.
;   4. Нажмите "Compile" (или Ctrl+F9, или кнопка с зелёной шестернёй).
;   5. Готовый установщик появится в папке Output\ рядом с этим .iss —
;      HelixFillDNA-Setup-1.0.0.exe.
;
; Установщик НЕ требует от пользователя ни Python, ни pip, ни командной
; строки — PyInstaller уже упаковал customtkinter/pyfaidx/весь рантайм
; внутрь dist\HelixFillDNA\, Inno Setup только копирует эти файлы в
; Program Files и создаёт ярлыки/пункт в "Установка и удаление программ".
; ============================================================================

#define MyAppName "HelixFillDNA"
; ИСПРАВЛЕНИЕ БАГА "релиз без установщика" (v1.0.1): раньше здесь было
; безусловное #define MyAppVersion "1.0.0", которое ВСЕГДА перезаписывало
; значение, переданное из CI через /DMyAppVersion=... — ISCC применяет /D
; ДО обработки самого скрипта, поэтому обычный #define ниже просто затирал
; его обратно на "1.0.0". #ifndef оставляет /D-значение в приоритете и
; использует "1.0.0" только как запасной вариант при ручной сборке
; двойным кликом (когда /D никто не передавал).
#ifndef MyAppVersion
  #define MyAppVersion "1.0.0"
#endif
#define MyAppPublisher "HelixFillDNA"
#define MyAppExeName "HelixFillDNA.exe"
; Путь к папке, которую собрал PyInstaller (dist\HelixFillDNA рядом с этим .iss)
#define SourceDir "dist\HelixFillDNA"
; Текстовая версия EULA, показываемая на экране принятия лицензии.
; Файл должен лежать рядом с этим .iss (обычная кодировка, plain text
; или .rtf — Inno Setup умеет и то, и то; markdown он не понимает).
#define LicenseFilePath "EULA_installer.txt"

[Setup]
AppId={{A1B2C3D4-1234-5678-9ABC-DEF012345678}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
; Ставим НЕ в системную папку по умолчанию для обычных пользователей без
; прав администратора можно сменить на {localappdata}\{#MyAppName}
PrivilegesRequired=lowest
OutputDir=Output
OutputBaseFilename=HelixFillDNA-Setup-{#MyAppVersion}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
; Иконка самого установщика (Setup.exe) и он же по умолчанию используется
; для ярлыков, если не переопределить IconFilename ниже в [Icons].
SetupIconFile=app_icon.ico
ArchitecturesInstallIn64BitMode=x64
; Показывает пользователю текст EULA с двумя радиокнопками
; ("Я принимаю условия" / "Я не принимаю условия") ДО ввода пути
; установки. Кнопка "Далее" неактивна, пока не выбран пункт принятия —
; Inno Setup делает это сам, никакого кода писать не нужно.
LicenseFile={#LicenseFilePath}

[Languages]
Name: "russian"; MessagesFile: "compiler:Languages\Russian.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
; Копируем ВСЁ содержимое папки, собранной PyInstaller'ом (exe, DLL,
; внутренние библиотеки, bin\ с bcftools/tabix/bgzip и т.д.) рекурсивно.
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; Иконка для ярлыков (меню Пуск / рабочий стол) — сам exe уже содержит
; иконку, зашитую PyInstaller'ом (см. icon= в .spec), но отдельный файл
; тоже кладём в {app}, чтобы IconFilename в [Icons] ниже мог на него сослаться.
Source: "app_icon.ico"; DestDir: "{app}"; Flags: ignoreversion
; Промт "обычные / продвинутые настройки": трафареты для сборки итогового
; файла. В обычном режиме приложение подставляет их само по выбранному
; источнику данных (FTDNA -> template_v3.txt, MyHeritage -> template_v5.txt),
; поэтому после установки они должны лежать в {app}\samples — берём их из
; папки samples\ рядом с этим .iss, а не из dist\ (PyInstaller кладёт свои
; datas то плоско, то в _internal\, и полагаться на это нельзя).
; skipifsourcedoesntexist: если трафарета нет у сборщика, установщик всё
; равно соберётся, а приложение покажет подсказку, куда положить файл.
Source: "samples\template_v3.txt"; DestDir: "{app}\samples"; Flags: ignoreversion skipifsourcedoesntexist
Source: "samples\template_v5.txt"; DestDir: "{app}\samples"; Flags: ignoreversion skipifsourcedoesntexist
Source: "samples\README.txt"; DestDir: "{app}\samples"; Flags: ignoreversion skipifsourcedoesntexist

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\app_icon.ico"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\app_icon.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#MyAppName}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; На всякий случай подчищаем папки, которые приложение создаёт само во
; время работы (donors/output/reference — см. main.py, там
; mkdir(parents=True, exist_ok=True)), если пользователь захочет удалить
; вместе с накопленными данными. Если хотите ОСТАВЛЯТЬ скачанные доноры/
; референс при удалении — просто удалите этот блок.
; Type: filesandordirs; Name: "{app}\donors"
; Type: filesandordirs; Name: "{app}\output"
; Type: filesandordirs; Name: "{app}\reference"
; Маленький файл с состоянием интерфейса (выбранный режим настроек
; "Обычные"/"Продвинутые") — создаётся программой, при удалении не нужен.
Type: files; Name: "{app}\ui_state.json"
