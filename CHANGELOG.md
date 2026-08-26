# Changelog

Формат основан на [Keep a Changelog](https://keepachangelog.com/ru/1.1.0/),
версионирование — по [Semantic Versioning](https://semver.org/lang/ru/).

## [Unreleased]

### Planned
- CI-сборка установщика на push git-тега (GitHub Actions).
- `PRIVACY.md` и сайт-витрина проекта.

## [1.0.0] - 2026-08-26

### Added
- GUI (CustomTkinter) и CLI-пайплайн для конвертации сырых данных
  FTDNA / MyHeritage / готового VCF в формат 23andMe.
- Автозагрузка референсного генома с проверкой целостности (SHA-256).
- Поддержка паролей Michigan Imputation Server со спецсимволами.
- Раздельное хранение донорских образцов по источникам и диагностика
  пересечения чипов.
- Автопредложение и автозагрузка донорских образцов из 1000 Genomes
  прямо из GUI, с прогрессом, отменой и повтором при сетевой ошибке.
- Автодетект несоответствия выбранного источника и содержимого файла
  (VCF / FTDNA / MyHeritage).
- Установщик для Windows: сборка exe через PyInstaller и мастер
  установки через Inno Setup (`HelixFillDNA-Setup-1.0.0.exe`), без
  необходимости прав администратора.
