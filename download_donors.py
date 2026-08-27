"""
download_donors.py v11 — Задача 1 промта "Автопредложение скачать доноров
через GUI": вызываемая функция download_donors_for_chip() без модульных
глобалов + автоматическая инвалидация устаревшего кэша при смене чипа
одного и того же источника + поддержка отмены на любом шаге.

Изменения v10 -> v11:

  download_donors_for_chip(positions_json, source, output_dir, htslib,
                            progress_cb=None, cancel_check=None) -> list[Path]
      Новая единая точка входа для вызова ИЗ ПРИЛОЖЕНИЯ (gui/app.py),
      без запуска отдельного процесса/argparse. Замена цепочки
      "main() читает --positions-json -> create_chip_positions_from_json
      -> process_chromosome в цикле -> запись chip_signature.txt" в виде
      обычной функции с прогрессом и отменой.

  Убраны модульные глобалы HTSLIB и _WORKING_SUFFIX_BY_MIRROR.
      Раньше main() присваивал `global HTSLIB`, и вся логика ниже неявно
      читала его из области видимости модуля — при повторных вызовах
      скачивания доноров в рамках ОДНОГО процесса (как это происходит в
      GUI, который является долгоживущим процессом, а не разовым CLI-
      запуском) это давало неявную зависимость от предыдущего запуска.
      Теперь htslib и working_suffix_by_mirror — обычные параметры/
      локальные переменные, которые каждый вызов create'ит заново.

  ⚠ Инвалидация устаревшего кэша при смене чипа (Задача 1, п.7).
      process_chromosome() и раньше, и сейчас пропускает скачивание
      хромосомы, если kgp_sub_{chrom}.vcf.gz уже лежит на диске —
      НЕЗАВИСИМО от того, какому чипу этот файл принадлежит. Раньше это
      было безопасно, потому что download_donors.py запускался вручную
      пользователем как отдельный CLI-скрипт (пользователь сам следил за
      output-dir/--donors-subdir). Как только скачивание доноров стало
      доступно из GUI по одной кнопке "Да, скачать" на одну и ту же
      donors/<source>/ директорию для РАЗНЫХ чипов одного источника
      (например, сначала FTDNA-чип A, потом FTDNA-чип B), баг
      "Invalid alleles / union из двух чипов" (см. main.py, Задача A/B)
      возвращается через чёрный ход: process_chromosome() тихо возьмёт
      donora от чипа A как "уже готового" для чипа B.

      Поэтому download_donors_for_chip() ПЕРЕД циклом по хромосомам
      сравнивает chip_signature.txt (если есть в output_dir) с
      сигнатурой, извлечённой из positions_json. При несовпадении —
      старые kgp_sub_*.vcf.gz/служебные файлы удаляются, и все 22
      хромосомы перекачиваются/перефильтровываются заново под новый чип.

  Отмена скачивания на любом шаге (Задача 1, п.5, вариант (а)).
      subprocess.run(...) (блокирующий, ждёт весь curl целиком) заменён
      на Popen + периодический poll() с проверкой cancel_check() — при
      отмене процесс curl получает terminate()/kill(), а не ждётся до
      конца текущего файла. download_with_urllib() тоже проверяет
      cancel_check() между чанками.

  CLI (`python download_donors.py ...`) сохранён без изменений в
  интерфейсе командной строки — новый флаг не добавлен, поведение
  --csv / --positions-json / --donors-subdir идентично v10 (для
  --positions-json путь теперь идёт через download_donors_for_chip, что
  является улучшением: инвалидация устаревшего кэша работает и из CLI).

=============================================================================
Изменения v11 -> v12 (промт "Ускорение и переиспользование скачивания
доноров", Часть 1 — ускорение):
=============================================================================

  1.1. Удалённая фильтрация без полного скачивания хромосомы
      (process_chromosome_remote() + _probe_bcftools_remote_support()).
      Раньше process_chromosome() ВСЕГДА качал весь VCF хромосомы
      (для chr1 — несколько ГБ на 2504 образца), а затем локально
      фильтровал его bcftools view -S/-R до ~30-40 тыс позиций чипа —
      то есть трафик тратился на файл целиком ради долей процента от
      него. Если установленный bcftools/htslib собран с libcurl (умеет
      делать HTTP Range-запросы), можно сразу вызвать
      `bcftools view -S eur20.txt -R pos.txt "<URL>" -Oz -o kgp_sub_N.vcf.gz`
      — htslib сам подтянет удалённый .tbi-индекс и скачает ТОЛЬКО те
      байтовые блоки, которые реально попадают под фильтр.

      Поддержка проверяется ОДИН раз перед циклом по хромосомам
      (_probe_bcftools_remote_support(), пробный `bcftools view -h` на
      самой маленькой хромосоме, chr21, с таймаутом) — если сборка без
      libcurl или проба не удалась, remote_capable=False, и весь прогон
      идёт по-старому (process_chromosome, полное скачивание).

      Даже при remote_capable=True КАЖДАЯ хромосома, если удалённая
      фильтрация для неё не удалась (сетевая ошибка, битое зеркало,
      таймаут), тихо откатывается на process_chromosome() — см.
      process_chromosome_auto(). Отказ одной хромосомы не роняет весь
      прогон и не отключает удалённый путь для остальных.

  1.2. aria2c как более быстрый способ ПОЛНОГО скачивания
      (download_with_aria2c()) — используется, когда aria2c установлен
      в системе, как приоритетный способ ДО curl/urllib в download_file().
      Несколько параллельных TCP-соединений на один файл (-x/-s), плюс
      докачка (-c) — как в оригинальном ручном гайде проекта
      (Гайд_по_конвертации.docx, Часть 2/5.3), только теперь встроено в
      сам download_donors.py, а не требует ручных команд.

  1.3. Параллельная обработка нескольких хромосом одновременно
      (download_donors_for_chip() теперь гоняет ThreadPoolExecutor с
      max_parallel_chromosomes воркерами вместо строгого `for chrom in
      range(1, 23)` одна за другой). Хромосомы полностью независимы
      (разные временные файлы, разные kgp_sub_N.vcf.gz), поэтому
      единственное ограничение — пропускная способность канала и число
      одновременных соединений, а не искусственная последовательность.
      cancel_check() проверяется и до запуска задачи, и внутри самих
      сетевых операций (через _run_cancelable) — отмена на любом шаге
      останавливает приём новых хромосом и быстро гасит уже запущенные.

  1.4. Докачка в download_with_urllib() (fallback, если НИ aria2c, НИ
      curl не установлены). Раньше при каждом вызове файл стирался и
      качался с нуля (dest.unlink() в начале функции) — на
      многогигабайтном VCF хромосомы это означало полную потерю
      прогресса при любом обрыве соединения. Теперь используется тот же
      приём, что и в main.py::_download_with_resume() для референсного
      генома — Range-запрос с текущим размером файла на диске.

  Цепочка приоритетов на хромосому:
      remote random access (1.1, если поддерживается)
        -> aria2c (1.2, если установлен)
        -> curl (докачка -C -)
        -> urllib (докачка через Range, 1.4)
      — каждый следующий способ используется, только если предыдущий
      недоступен в системе или не сработал.

  Новые необязательные параметры download_donors_for_chip():
      remote_filter: bool = True — можно принудительно выключить
          удалённую фильтрацию (например, для отладки или если она
          заведомо не поддерживается инфраструктурой).
      max_parallel_chromosomes: int = DEFAULT_PARALLEL_CHROMOSOMES (3) —
          сколько хромосом обрабатывать одновременно.
      Оба параметра также доступны в CLI как --no-remote-filter и
      --parallel-chromosomes (см. main()). Поведение по умолчанию без
      явных флагов уже использует ускорения — старое строго
      последовательное полное скачивание остаётся доступно через
      --no-remote-filter --parallel-chromosomes 1.

  Задача B/D (проверка и надёжность реюза доноров между людьми одного
  чипа, chip_signature_broad) в этой версии НЕ затронута — это Часть 2
  того же промта, реализуется отдельно.

=============================================================================
Изменения v12 -> v13 (промт "Диагностика реальной удалённой фильтрации +
устойчивая настройка CA-сертификатов libcurl"):
=============================================================================

Предыстория: на реальной Windows-машине _probe_bcftools_remote_support()
изначально падала с SSL-ошибкой (сборка bcftools/libcurl не использует
системное хранилище сертификатов Windows), а после ручной установки
CURL_CA_BUNDLE — второй проблемой оказался собственный curl.exe в bin/
(конфликтует с системным curl.exe, если bin_dir добавлен в начало PATH).
Оба фикса раньше делались вручную в консоли и не переживали новую сессию
терминала/запуск из GUI — теперь это делает core/network_utils.py
(ensure_network_ready()), вызываемый один раз при старте в main.py и
gui/app.py (см. изменения там).

  diagnose_remote_filter(htslib, output_dir, test_chrom=21) -> dict
      Часть 1 промта: раньше единственной проверкой удалённого доступа
      был _probe_bcftools_remote_support() — лёгкий `bcftools view -h`
      БЕЗ ТЕЛА VCF. Он подтверждает, что libcurl технически может открыть
      URL, но не проверяет, что `-S eur20.txt -R pos.txt` (та же команда,
      что реально выполняет process_chromosome_remote() на каждой из 22
      хромосом) отдаёт непустой валидный результат — например, если
      удалённый .tbi индекс скачался, но битый, или Range-запросы к телу
      файла обрываются на середине (в отличие от короткого запроса
      заголовка), проба -h всё равно может пройти успешно, а реальная
      фильтрация — нет.

      Функция строит минимальный eur20.txt (переиспользует
      create_eur20_list()) и pos.txt с несколькими реальными позициями
      test_chrom, вызывает ровно ту же команду bcftools, что и
      process_chromosome_remote() (не process_chromosome_remote()
      напрямую — она пишет kgp_sub_{chrom}.vcf.gz в output_dir, а
      диагностика должна быть безопасна для запуска ДО реального прогона
      донор-скачивания и не портить/не путать с уже существующим кэшем
      доноров, поэтому пишет во временный файл и сама его удаляет), во
      временный .vcf.gz, считает записи через `bcftools view -H | wc`
      (эквивалент подсчёта строк из stdout), возвращает отчёт и НИКОГДА
      не бросает исключение — это диагностика, а не критический путь.

  CLI-флаг --diagnose-remote (в download_donors.py::main()) — печатает
      отчёт diagnose_remote_filter() в консоль и завершается с кодом
      0 (успех) / 1 (неудача), не запуская основное скачивание.

  which_curl_ignoring_dir()/warn_if_conflicting_curl() из
      core/network_utils.py используются здесь вместо голого
      shutil.which("curl") в _downloader_chain() — гарантирует, что
      найденный "curl" не окажется потенциально битым бандловым
      curl.exe из --bin-dir (см. докстринг network_utils.py).

=============================================================================
Изменения v13 -> v14 (промт "точечный патч remote-фильтрации для крупных
хромосом"):
=============================================================================

Предыстория: реальный прогон download_donors_for_chip() на 758 990
позициях чипа (Путь A, verify_path_a.py) показал, что удалённая
фильтрация (Часть 1.1, process_chromosome_remote()), уверенно
подтверждённая диагностикой (--diagnose-remote) на маленькой chr21,
на КРУПНЫХ хромосомах (в наблюдавшемся прогоне — chr1, chr2, chr3)
систематически проваливается на всех зеркалах/суффиксах с ошибками
вида `[E::bgzf_read_block] Failed to read BGZF header at offset ...`
и `[E::hts_itr_next] Failed to seek to offset ...: Invalid seek` —
Range-доступ вглубь большого файла у соответствующих зеркал нестабилен,
в отличие от коротких запросов к маленьким файлам/заголовкам.

  _is_transient_range_error(stderr) -> bool + _TRANSIENT_RANGE_ERROR_MARKERS
      Отличает "файла с этим суффиксом нет" (404/No such file or
      directory — нормальная ситуация, есть смысл пробовать другой
      суффикс на том же зеркале) от системной ошибки Range-доступа
      зеркала к ЭТОМУ файлу (bgzf_read_block/Invalid seek — повторный
      суффикс на том же зеркале почти наверняка провалится так же).

  process_chromosome_remote() при обнаружении системной ошибки сразу
      переходит к следующему ЗЕРКАЛУ, не тратя попытку на оставшийся
      суффикс текущего. Также теперь логирует суммарное время,
      потраченное на remote-попытки для каждой хромосомы (обнаружено:
      единицы секунд на провал, не таймаут — узкое место не в переборе
      самом по себе, а в последующем полном скачивании).

  DEFAULT_REMOTE_SKIP_LARGE_CHROMS = {1..8} + новый параметр
      process_chromosome_auto(..., skip_remote: bool = False) и
      download_donors_for_chip(..., remote_skip_large_chroms=...) —
      для хромосом из этого набора remote-путь не пробуется вовсе,
      сразу используется полное скачивание. {1,2,3} — подтверждено
      реальным прогоном, {4..8} — экстраполяция по размеру файлов
      1000 Genomes phase3 (тоже в числе самых крупных), требует
      подтверждения отдельным прогоном. Параметр полностью управляем:
      None/set() возвращает поведение v13 (пробовать remote для всех).
      Новый CLI-флаг --remote-skip-chroms (по умолчанию "1,2,3,4,5,6,7,8",
      пустая строка = поведение v13).

  Не тронуто: перебор для случая "файла нет", cancel_check()/progress_cb()
  на всех путях, формат chip_signature.txt/positions.json, защита от
  Invalid alleles (main.py, Задачи A/B/C) — вне scope патча.

=============================================================================
Изменения v14 -> v15 (мелкий патч "суффикс-специфичные временные файлы"):
=============================================================================

Подтверждено реальным прогоном ПОСЛЕ v14 (patch из этого файла сработал:
chr1/chr2/chr3 больше не перебирают 6 комбинаций зеркало×суффикс, remote-
путь для них пропускается сразу, как и задумано). Обнаружена отдельная,
не связанная с v14, мелкая неэффективность (см. download_chromosome_vcf()):
временный файл скачивания хромосомы был ОДИН на оба кандидата суффикса
(v5b/v5a), поэтому уже частично скачанные байты под один суффикс
подставлялись под попытку резюмирования по URL другого — на практике
это означало ~10-15 лишних секунд на крупную хромосому (2 попытки по 3с
ожидания под заведомо неверным суффиксом v5b, прежде чем дойти до
верного v5a), но НЕ приводило к потере данных или порче докачки.

  _suffix_temp_path(dest, suffix) — новая вспомогательная функция:
      строит суффикс-специфичное имя временного файла. download_chromosome_vcf()
      теперь качает/резюмирует каждый суффикс в свой файл и переименовывает
      его в исходный `dest` только при подтверждённом успехе (файл + индекс
      прошли verify_file()) — остальной код (process_chromosome() и далее)
      продолжает работать с прежним `dest`, интерфейс не изменился.

  Одноразовая миграция файлов старого (суффикс-неспецифичного) формата,
      оставшихся на диске с прошлых запусков (до этого патча): такой файл
      переносится под последний кандидат VCF_SUFFIX_CANDIDATES (v5a — по
      данным реального прогона именно он существует для крупных хромосом
      на всех трёх зеркалах), чтобы уже скачанные гигабайты не тратились
      впустую под заведомо неверный суффикс при первой же попытке.

=============================================================================
Изменения v15 -> v16 (промт "Доноры для VCF-источника: понятная отмена +
общий кэш сырых хромосом"):
=============================================================================

Предыстория (см. main.py/gui/app.py для UX-части того же промта): каждый
источник (ftdna/myheritage/vcf) хранит СВОЙ отфильтрованный кэш доноров
(donors/<source>/<panel>/kgp_sub_*.vcf.gz) — это правильно и не меняется
(позиции фильтрации разные). Но ПОЛНАЯ, ещё не отфильтрованная хромосома
1000 Genomes (ALL.chr{N}...vcf.gz, по несколько ГБ каждая), которую
process_chromosome() качает при неудаче удалённой фильтрации (Часть 1.1) —
это ОДИН И ТОТ ЖЕ файл независимо от источника/чипа (различаются только
позиции, по которым потом идёт локальная bcftools-фильтрация). Раньше этот
временный файл (`ALL.chr{chrom}.download.vcf.gz`) удалялся в конце
process_chromosome() — значит, второй источник (например, переключение
FTDNA -> VCF для того же человека, или просто другой чип FTDNA) заново
качал те же самые многогигабайтные файлы с зеркал 1000 Genomes.

  raw_cache_dir: Optional[Path] — новый необязательный параметр в
      download_chromosome_vcf() / process_chromosome() /
      process_chromosome_auto() / download_donors_for_chip(). По
      умолчанию None (обратная совместимость: поведение полностью
      совпадает с v15 — временный файл хромосомы по-прежнему исчезает из
      output_dir после фильтрации, просто теперь без общего кэша он
      никуда не копируется).

  Реализация сделана внутри download_chromosome_vcf() (а не
  process_chromosome()), потому что именно там в момент успеха известен
  конкретный суффикс (v5a/v5b), под которым файл был скачан — process_
  chromosome() снаружи этой информации не имеет (dest уже переименован).
  Функционально результат тот же, что описан в промте (полный файл
  хромосомы переиспользуется между источниками), только физически кладётся
  на уровень ниже:
    - _load_from_raw_cache(): ПЕРЕД обращением к зеркалам проверяет
      donors/_raw_chromosomes/<genome_build>/ALL.chr{N}.<suffix>.vcf.gz(.tbi)
      для каждого suffix из VCF_SUFFIX_CANDIDATES — если валидный файл
      найден, он захардлинкивается (или копируется, если хардлинк
      невозможен — другой раздел диска) в `dest`/`dest.tbi`, и сеть вообще
      не используется.
    - _store_in_raw_cache(): сразу после успешного скачивания (до
      локальной bcftools-фильтрации) хардлинкает (или копирует) готовый
      `dest`/`dest.tbi` в общий кэш под тем же suffix-специфичным именем —
      так следующий источник/чип найдёт его через _load_from_raw_cache().
    - Хардлинк (os.link) предпочитается копированию, если raw_cache_dir и
      output_dir на одном разделе диска — не удваивает место на диске за
      многогигабайтные файлы; при ошибке (другой раздел/ФС без хардлинков)
      тихо откатывается на shutil.copy2().

  main.py/gui/app.py передают путь donors/_raw_chromosomes/<genome_build>/
  (единая точка — main.py::raw_chromosome_cache_dir()) только если
  пользователь явно включил соответствующий чекбокс/флаг — по умолчанию
  выключено, экономия места на диске приоритетнее незапрошенного
  использования десятков ГБ.

  Не затронуто: process_chromosome_remote()/удалённая фильтрация (Часть
  1.1/v14) — общий кэш актуален только для пути ПОЛНОГО скачивания.
"""
from __future__ import annotations
import concurrent.futures
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Set, Tuple

from core.network_utils import (
    ensure_network_ready, which_curl_ignoring_dir, warn_if_conflicting_curl,
)

MIRRORS = [
    "https://1000genomes.s3.amazonaws.com/release/20130502/",
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    "http://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/",
]
SAMPLES_FILENAME = "integrated_call_samples_v3.20130502.ALL.panel"

# ---------------------------------------------------------------------------
# Промт "Monomorphic sites / настраиваемое количество EUR-доноров":
# раньше число доноров было жёстко зашито как 20 (create_eur20_list()
# всегда брала первые 20 строк с panel[2]=='EUR'). На маленькой выборке
# многие сайты, полиморфные в популяции в целом, случайно оказываются
# мономорфными во всех 20 взятых образцах — MIS такие сайты исключает из
# QC ("Monomorphic sites"), что снижает итоговое покрытие/call rate.
#
# EUR_SAMPLE_COUNT_ALL — сигнальное значение "взять всю доступную EUR-
# подвыборку панели вместо фиксированных 20" (по умолчанию именно оно —
# см. DEFAULT_EUR_SAMPLE_COUNT). Это НЕ подбор конкретных образцов под
# конкретные позиции чипа (такой "умный" greedy-подбор сознательно не
# реализован — он оптимизирует итоговые метрики QC под конкретный чип
# точечным выбором семпловой подвыборки, а не даёт более представительную
# случайную выборку, что имеет разные последствия для достоверности
# результата) — только честное увеличение размера случайной выборки.
EUR_SAMPLE_COUNT_ALL = None
DEFAULT_EUR_SAMPLE_COUNT = EUR_SAMPLE_COUNT_ALL
# Разумный потолок для валидации в CLI/GUI — примерный размер EUR-
# подвыборки 1000 Genomes phase3 (GRCh37); настоящий предел зависит от
# факта попавших в файл панели строк и берётся оттуда же.
MAX_EUR_SAMPLE_COUNT = 503
VCF_SUFFIX_CANDIDATES = ["v5b", "v5a"]
VCF_TEMPLATE = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_{suffix}.20130502.genotypes.vcf.gz"

# ---------------------------------------------------------------------------
# Промт "TopMed/HRC" (Промт_TopMed_HRC_v2.md), п.3: константы для GRCh38-
# релиза 1000 Genomes — используется панелью TopMed.
#
# ⚠ Проверено при реализации по нескольким независимым источникам
# (документация GATK/UCSC-мирроринг/публикации по alt-aware выравниванию
# на GRCh38, а также прямые команды bcftools view на этот URL из чужих
# постов на biostars) — датасет "20190312_biallelic_SNV_and_INDEL"
# (интегрированный фазированный SNV+INDEL коллсет 1000 Genomes на
# GRCh38, ~2548 образцов). В ОТЛИЧИЕ от GRCh37-релиза (3 зеркала — S3,
# EBI FTP, NCBI FTP — и два кандидата суффикса v5a/v5b) для этой
# конкретной GRCh38-директории подтвердить более одного зеркала не
# удалось: S3-бакет 1000genomes и NCBI ftp-trace хостят другие GRCh38-
# датасеты (сам референсный .fa, high-coverage 3202-сэмплов и т.д.), но
# не подтверждённо — именно эту директорию. Поэтому GRCH38_MIRRORS
# сейчас содержит ОДНО зеркало — официальный EBI FTP, откуда данные
# раздаются и на который ссылается сама документация 1000 Genomes.
# Это снижает отказоустойчивость по сравнению с GRCh37-путём (при
# недоступности единственного зеркала откат некуда) — если для боевого
# использования TopMed понадобится больше зеркал, их нужно будет найти
# и проверить отдельно перед добавлением сюда.
#
# Суффикс здесь тоже, в отличие от GRCh37, ОДИН подтверждённый вариант
# ("v2a_27022019") — списочная структура GRCH38_VCF_SUFFIX_CANDIDATES
# сохранена ради единообразия с download_chromosome_vcf()/
# process_chromosome_remote() (которые перебирают список суффиксов),
# но переборе там фактически не из чего выбирать.
GRCH38_MIRRORS = [
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/",
]
GRCH38_VCF_SUFFIX_CANDIDATES = ["v2a_27022019"]
GRCH38_VCF_TEMPLATE = "ALL.chr{chrom}.shapeit2_integrated_snvindels_{suffix}.GRCh38.phased.vcf.gz"

# Единая точка выбора набора констант по genome_build ("grch37" | "grch38")
# — все функции, читающие MIRRORS/VCF_TEMPLATE/VCF_SUFFIX_CANDIDATES
# напрямую, теперь принимают параметр genome_build и обращаются через эти
# словари, а не к голым модульным константам. Старые модульные константы
# (MIRRORS/VCF_SUFFIX_CANDIDATES/VCF_TEMPLATE без суффикса _BY_BUILD)
# оставлены как есть — это конфигурация GRCh37 по умолчанию, обратная
# совместимость с кодом (в т.ч. вне этого модуля), который может к ним
# обращаться напрямую.
MIRRORS_BY_BUILD: dict[str, list[str]] = {
    "grch37": MIRRORS,
    "grch38": GRCH38_MIRRORS,
}
VCF_SUFFIX_CANDIDATES_BY_BUILD: dict[str, list[str]] = {
    "grch37": VCF_SUFFIX_CANDIDATES,
    "grch38": GRCH38_VCF_SUFFIX_CANDIDATES,
}
VCF_TEMPLATE_BY_BUILD: dict[str, str] = {
    "grch37": VCF_TEMPLATE,
    "grch38": GRCH38_VCF_TEMPLATE,
}
DEFAULT_GENOME_BUILD = "grch37"


def _mirrors_for_build(genome_build: str) -> list[str]:
    try:
        return MIRRORS_BY_BUILD[genome_build]
    except KeyError:
        raise RuntimeError(
            f"Неизвестная генетическая сборка: {genome_build!r}. "
            f"Доступные: {', '.join(MIRRORS_BY_BUILD.keys())}"
        )


def _vcf_suffix_candidates_for_build(genome_build: str) -> list[str]:
    return VCF_SUFFIX_CANDIDATES_BY_BUILD.get(genome_build, VCF_SUFFIX_CANDIDATES)


def _vcf_template_for_build(genome_build: str) -> str:
    return VCF_TEMPLATE_BY_BUILD.get(genome_build, VCF_TEMPLATE)


# Промт "TopMed/HRC" (Промт_TopMed_HRC_v2.md), п.3, продолжение: для
# GRCh38-релиза 1000 Genomes значение колонки CHROM внутри самого VCF
# обычно записано с префиксом "chr" (chr1, chr2, ..., chrX), тогда как
# GRCh37-релиз (phase3) — без префикса (1, 2, ..., X). Позиции чипа
# (chip_positions) всегда хранятся в каноническом виде БЕЗ префикса (см.
# adapters/*.py::_normalize_chrom()) — при записи pos_chr{chrom}.txt для
# bcftools view -R нужно подставить префикс, который реально стоит в
# CHROM-колонке той сборки, которую мы читаем, иначе -R не найдёт ни
# одной позиции (тихо вернёт пустой результат, не ошибку). Согласовано
# с REFERENCE_PANELS["<panel>"]["chrom_prefix"] в main.py, но продублировано
# здесь как отдельная маленькая таблица — download_donors.py не импортирует
# main.py (нет такой зависимости и раньше, вносить её ради одной строки
# не хочется), а сама таблица достаточно стабильна (привязана к
# genome_build, а не к panel), чтобы дублирование не создавало риска
# рассинхронизации.
GENOME_BUILD_CHROM_PREFIX: dict[str, str] = {
    "grch37": "",
    "grch38": "chr",
}


def _chrom_label_for_build(chrom: int, genome_build: str) -> str:
    """Строка для колонки CHROM в pos_chr{chrom}.txt под данную сборку —
    chrom здесь всегда int (1..22), префикс "chr" подставляется только
    для GRCh38 (см. GENOME_BUILD_CHROM_PREFIX)."""
    prefix = GENOME_BUILD_CHROM_PREFIX.get(genome_build, "")
    return f"{prefix}{chrom}"


MAX_RETRIES = 3
IS_WINDOWS = os.name == "nt"

# --- Часть 1: ускорение скачивания (v12) ------------------------------
# Проба удалённого доступа bcftools делает лёгкий запрос (только заголовок
# VCF, без тела) на самую маленькую хромосому — секунды даже на медленном
# канале, но на случай совсем плохой сети/мёртвого зеркала всё равно нужен
# таймаут, чтобы не подвесить весь запуск.
REMOTE_PROBE_TIMEOUT = 20  # секунд
REMOTE_PROBE_CHROM = 21    # самая маленькая аутосома 1000 Genomes phase3

# Таймаут на одну попытку удалённой фильтрации одной хромосомы. Это не
# скачивание всего файла (см. 1.1), а точечные Range-запросы, но при очень
# большом количестве позиций чипа/плохой сети хочет разумный потолок,
# чтобы process_chromosome_auto() вовремя откатился на полное скачивание,
# а не завис на удалённом пути навсегда.
REMOTE_CHROM_TIMEOUT = 900  # секунд (15 минут) на попытку

ARIA2C_CONNECTIONS = 8  # столько же, сколько в ручном гайде проекта (aria2c -x 8 -s 8)

# Сколько хромосом обрабатывать одновременно (Часть 1.3). Хромосомы
# независимы, единственное реальное ограничение — пропускная способность
# канала и число одновременных соединений; 3 — консервативный дефолт,
# балансирующий скорость и не создающий слишком много одновременных
# subprocess/сетевых соединений на слабых машинах.
DEFAULT_PARALLEL_CHROMOSOMES = 3

# v14, Шаг 2 промта "точечный патч remote-фильтрации для крупных хромосом":
# хромосомы, для которых remote-путь (1.1) по умолчанию пропускается и
# сразу используется полное скачивание (process_chromosome(), которое
# само уже ускорено через aria2c/curl/докачку — см. 1.2/1.4).
#
# ⚠ Подтверждено реальным прогоном на 758 990 позициях чипа: chr1, chr2,
# chr3 систематически проваливают remote-фильтрацию на ВСЕХ зеркалах и
# суффиксах (bgzf_read_block/Invalid seek — см. _is_transient_range_error),
# каждая попытка быстрая (единицы секунд), но 6 комбинаций зеркало×суффикс
# на хромосому всё равно суммарно тратят время и сетевые запросы впустую
# перед неизбежным откатом. chr4-8 в список включены по аналогии (тоже
# входят в число самых крупных файлов 1000 Genomes phase3, счёт на ГБ) —
# ⚠ ЭТО ЭКСТРАПОЛЯЦИЯ, не подтверждена отдельным прогоном; при появлении
# данных по chr4-8 (или по другой сети/зеркалам, где remote-путь может
# отработать иначе) значение стоит пересчитать.
#
# Параметр полностью управляем вызывающим кодом (remote_skip_large_chroms
# в download_donors_for_chip) — можно передать None/пустое множество,
# чтобы принудительно пробовать remote-путь для всех хромосом как в v13,
# или свой список по факту наблюдений на конкретной сети.
# v17 (патч "убрать бесполезный перебор suffix/зеркал для chr9-11"):
# реальный прогон (лог пользователя, myheritage/hrc, 584 023 позиции чипа)
# показал, что chr9, chr10, chr11 ТОЖЕ систематически проваливают
# remote-фильтрацию на всех зеркалах/суффиксах с теми же самыми
# симптомами (bgzf_read_block/Invalid seek/таймаут), что и chr1-8 — но, в
# отличие от них, КАЖДАЯ такая попытка не быстрая (секунды), а тратит
# ПОЛНЫЙ REMOTE_CHROM_TIMEOUT (900с) на зеркалах, где происходит таймаут
# соединения, а не явная ошибка (см. "chr10: на remote-попытки потрачено
# 1040.1с", "chr11: ... 1019.6с" в логе) — то есть почти по 17 минут
# впустую на каждую хромосому перед неизбежным откатом на полное
# скачивание. Добавлены в список пропуска по факту наблюдения.
#
# v18 (реальный прогон, ftdna/hrc, 755 902 позиции чипа): chr14-chr18
# ТОЖЕ систематически проваливают remote-фильтрацию на всех зеркалах и
# суффиксах — та же картина, что и раньше с chr1-11 (bgzf_read_block/
# Invalid seek/BCF read error/таймаут), с теми же дорогими накладными
# расходами: "chr14: на remote-попытки потрачено 2701.8с" (45 минут!),
# "chr17: ... 1144.2с" (19 минут) — прежде чем откатиться на полное
# скачивание, которое затем сработало быстро и надёжно. chr12, chr13,
# chr19-22 в этом прогоне отдельно не проверялись (лог не содержал их
# remote-попыток) — намеренно НЕ добавлены в список пропуска без
# подтверждения, по тому же принципу "только по факту наблюдения", что
# и раньше.
DEFAULT_REMOTE_SKIP_LARGE_CHROMS: Set[int] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18}

# Файлы, привязанные к конкретному чипу/сигнатуре — удаляются целиком
# при обнаружении несовпадения chip_signature.txt (Задача 1, п.7).
# Промт "настраиваемое количество EUR-доноров": имя файла со списком
# образцов теперь зависит от eur_sample_count ("eur20.txt", "eur120.txt",
# "eur503.txt", ...) — фиксированный кортеж больше не может перечислить
# все варианты явно, такие файлы чистятся отдельно через
# glob("eur*.txt") в _invalidate_stale_donor_cache() ниже.
_SIGNATURE_SCOPED_STATIC_FILES = ("ftdna_pos.txt", "chip_signature.txt")
# Оставлено для обратной совместимости с любым внешним кодом, который мог
# импортировать старое имя константы напрямую.
_SIGNATURE_SCOPED_FILES = _SIGNATURE_SCOPED_STATIC_FILES


class DownloadCancelled(RuntimeError):
    """Скачивание было прервано пользователем через cancel_check()."""


class HtslibTools:
    def __init__(self, bin_dir: Optional[Path]):
        self.bin_dir = bin_dir
        self.bgzip_path = self._find("bgzip")
        self.tabix_path = self._find("tabix")
        self.bcftools_path = self._find("bcftools")

    def _find(self, name: str) -> Optional[str]:
        exe = name + (".exe" if IS_WINDOWS else "")
        if self.bin_dir:
            candidate = self.bin_dir / exe
            if candidate.is_file():
                return str(candidate)
        return shutil.which(name)

    @property
    def has_bgzip(self) -> bool: return self.bgzip_path is not None
    @property
    def has_tabix(self) -> bool: return self.tabix_path is not None
    @property
    def has_bcftools(self) -> bool: return self.bcftools_path is not None


def check_dependencies(htslib: HtslibTools) -> None:
    """Проверяет наличие bcftools/tabix. htslib передаётся параметром —
    больше не читается из модульного глобала (см. докстринг файла)."""
    if not htslib.has_bcftools:
        print("\n" + "=" * 70)
        print("ОШИБКА: bcftools не найден")
        print("=" * 70)
        print("Укажите папку с бинарниками через --bin-dir или добавьте её в PATH.")
        sys.exit(1)
    if not htslib.has_tabix:
        print("⚠ tabix не найден — индексация будет невозможна")
    print("✓ bcftools найден")
    if htslib.has_tabix:
        print("✓ tabix найден")


class ProgressBar:
    def __init__(self, total: int, label: str = "", width: int = 40):
        self.total, self.current, self.label, self.width = total, 0, label, width

    def update(self, n: int = 1):
        self.current += n
        pct = self.current / self.total * 100 if self.total > 0 else 0
        filled = int(self.width * self.current / self.total) if self.total > 0 else self.width
        bar = '█' * filled + '░' * (self.width - filled)
        print(f'\r{self.label}: [{bar}] {pct:.1f}% ({self.current}/{self.total})', end='', flush=True)
        if self.current >= self.total:
            print()

    def finish(self):
        if self.current < self.total:
            self.current = self.total
            self.update(0)


# ---------------------------------------------------------------------------
# Скачивание одного файла — с поддержкой отмены (Задача 1, п.5, вариант (а))
# ---------------------------------------------------------------------------
# Промт "скачивание должно быть видно в приложении, а не в отдельных
# окнах cmd": в оконной сборке PyInstaller (console=False) у GUI-процесса
# нет своей консоли, поэтому Windows создаёт НОВОЕ консольное окно для
# каждого запускаемого консольного приложения (aria2c.exe/curl.exe) — на
# экране появлялись отдельные чёрные окна с прогрессом, а сам прогресс
# при этом не попадал в лог приложения. CREATE_NO_WINDOW подавляет
# создание такого окна, а перехват stdout/stderr позволяет показывать
# прогресс скачивания прямо в логе приложения (см. _pump_progress).
_CREATE_NO_WINDOW = 0x08000000 if os.name == "nt" else 0

# Прогресс aria2c/curl обновляется возвратом каретки ("\r") много раз в
# секунду. Печатать каждое такое обновление в лог нельзя — это тысячи
# строк на файл, лог станет нечитаемым и начнёт тормозить GUI. Поэтому
# показываем не чаще, чем раз в _PROGRESS_EMIT_INTERVAL секунд.
_PROGRESS_EMIT_INTERVAL = 2.0


def _pump_progress(stream, label: str, emit: Callable[[str], None],
                    interval: float = _PROGRESS_EMIT_INTERVAL) -> None:
    """
    Читает вывод дочернего процесса (aria2c/curl) и отдаёт его в лог
    приложения через emit(), прореживая частые обновления прогресса.

    Читаем побайтово-порциями и делим не только по "\n", но и по "\r":
    строка прогресса у aria2c/curl перерисовывается возвратом каретки и
    без такого деления пришла бы одним бесконечно растущим куском.

    Никогда не бросает исключений наружу: поток-читатель не должен ронять
    скачивание, что бы ни пришло из дочернего процесса.
    """
    last_emit = 0.0
    buf = ""
    try:
        while True:
            data = stream.read(256)
            if not data:
                break
            if isinstance(data, bytes):
                data = data.decode("utf-8", errors="replace")
            buf += data
            parts = buf.replace("\r", "\n").split("\n")
            buf = parts.pop()
            for line in parts:
                line = line.strip()
                if not line:
                    continue
                now = time.monotonic()
                # Строки прогресса ("[#abc 1.0GiB/1.2GiB(83%) ...]")
                # прореживаем по времени, а всё остальное (ошибки,
                # предупреждения aria2c) показываем сразу.
                is_progress = line.startswith("[#") or "ETA:" in line
                if is_progress and (now - last_emit) < interval:
                    continue
                if is_progress:
                    last_emit = now
                emit(f"    {label}: {line}")
    except Exception:
        pass
    finally:
        try:
            stream.close()
        except Exception:
            pass


def _run_cancelable(cmd: list[str], cancel_check: Optional[Callable[[], bool]] = None,
                     poll_interval: float = 0.5,
                     progress_label: Optional[str] = None,
                     **popen_kwargs) -> int:
    """
    Запускает subprocess через Popen и периодически опрашивает его,
    вместо subprocess.run(), который блокирует поток до полного
    завершения процесса. Это единственный способ реально прервать уже
    начавшееся скачивание одной хромосомы, а не только между хромосомами
    (см. докстринг файла, "Отмена скачивания на любом шаге").

    Бросает DownloadCancelled, если cancel_check() вернул True — процесс
    при этом получает terminate(), а при необходимости kill().
    """
    # Подавляем создание отдельного консольного окна на Windows и, если
    # запрошен показ прогресса, перехватываем вывод дочернего процесса.
    if _CREATE_NO_WINDOW:
        popen_kwargs["creationflags"] = (
            popen_kwargs.get("creationflags", 0) | _CREATE_NO_WINDOW
        )
    pump: Optional[threading.Thread] = None
    if progress_label and "stdout" not in popen_kwargs:
        popen_kwargs["stdout"] = subprocess.PIPE
        popen_kwargs["stderr"] = subprocess.STDOUT

    proc = subprocess.Popen(cmd, **popen_kwargs)
    if progress_label and proc.stdout is not None:
        pump = threading.Thread(
            target=_pump_progress,
            args=(proc.stdout, progress_label, lambda s: print(s, flush=True)),
            daemon=True,
        )
        pump.start()
    try:
        while True:
            try:
                return proc.wait(timeout=poll_interval)
            except subprocess.TimeoutExpired:
                if cancel_check and cancel_check():
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                    raise DownloadCancelled("Скачивание отменено пользователем")
    except DownloadCancelled:
        raise
    except BaseException:
        # Любая другая причина прерывания (например, KeyboardInterrupt) —
        # тоже гасим дочерний процесс, а не оставляем его висеть.
        proc.kill()
        raise


def download_with_aria2c(url: str, dest: Path,
                          cancel_check: Optional[Callable[[], bool]] = None) -> bool:
    """
    Часть 1.2: скачивание с несколькими параллельными TCP-соединениями на
    один файл (aria2c -x/-s), как в ручном гайде проекта (aria2c -x 8 -s 8
    -c) — заметно быстрее одиночного curl/urllib на больших VCF-файлах
    1000 Genomes, особенно на каналах с высокой задержкой. Поддерживает
    докачку (-c), как и download_with_curl().

    Возвращает False (не бросает), если aria2c не установлен в системе —
    вызывающая сторона (download_file) просто пробует следующий способ.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "aria2c",
        "-x", str(ARIA2C_CONNECTIONS), "-s", str(ARIA2C_CONNECTIONS),
        "-c", "--allow-overwrite=true", "--console-log-level=warn",
        "--summary-interval=0",
        "-d", str(dest.parent), "-o", dest.name,
        url,
    ]
    try:
        ret = _run_cancelable(cmd, cancel_check=cancel_check,
                               progress_label=f"aria2c {dest.name}")
        return ret == 0
    except FileNotFoundError:
        return False


def download_with_curl(url: str, dest: Path,
                        cancel_check: Optional[Callable[[], bool]] = None) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-L", "--fail", "--retry", "0", "-C", "-", "--progress-bar", "-o", str(dest), url]
    try:
        ret = _run_cancelable(cmd, cancel_check=cancel_check,
                               progress_label=f"curl {dest.name}")
        return ret == 0
    except FileNotFoundError:
        return False


def download_with_urllib(url: str, dest: Path,
                          cancel_check: Optional[Callable[[], bool]] = None) -> bool:
    """
    Часть 1.4: скачивание через встроенный urllib — последний fallback,
    если НИ aria2c, НИ curl не установлены в системе. Раньше при каждом
    вызове файл стирался и качался заново (dest.unlink() в начале
    функции), что на многогигабайтном VCF хромосомы означало полную
    потерю прогресса при любом обрыве соединения. Теперь поддерживает
    докачку через Range-запрос — та же идея, что и в
    main.py::_download_with_resume() для референсного генома.
    """
    import urllib.request
    import urllib.error

    existing_size = dest.stat().st_size if dest.exists() else 0
    headers = {'User-Agent': 'Mozilla/5.0'}
    if existing_size > 0:
        headers['Range'] = f'bytes={existing_size}-'
        print(f"  Найден частично скачанный файл ({existing_size / 1024**2:.1f} МБ), докачиваю...")

    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=120) as response:
            resumed = getattr(response, "status", 200) == 206
            if existing_size > 0 and not resumed:
                # Сервер не поддержал Range и вернул файл целиком (200) —
                # начинаем с нуля, чтобы не задвоить уже скачанные байты.
                existing_size = 0
            content_length = int(response.headers.get('Content-Length', 0))
            total = existing_size + content_length
            bar = ProgressBar(total, f"  {dest.name}")
            bar.current = existing_size
            mode = 'ab' if (resumed and existing_size > 0) else 'wb'
            with open(dest, mode) as f:
                while True:
                    if cancel_check and cancel_check():
                        raise DownloadCancelled("Скачивание отменено пользователем")
                    chunk = response.read(65536)
                    if not chunk:
                        break
                    f.write(chunk)
                    bar.update(len(chunk))
            bar.finish()
        return True
    except DownloadCancelled:
        raise
    except urllib.error.HTTPError as e:
        if e.code == 416:
            # Requested Range Not Satisfiable — обычно значит, что файл
            # уже скачан полностью (существующий размер >= размера на сервере).
            print(f"\n  ✓ {dest.name}: похоже, уже скачан полностью (416)")
            return True
        print(f"\n  ✗ urllib ошибка: {e.code} {e.reason}")
        return False
    except Exception as e:
        print(f"\n  ✗ urllib ошибка: {e}")
        return False


def _downloader_chain(bin_dir: Optional[Path] = None) -> list[Callable[..., bool]]:
    """
    Часть 1.2/1.4: цепочка способов ПОЛНОГО скачивания файла, по убыванию
    скорости — aria2c (несколько соединений на файл, если установлен) ->
    curl (если установлен) -> urllib (всегда доступен, с докачкой).
    Формируется заново на каждый вызов download_file(), а не кэшируется
    модульно — shutil.which() достаточно дёшев, а актуальность (вдруг
    aria2c поставили/удалили между вызовами в рамках одного долгоживущего
    GUI-процесса) важнее микрооптимизации.

    bin_dir: если задан (обычно --bin-dir с bcftools/tabix), поиск curl
        игнорирует его через which_curl_ignoring_dir() — на Windows там
        может лежать собственный curl.exe в комплекте с htslib-бандлом,
        у которого зашит нерабочий относительный путь к сертификатам
        (см. докстринг core/network_utils.py и v13-заметку в начале
        этого файла) — нам нужен системный curl, а не он.
    """
    chain: list[Callable[..., bool]] = []
    if shutil.which("aria2c"):
        chain.append(download_with_aria2c)
    if which_curl_ignoring_dir(bin_dir):
        chain.append(download_with_curl)
    chain.append(download_with_urllib)
    return chain


def download_file(url: str, dest: Path,
                   cancel_check: Optional[Callable[[], bool]] = None,
                   bin_dir: Optional[Path] = None) -> bool:
    """
    Пробует способы скачивания по очереди (см. _downloader_chain()) —
    переход к следующему способу происходит, если предыдущий недоступен
    в системе ИЛИ явно вернул неудачу. Файл на диске между попытками не
    трогаем: каждый способ сам умеет докачивать уже частично скачанные
    байты (-c у aria2c/curl, Range у urllib), так что смена инструмента
    на середине не теряет прогресс.
    """
    for downloader in _downloader_chain(bin_dir=bin_dir):
        if downloader(url, dest, cancel_check=cancel_check):
            return True
    return False


def verify_file(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    if path.name.endswith('.gz') or path.name.endswith('.tbi'):
        try:
            with gzip.open(path, 'rb') as f:
                while True:
                    if not f.read(65536): break
            return True
        except Exception:
            return False
    try:
        with path.open('r', encoding='utf-8') as f:
            f.read(1024)
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# v18: быстрая проверка целостности BGZF без полной распаковки — используется
# ТОЛЬКО для _load_from_raw_cache() (общий кэш сырых хромосом).
# ---------------------------------------------------------------------------
# Стандартный 28-байтовый BGZF EOF-маркер (пустой BGZF-блок, которым htslib
# завершает КАЖДЫЙ корректно дозаписанный .gz/.tbi файл — задокументирован
# в спецификации SAM/BGZF). Его отсутствие в последних 28 байтах надёжно
# указывает на оборванную/неполную докачку — ровно то, что нам и нужно
# отловить, не читая файл целиком.
_BGZF_EOF_MARKER = bytes([
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
])


def verify_file_fast(path: Path) -> bool:
    """
    Облегчённая проверка BGZF-файла (.gz/.tbi) — вместо распаковки
    гигабайтного файла целиком (verify_file(), дорого при КАЖДОМ
    обращении к общему кэшу сырых хромосом — на chr1 это реально
    заметное время на каждый запуск, даже когда сеть вообще не
    участвует) читает только первые несколько байт (gzip/BGZF-сигнатура)
    и последние 28 байт (стандартный BGZF EOF-маркер).

    Это НЕ полноценная проверка CRC каждого блока (для этого пришлось бы
    распаковать файл — то, чего мы и избегаем), но она надёжно отличает
    полностью и корректно дозаписанный файл (какими htslib/bgzip всегда
    завершают вывод) от оборванной докачки/усечённого файла — а именно
    это единственный реалистичный сценарий порчи файла, лежащего в
    общем кэше сырых хромосом (см. докстринг файла про докачку).

    Не используется как замена verify_file() везде — для файлов сразу
    после СЕТЕВОГО скачивания (где риск порчи середины потока реальнее)
    полная проверка остаётся в силе, см. download_chromosome_vcf().
    """
    try:
        size = path.stat().st_size
        if size < 28:
            return False
        with path.open('rb') as f:
            head = f.read(4)
            if head[:2] != b'\x1f\x8b':
                # даже не gzip — точно не наш случай, откатываемся на
                # полную проверку, чтобы не давать ложноположительный результат
                return verify_file(path)
            f.seek(-28, os.SEEK_END)
            tail = f.read(28)
        return tail == _BGZF_EOF_MARKER
    except OSError:
        return False


def is_404_error(dest: Path, download_ok: bool) -> bool:
    return (not download_ok) and (not dest.exists() or dest.stat().st_size == 0)


# ---------------------------------------------------------------------------
# Промт "проверка 'донор не пустой' после скачивания/фильтрации".
#
# Подтверждено реальным прогоном: kgp_sub_*.vcf.gz для panel="topmed"
# получались практически пустыми (корректный VCF-заголовок + список
# образцов, но 0 строк с вариантами) из-за нестабильного HTTPS-соединения
# при удалённой Range-фильтрации (process_chromosome_remote()) — в
# конкретном случае мешал активный VPN у пользователя в сочетании с
# отсутствующими в окружении CURL_CA_BUNDLE/SSL_CERT_FILE. Файл при этом
# получал ненулевой размер и bcftools возвращал код 0 (сам заголовок
# скачался и записался нормально), поэтому process_chromosome()/
# process_chromosome_remote() считали хромосому успешно обработанной, а
# check_donor_cache() потом принимал такой пустой кэш как валидный на
# каждом следующем запуске (проверяется только chip_signature.txt, не
# содержимое файлов).
# ---------------------------------------------------------------------------
def _count_vcf_records(bcftools_path: str, vcf_path: Path) -> int:
    """
    Считает число строк с вариантами в VCF (эквивалент
    `bcftools view -H <file> | wc -l`, как уже делает
    diagnose_remote_filter()).

    Возвращает -1, если саму проверку выполнить не удалось (bcftools
    отсутствует/упал/файл не открылся) — это НЕ то же самое, что "0
    записей": вызывающий код должен трактовать -1 как "проверка
    неубедительна, не отбраковываем файл только на этом основании", а 0
    — как подтверждённо пустой донор.
    """
    try:
        result = subprocess.run(
            [bcftools_path, "view", "-H", str(vcf_path)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            return -1
        return len([l for l in result.stdout.splitlines() if l.strip()])
    except Exception:
        return -1


_EMPTY_DONOR_HINT = (
    "вероятно, обрыв HTTPS-соединения при скачивании/фильтрации данных "
    "1000 Genomes — активный VPN или прокси нередко мешает стабильности "
    "такого соединения. Попробуйте отключить VPN/прокси и повторить."
)


def _vcf_has_any_record(bcftools_path: str, vcf_path: Path) -> Optional[bool]:
    """
    Дешёвая проверка "есть ли в VCF хотя бы одна строка с вариантом" —
    в отличие от _count_vcf_records() (который читает ВЕСЬ вывод
    bcftools целиком), читает только самую первую строку из потока и
    сразу останавливает процесс. Полный подсчёт по многогигабайтному
    файлу (полная нефильтрованная хромосома в общем кэше, десятки
    миллионов строк) был бы неоправданно дорогим для одной лишь проверки
    "пусто/не пусто" — а именно она и нужна для валидации закэшированных
    файлов ПЕРЕД их переиспользованием.

    Возвращает True/False, либо None, если саму проверку выполнить не
    удалось (bcftools не найден/ошибка запуска процесса) — трактуется
    вызывающим кодом так же, как -1 у _count_vcf_records(): "неубедительно,
    не считаем это признаком пустого файла".
    """
    try:
        proc = subprocess.Popen(
            [bcftools_path, "view", "-H", str(vcf_path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
    except Exception:
        return None
    try:
        first_line = proc.stdout.readline() if proc.stdout else ""
        return bool(first_line and first_line.strip())
    except Exception:
        return None
    finally:
        try:
            if proc.stdout:
                proc.stdout.close()
        except Exception:
            pass
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass


def _vcf_first_chrom(bcftools_path: str, vcf_path: Path) -> Optional[str]:
    """
    Дешёво возвращает значение колонки CHROM самой первой строки VCF (или
    None, если файл пуст/bcftools недоступен) — читает только первую
    строку вывода, не весь файл. Используется, чтобы определить, в каком
    виде донорский файл реально называет свою хромосому ("1" или "chr1"),
    не полагаясь на предположения о конкретном датасете 1000 Genomes.
    """
    try:
        proc = subprocess.Popen(
            [bcftools_path, "view", "-H", str(vcf_path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
    except Exception:
        return None
    try:
        first_line = proc.stdout.readline() if proc.stdout else ""
        if not first_line.strip():
            return None
        return first_line.split("\t", 1)[0]
    except Exception:
        return None
    finally:
        try:
            if proc.stdout:
                proc.stdout.close()
        except Exception:
            pass
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Промт "выровнять CHROM доноров с sample.vcf.gz для merge".
#
# Подтверждено реальным прогоном: GRCh38-датасет 1000 Genomes EBI
# (20190312_biallelic_SNV_and_INDEL), с которого качаются доноры для
# TopMed, хранит CHROM БЕЗ префикса "chr" ("##contig=<ID=1>"), тогда как
# main.py::build_vcf()/split_autosomes() пишут sample.vcf.gz этой же
# панели С префиксом "chr" (REFERENCE_PANELS["topmed"]["chrom_prefix"]).
# bcftools merge/concat сравнивают CHROM буквально — при несовпадении
# донор и sample.vcf.gz физически не пересекутся ни по одной позиции при
# слиянии, и main.py::_post_merge_intersect() на Этапе 6 отфильтрует всё
# подчистую ("N → 0 позиций"), даже когда сами доноры уже не пустые (см.
# фикс regions-файла в process_chromosome()/process_chromosome_remote()
# выше — тот фикс решает проблему ФИЛЬТРАЦИИ доноров по позициям чипа,
# этот — проблему их СОВМЕСТИМОСТИ с sample.vcf.gz при последующем merge;
# это два независимых проявления одного и того же ошибочного
# предположения об именовании контигов конкретного датасета).
# ---------------------------------------------------------------------------
def _ensure_donor_chrom_prefix(
    htslib: HtslibTools, out_file: Path, chrom: int, genome_build: str,
) -> bool:
    """
    Приводит CHROM готового донорского файла (kgp_sub_{chrom}.vcf.gz) к
    виду, который GENOME_BUILD_CHROM_PREFIX[genome_build] ожидает от
    sample.vcf.gz этой же сборки — через `bcftools annotate
    --rename-chrs` с переиндексацией результата. Если файл уже в нужном
    виде (типичный случай для GRCh37/HRC, где префикса нет вовсе, и для
    любого GRCh38-датасета, который в будущем будет называть контиги как
    ожидается) — ничего не делает, просто возвращает True.

    Возвращает False только если переименование реально потребовалось,
    но не удалось (bcftools annotate/переиндексация завершились с
    ошибкой) — вызывающий код должен считать это неудачей хромосомы, так
    как несовпадающий CHROM всё равно сделал бы донора бесполезным при
    дальнейшем merge.
    """
    desired_prefix = GENOME_BUILD_CHROM_PREFIX.get(genome_build, "")
    if not desired_prefix:
        return True  # GRCh37/HRC — CHROM без префикса ожидается и так

    current_chrom = _vcf_first_chrom(htslib.bcftools_path, out_file)
    if current_chrom is None:
        # Не удалось определить (bcftools недоступен и т.п.) — не пытаемся
        # переименовывать то, что не можем прочитать; проверка "донор не
        # пустой" (record_count) уже отработала до вызова этой функции.
        return True

    desired_chrom = f"{desired_prefix}{chrom}"
    if current_chrom == desired_chrom:
        return True  # уже в нужном виде

    print(
        f"  ℹ chr{chrom}: донор называет свою хромосому {current_chrom!r}, "
        f"а sample.vcf.gz этой панели ожидает {desired_chrom!r} — привожу "
        f"CHROM донора в соответствие (bcftools annotate --rename-chrs), "
        f"иначе merge с sample.vcf.gz не найдёт ни одного совпадения."
    )

    rename_map = out_file.parent / f"_rename_chrs_{chrom}.txt"
    renamed_out = out_file.parent / f"kgp_sub_{chrom}.renamed.vcf.gz"
    try:
        rename_map.write_text(f"{current_chrom}\t{desired_chrom}\n", encoding="utf-8")
        cmd = [
            htslib.bcftools_path, "annotate",
            "--rename-chrs", str(rename_map),
            str(out_file), "-Oz", "-o", str(renamed_out),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(
                f"  ✗ chr{chrom}: не удалось переименовать CHROM донора "
                f"({(result.stderr or '').strip()[-300:] or 'нет вывода stderr'})"
            )
            renamed_out.unlink(missing_ok=True)
            return False

        renamed_out.replace(out_file)
        if htslib.has_tabix:
            out_file.with_suffix(out_file.suffix + ".tbi").unlink(missing_ok=True)
            tbi_res = subprocess.run(
                [htslib.tabix_path, "-p", "vcf", str(out_file)],
                capture_output=True, text=True,
            )
            if tbi_res.returncode != 0:
                print(f"  ✗ chr{chrom}: не удалось переиндексировать донора после переименования CHROM")
                return False
        print(f"  ✓ chr{chrom}: CHROM донора приведён к {desired_chrom!r}")
        return True
    finally:
        rename_map.unlink(missing_ok=True)
        renamed_out.unlink(missing_ok=True)


def _suffix_temp_path(dest: Path, suffix: str) -> Path:
    """
    v15 (мелкий патч по итогам реального прогона, см. докстринг
    download_chromosome_vcf ниже): строит суффикс-специфичный путь
    временного файла скачивания хромосомы, например
    "ALL.chr1.download.vcf.gz" -> "ALL.chr1.download.v5a.vcf.gz".

    Вставляем суффикс перед ".vcf.gz" вручную (а не через dest.stem/
    dest.suffix), потому что у ".vcf.gz" два "суффикса" с точки зрения
    pathlib (".gz" — dest.suffix, ".vcf" остаётся в dest.stem) — обычная
    конкатенация дала бы "ALL.chr1.download.vcf.v5a.gz", что и работало
    бы, но менее читаемо в логах/на диске, чем place-суффикса перед
    ".vcf.gz" целиком.
    """
    name = dest.name
    if name.endswith(".vcf.gz"):
        base = name[: -len(".vcf.gz")]
        return dest.with_name(f"{base}.{suffix}.vcf.gz")
    return dest.with_name(f"{name}.{suffix}")


# ---------------------------------------------------------------------------
# v16: общий кэш «сырых» (ещё не отфильтрованных) хромосом 1000 Genomes,
# переиспользуемый между всеми источниками/чипами одной генетической
# сборки (genome_build) — см. докстринг файла.
# ---------------------------------------------------------------------------
def _link_or_copy(src: Path, dst: Path) -> None:
    """
    Хардлинк src -> dst, если возможно (тот же раздел диска — не тратит
    место на многогигабайтные файлы), иначе обычное копирование
    (например, raw_cache_dir и output_dir на разных ФС/разделах).
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


# ---------------------------------------------------------------------------
# v19 ("ещё ускорение — кэш подвыборки 20 образцов"): полная хромосома
# 1000 Genomes содержит 2504 образца, из которых пайплайну нужны только
# 20 (eur20.txt) — раньше отбор этих 20 образцов и фильтрация по позициям
# конкретного чипа делались ОДНОЙ командой bcftools (-S + -R), что честно
# читает входной файл только один раз, НО каждый новый чип/источник
# заново читает и распаковывает ВЕСЬ файл хромосомы (сотни МБ - единицы
# ГБ), хотя после отбора образцов от него остаётся ~0.8% данных.
#
# Список 20 образцов детерминирован для конкретного релиза 1000 Genomes
# (create_eur20_list() всегда берёт первые 20 EUR из одного и того же
# файла панели) — то есть НЕ зависит от чипа/источника/человека, только
# от genome_build. Поэтому подвыборку можно построить один раз и
# переиспользовать для ЛЮБОГО следующего чипа/источника той же сборки —
# ровно та же идея, что уже применена к самим сырым хромосомам
# (_load_from_raw_cache/_store_in_raw_cache), только на следующем уровне.
#
# -m2 -M2 -v snps (биаллельные SNP) применяются уже на этапе подвыборки,
# а не после фильтрации по позициям — итоговый набор записей от этого не
# меняется (это независимые AND-условия, порядок не важен), а кэшируемый
# файл получается ещё компактнее.
# ---------------------------------------------------------------------------
def _eur_subset_cache_path(cache_dir: Path, chrom: int, eur_sample_count: int) -> Path:
    """
    Промт "настраиваемое количество EUR-доноров": имя файла кэша
    подвыборки теперь включает фактическое число образцов
    (f"EUR{N}.chr{chrom}.vcf.gz") — без этого подвыборка на 20 образцов
    молча использовалась бы как "уже готовая" и для запроса с другим
    eur_sample_count, давая неверный (заниженный/иной) набор доноров без
    единой ошибки в логе.
    """
    return Path(cache_dir) / f"EUR{eur_sample_count}.chr{chrom}.vcf.gz"


def _ensure_eur_subset(
    tmp_vcf: Path,
    eur_file: Path,
    chrom: int,
    htslib: HtslibTools,
    output_dir: Path,
    bcftools_threads: int = 1,
    raw_cache_dir: Optional[Path] = None,
    eur_sample_count: int = 20,
) -> Path:
    """
    Возвращает путь к VCF со всеми позициями хромосомы, но только
    eur_sample_count образцами (без фильтрации по чипу) — либо взятому
    из общего кэша (если он там уже есть от предыдущего чипа/источника
    той же сборки И того же eur_sample_count), либо только что
    построенному из tmp_vcf.

    eur_sample_count (промт "настраиваемое количество EUR-доноров"):
    сколько образцов реально в eur_file — используется только для имени
    файла кэша подвыборки (_eur_subset_cache_path()), чтобы разные
    размеры выборки не путали друг друга на диске. Сам список образцов в
    eur_file вызывающий код уже подготовил под нужный размер — здесь он
    не обрезается и не проверяется повторно.

    Если raw_cache_dir не задан — подвыборка строится во временный файл
    в output_dir и НЕ кэшируется (не с чем сравнивать в следующий раз,
    так как без общего кэша сырых хромосом каждый прогон и так качает
    хромосому заново — двухшаговая схема здесь бы только добавляла
    лишний проход по диску без пользы, поэтому вызывающий код в этом
    случае должен использовать старую однокомандную фильтрацию).
    """
    if raw_cache_dir is not None:
        cached = _eur_subset_cache_path(Path(raw_cache_dir), chrom, eur_sample_count)
        cached_tbi = cached.with_suffix(cached.suffix + ".tbi")
        if (cached.exists() and cached_tbi.exists()
                and verify_file_fast(cached) and verify_file_fast(cached_tbi)):
            # Промт "проверка 'донор не пустой'", продолжение: verify_file_fast()
            # проверяет только структурную целостность BGZF (заголовок +
            # корректный EOF-маркер) — файл, который однажды был построен
            # ИЗ ПУСТОГО/ОБОРВАННОГО источника (нестабильная сеть/VPN во
            # время предыдущего прогона), пройдёт эту проверку и будет
            # переиспользоваться НАВСЕГДА как якобы валидный кэш, даже
            # если в нём физически нет ни одной строки с вариантом — это
            # ровно то, что наблюдалось в реальном прогоне (все 22
            # хромосомы синхронно давали 0 записей именно потому, что
            # подвыборка бралась из одного и того же испорченного кэша
            # каждый раз). Дешёвая проверка "есть хотя бы одна запись"
            # ловит это, не читая многомиллионностроковый файл целиком.
            has_record = _vcf_has_any_record(htslib.bcftools_path, cached)
            if has_record is False:
                print(
                    f"  ⚠ chr{chrom}: закэшированная подвыборка ({cached.name}) "
                    f"структурно цела, но пуста (0 записей) — похоже, кэш был "
                    f"построен из оборванного источника в прошлый раз. Удаляю "
                    f"испорченный кэш и строю подвыборку заново."
                )
                cached.unlink(missing_ok=True)
                cached_tbi.unlink(missing_ok=True)
            else:
                # has_record is True, либо None (сама проверка не удалась —
                # не отбраковываем кэш только на этом основании).
                print(f"  ✓ chr{chrom}: подвыборка {eur_sample_count} образцов взята из кэша ({cached.name})")
                return cached

    subset_vcf = output_dir / f"EUR{eur_sample_count}.chr{chrom}.build.vcf.gz"
    cmd = [
        htslib.bcftools_path, "view",
        "-S", str(eur_file), "--force-samples",
        "-m2", "-M2", "-v", "snps",
        "--threads", str(max(1, bcftools_threads)),
        str(tmp_vcf),
        "-Oz", "-o", str(subset_vcf),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        subset_vcf.unlink(missing_ok=True)
        raise RuntimeError(
            f"bcftools view (подвыборка {eur_sample_count} образцов) для chr{chrom} "
            f"завершился с ошибкой:\n{result.stderr}"
        )

    if htslib.has_tabix:
        subprocess.run(
            [htslib.tabix_path, "-p", "vcf", str(subset_vcf)],
            check=True, capture_output=True,
        )

    if raw_cache_dir is not None:
        cached = _eur_subset_cache_path(Path(raw_cache_dir), chrom, eur_sample_count)
        cached_tbi = cached.with_suffix(cached.suffix + ".tbi")
        subset_tbi = subset_vcf.with_suffix(subset_vcf.suffix + ".tbi")
        try:
            _link_or_copy(subset_vcf, cached)
            if subset_tbi.exists():
                _link_or_copy(subset_tbi, cached_tbi)
            print(f"  ✓ chr{chrom}: подвыборка {eur_sample_count} образцов сохранена в кэш ({cached})")
            # Данные уже в кэше (хардлинк — тот же inode, либо копия) —
            # временная копия в output_dir больше не нужна и не должна
            # оставаться мусором после завершения обработки хромосомы.
            subset_vcf.unlink(missing_ok=True)
            subset_tbi.unlink(missing_ok=True)
        except OSError as e:
            print(f"  ⚠ Не удалось сохранить подвыборку {eur_sample_count} образцов в кэш ({raw_cache_dir}): {e}")
            # Не удалось закэшировать — используем то, что уже построено
            # локально, вызывающий код сам вычистит его как обычный
            # временный файл в output_dir.
            return subset_vcf
        return cached

    return subset_vcf


def _raw_cache_has_chrom(raw_cache_dir: Path, chrom: int, genome_build: str = DEFAULT_GENOME_BUILD) -> bool:
    """
    v17: быстрая ПРОВЕРКА (без хардлинка/копирования — этим занимается
    _load_from_raw_cache()), есть ли валидная нефильтрованная хромосома в
    общем кэше сырых хромосом. Если да — пробовать remote-фильтрацию для
    неё вообще не имеет смысла: process_chromosome() (полное скачивание)
    и так почти мгновенно возьмёт файл из этого же кэша через
    download_chromosome_vcf() -> _load_from_raw_cache(), а remote-путь
    может потратить впустую до REMOTE_CHROM_TIMEOUT секунд на каждое
    зеркало (см. докстринг DEFAULT_REMOTE_SKIP_LARGE_CHROMS про
    chr9-11 — до ~17 минут на хромосому).

    genome_build (промт "TopMed/HRC", п.3): определяет, по какому шаблону
    имени файла искать в кэше — raw_cache_dir уже раздельный по сборке
    (см. main.py::raw_chromosome_cache_dir()), но имя файла внутри него
    всё равно должно соответствовать шаблону/суффиксам именно этой
    сборки, иначе поиск ничего не найдёт даже при валидном кэше.
    """
    raw_cache_dir = Path(raw_cache_dir)
    template = _vcf_template_for_build(genome_build)
    for suffix in _vcf_suffix_candidates_for_build(genome_build):
        cached_vcf = raw_cache_dir / template.format(chrom=chrom, suffix=suffix)
        cached_tbi = cached_vcf.with_suffix(cached_vcf.suffix + ".tbi")
        if cached_vcf.exists() and cached_tbi.exists() and verify_file_fast(cached_vcf) and verify_file_fast(cached_tbi):
            return True
    return False


def _load_from_raw_cache(
    raw_cache_dir: Path, chrom: int, dest: Path, tbi_dest: Path,
    genome_build: str = DEFAULT_GENOME_BUILD,
    bcftools_path: Optional[str] = None,
) -> Optional[str]:
    """
    Проверяет donors/_raw_chromosomes/<genome_build>/ на наличие уже
    скачанной (кем-то другим — другим источником/чипом) полной хромосомы.
    Если найдена валидная пара (файл + .tbi) под любым из суффиксов,
    актуальных для genome_build — линкует/копирует её в dest/tbi_dest и
    возвращает найденный суффикс. Иначе — None (сеть/зеркала как раньше).

    bcftools_path (промт "проверка 'донор не пустой'", продолжение):
    если передан, дополнительно дешёво проверяет (_vcf_has_any_record —
    только первая строка вывода, без полного прохода по многомиллионным
    записям), что закэшированный файл не пуст, ПЕРЕД тем как отдать его
    вызывающему коду как "готово". verify_file_fast() выше проверяет
    только структурную целостность BGZF (корректный EOF-маркер) — файл,
    однажды построенный из оборванного/неполного скачивания (нестабильная
    сеть/VPN), проходит эту проверку, но может физически не содержать ни
    одной записи, и тогда переиспользовался бы НАВСЕГДА как валидный. При
    обнаружении такого — испорченные файлы кэша удаляются, и функция
    пробует следующий суффикс/зеркало вместо немедленного возврата.
    Если bcftools_path не передан (обратная совместимость) — проверка
    пропускается, поведение как раньше.
    """
    template = _vcf_template_for_build(genome_build)
    for suffix in _vcf_suffix_candidates_for_build(genome_build):
        cached_vcf = raw_cache_dir / template.format(chrom=chrom, suffix=suffix)
        cached_tbi = cached_vcf.with_suffix(cached_vcf.suffix + ".tbi")
        if cached_vcf.exists() and cached_tbi.exists() and verify_file_fast(cached_vcf) and verify_file_fast(cached_tbi):
            if bcftools_path:
                has_record = _vcf_has_any_record(bcftools_path, cached_vcf)
                if has_record is False:
                    print(
                        f"  ⚠ chr{chrom}: закэшированная полная хромосома "
                        f"({cached_vcf.name}) структурно цела, но пуста (0 "
                        f"записей) — похоже, кэш был построен из оборванного "
                        f"скачивания в прошлый раз. Удаляю испорченный кэш и "
                        f"пробую следующий вариант."
                    )
                    cached_vcf.unlink(missing_ok=True)
                    cached_tbi.unlink(missing_ok=True)
                    continue
            _link_or_copy(cached_vcf, dest)
            _link_or_copy(cached_tbi, tbi_dest)
            return suffix
    return None


def _store_in_raw_cache(
    raw_cache_dir: Path, chrom: int, suffix: str, dest: Path, tbi_dest: Path,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> None:
    """
    Сразу после успешного скачивания сохраняет копию (хардлинк, если
    возможно) полной хромосомы в общий кэш — для переиспользования
    следующим источником/чипом той же генетической сборки. Не бросает
    исключений наружу: сохранение в общий кэш — оптимизация, а не
    критический путь, поломка диска/прав доступа не должна ронять
    обработку текущей хромосомы.
    """
    try:
        raw_cache_dir = Path(raw_cache_dir)
        raw_cache_dir.mkdir(parents=True, exist_ok=True)
        template = _vcf_template_for_build(genome_build)
        cached_vcf = raw_cache_dir / template.format(chrom=chrom, suffix=suffix)
        cached_tbi = cached_vcf.with_suffix(cached_vcf.suffix + ".tbi")
        _link_or_copy(dest, cached_vcf)
        _link_or_copy(tbi_dest, cached_tbi)
    except OSError as e:
        print(f"  ⚠ Не удалось сохранить chr{chrom} в общий кэш сырых хромосом ({raw_cache_dir}): {e}")


def download_chromosome_vcf(
    chrom: int,
    dest: Path,
    htslib: HtslibTools,
    working_suffix_by_mirror: dict[str, str],
    cancel_check: Optional[Callable[[], bool]] = None,
    raw_cache_dir: Optional[Path] = None,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> bool:
    """
    genome_build (промт "TopMed/HRC", п.3): выбирает набор зеркал/
    шаблона имени файла/кандидатов суффикса через _mirrors_for_build()/
    _vcf_template_for_build()/_vcf_suffix_candidates_for_build() вместо
    жёстко зашитых модульных констант MIRRORS/VCF_TEMPLATE/
    VCF_SUFFIX_CANDIDATES (которые остаются алиасами GRCh37-конфигурации
    для обратной совместимости, см. их определение выше). Также передаётся
    в _load_from_raw_cache()/_store_in_raw_cache(), чтобы общий кэш сырых
    хромосом не путал файлы разных сборок под одинаковыми номерами
    хромосом (у него и так отдельная директория на диске на каждый
    genome_build — см. main.py::raw_chromosome_cache_dir() — но сама
    функция кэша должна знать, по какому шаблону/суффиксам искать).

    working_suffix_by_mirror передаётся вызывающим кодом явно (обычный
    dict, живущий ровно в рамках одного вызова download_donors_for_chip
    / одного прогона main()) — раньше это был модульный глобал
    _WORKING_SUFFIX_BY_MIRROR, из-за чего повторные вызовы в одном
    процессе неявно "помнили" рабочий суффикс с прошлого запуска.

    =========================================================================
    v15 — мелкий патч "суффикс-специфичные временные файлы" (по итогам
    реального прогона после v14, см. тест.txt/резюме пользователя):
    =========================================================================

    Раньше временный файл скачивания (`dest`, например
    "ALL.chr1.download.vcf.gz") был ОДИН на оба кандидата суффикса
    (v5b/v5a). При переборе суффиксов на одном зеркале уже частично
    скачанные байты от предыдущей попытки подставлялись под попытку
    докачки/resume по URL СЛЕДУЮЩЕГО (другого!) суффикса — в реальном
    прогоне это означало, что многогигабайтный частично скачанный файл
    под v5a пытался резюмироваться против v5b (которого для крупных
    хромосом не существует на большинстве зеркал), получал 404,
    `is_404_error()` не распознавал это как "файла нет" (потому что
    локальный файл был не пуст — там лежали чужие для этого URL байты),
    и код тратил 2 попытки × 3 секунды ожидания впустую, прежде чем
    дойти до суффикса v5a, где резюмирование наконец шло против верного
    URL и завершалось успешно. Данные при этом не портились — только
    терялось ~10-15 секунд на хромосому.

    Исправление: у каждого суффикса теперь свой временный файл
    (`_suffix_temp_path(dest, suffix)`), поэтому резюмирование всегда
    идёт против ТОГО ЖЕ URL, что и создал файл. По завершении успешного
    суффикса этот файл переименовывается в исходный `dest` — остальной
    код (process_chromosome() и далее) продолжает работать с `dest`,
    как и раньше, без изменений интерфейса.

    Миграция уже имевшихся частично скачанных файлов старого (суффикс-
    неспецифичного) формата: если на диске остался `dest` с прошлого
    запуска (до этого патча) — переносим его под наиболее вероятный
    суффикс (последний в VCF_SUFFIX_CANDIDATES, т.е. "v5a" — по данным
    реального прогона именно он существует для крупных хромосом на всех
    трёх зеркалах), а не оставляем как мёртвый вес на диске / не тратим
    его впустую под заведомо неверный суффикс, как было раньше.
    """
    tbi_dest = dest.with_suffix(dest.suffix + ".tbi")

    if dest.exists() and dest.stat().st_size > 0 and tbi_dest.exists() and tbi_dest.stat().st_size > 0:
        if verify_file(dest) and verify_file(tbi_dest):
            print(f"  ✓ {dest.name} и его индекс уже существуют")
            return True
        else:
            print(f"  ⚠ Файл или индекс повреждены, удаляю и перекачиваю")
            dest.unlink(missing_ok=True)
            tbi_dest.unlink(missing_ok=True)

    # v15: одноразовая миграция частично скачанного файла старого формата
    # (без суффикса в имени) — см. докстринг выше. Выполняется здесь, а
    # не внутри цикла по зеркалам/суффиксам, потому что имя файла не
    # зависит от зеркала — достаточно перенести его один раз на весь
    # вызов функции.
    if dest.exists() and dest.stat().st_size > 0:
        likely_suffix = _vcf_suffix_candidates_for_build(genome_build)[-1]
        migrated_dest = _suffix_temp_path(dest, likely_suffix)
        if not migrated_dest.exists():
            print(
                f"  ℹ Найден частично скачанный файл старого формата "
                f"({dest.stat().st_size / 1024**2:.1f} МБ) — переношу под "
                f"суффикс {likely_suffix!r} (наиболее вероятный по прошлым "
                f"прогонам) вместо повторной попытки докачки под заведомо "
                f"неверным суффиксом."
            )
            dest.rename(migrated_dest)
            legacy_tbi = tbi_dest
            if legacy_tbi.exists():
                migrated_tbi = migrated_dest.with_suffix(migrated_dest.suffix + ".tbi")
                if not migrated_tbi.exists():
                    legacy_tbi.rename(migrated_tbi)
                else:
                    legacy_tbi.unlink(missing_ok=True)

    # v16: общий кэш сырых хромосом — проверяем ДО обращения к зеркалам.
    # Если полная хромосома уже была скачана для другого источника/чипа
    # (той же генетической сборки), берём её отсюда без единого сетевого
    # запроса. См. докстринг файла ("Изменения v15 -> v16").
    if raw_cache_dir is not None:
        raw_cache_dir = Path(raw_cache_dir)
        found_suffix = _load_from_raw_cache(
            raw_cache_dir, chrom, dest, tbi_dest, genome_build=genome_build,
            bcftools_path=htslib.bcftools_path,
        )
        if found_suffix:
            print(
                f"  ✓ chr{chrom}: взят из общего кэша сырых хромосом "
                f"({raw_cache_dir}, суффикс {found_suffix}) — без обращения "
                f"к зеркалам 1000 Genomes"
            )
            working_suffix_by_mirror.setdefault("_raw_cache", found_suffix)
            return True

    build_mirrors = _mirrors_for_build(genome_build)
    build_template = _vcf_template_for_build(genome_build)
    build_suffix_candidates = _vcf_suffix_candidates_for_build(genome_build)

    for mirror in build_mirrors:
        if cancel_check and cancel_check():
            raise DownloadCancelled("Скачивание отменено пользователем")

        mirror_name = mirror.split('/')[2]
        known_suffix = working_suffix_by_mirror.get(mirror)
        suffixes_to_try = [known_suffix] if known_suffix else build_suffix_candidates

        for suffix in suffixes_to_try:
            filename = build_template.format(chrom=chrom, suffix=suffix)
            url = mirror + filename
            # v15: суффикс-специфичный временный файл — см. докстринг выше.
            suffix_dest = _suffix_temp_path(dest, suffix)
            suffix_tbi_dest = suffix_dest.with_suffix(suffix_dest.suffix + ".tbi")
            print(f"\n  Зеркало: {mirror_name} "
                  f"{'(HTTPS)' if 'amazonaws' in mirror_name else ''} [суффикс {suffix}]")

            for attempt in range(1, MAX_RETRIES + 1):
                if cancel_check and cancel_check():
                    raise DownloadCancelled("Скачивание отменено пользователем")

                print(f"  Попытка {attempt}/{MAX_RETRIES}...")
                ok = download_file(url, suffix_dest, cancel_check=cancel_check, bin_dir=htslib.bin_dir)

                if ok and verify_file(suffix_dest):
                    print(f"  Скачиваю индекс...")
                    tbi_url = url + ".tbi"
                    tbi_ok = download_file(tbi_url, suffix_tbi_dest, cancel_check=cancel_check, bin_dir=htslib.bin_dir)

                    if tbi_ok and verify_file(suffix_tbi_dest):
                        print(f"  ✓ {dest.name} и индекс скачаны (суффикс {suffix})")
                        suffix_dest.rename(dest)
                        suffix_tbi_dest.rename(tbi_dest)
                        working_suffix_by_mirror[mirror] = suffix
                        if raw_cache_dir is not None:
                            _store_in_raw_cache(raw_cache_dir, chrom, suffix, dest, tbi_dest, genome_build=genome_build)
                        return True
                    else:
                        print(f"  ⚠ Индекс не скачался, создаю локально через tabix...")
                        suffix_tbi_dest.unlink(missing_ok=True)
                        if htslib.has_tabix:
                            res = subprocess.run(
                                [htslib.tabix_path, "-p", "vcf", "-f", str(suffix_dest)],
                                capture_output=True
                            )
                            if res.returncode == 0 and suffix_tbi_dest.exists():
                                print(f"  ✓ Индекс создан локально")
                                suffix_dest.rename(dest)
                                suffix_tbi_dest.rename(tbi_dest)
                                working_suffix_by_mirror[mirror] = suffix
                                if raw_cache_dir is not None:
                                    _store_in_raw_cache(raw_cache_dir, chrom, suffix, dest, tbi_dest, genome_build=genome_build)
                                return True
                            else:
                                print(f"  ✗ Не удалось создать индекс локально")
                        else:
                            print(f"  ✗ Индекс не скачался, а tabix недоступен")

                        suffix_dest.unlink(missing_ok=True)
                        suffix_tbi_dest.unlink(missing_ok=True)
                        break

                elif is_404_error(suffix_dest, ok):
                    print(f"  ✗ Файл с суффиксом {suffix} не найден (404)")
                    break
                else:
                    print(f"  ✗ Ошибка скачивания")

                if attempt < MAX_RETRIES:
                    print(f"  Жду 3 секунды...")
                    time.sleep(3)

            suffix_dest.unlink(missing_ok=True)
            suffix_tbi_dest.unlink(missing_ok=True)
        print(f"  ✗ Зеркало {mirror_name} не помогло")

    print(f"  ✗ Не удалось скачать chr{chrom}")
    return False


def create_eur_donor_list(
    output_dir: Path,
    bin_dir: Optional[Path] = None,
    eur_sample_count: Optional[int] = DEFAULT_EUR_SAMPLE_COUNT,
) -> list[str]:
    """
    ⚠ Исправление (проверка на соответствие Задаче 1): раньше при ошибке
    скачивания панели образцов здесь вызывался sys.exit(1). Это было
    безопасно для одноразового CLI-скрипта, но download_donors_for_chip()
    вызывает эту функцию напрямую из ДОЛГОЖИВУЩЕГО GUI-процесса —
    sys.exit(1) в этом контексте убивает всё приложение целиком, а не
    только скачивание доноров, и ломает диалог "Повторить?" (Задача 1,
    п.6), потому что RuntimeError/DownloadCancelled просто не успевают
    добраться до except в _run_stages_1_6(). Теперь бросаем RuntimeError —
    main() (CLI) как и раньше ловит его и завершается с ненулевым кодом
    через sys.exit(str(e)), а download_donors_for_chip() пробрасывает его
    вызывающему коду как обычное исключение.

    eur_sample_count (промт "Monomorphic sites / настраиваемое количество
        EUR-доноров"): сколько строк с panel[2]=='EUR' взять из файла
        панели образцов.
          None (по умолчанию, EUR_SAMPLE_COUNT_ALL) — взять ВСЮ доступную
              EUR-подвыборку панели (обычно порядка 500 человек в
              phase3) — уменьшает долю Monomorphic sites на QC MIS по
              сравнению со старым фиксированным поведением "первые 20",
              ценой большего объёма скачивания/фильтрации на каждую
              хромосому.
          int — явный лимит (например, 20 — прежнее поведение, 100 —
              компромисс по трафику/скорости). Если в файле панели EUR-
              строк меньше запрошенного числа — берутся все найденные, с
              предупреждением в лог (это не ошибка).

    Имя файла со списком образцов на диске зависит от фактического числа
    отобранных образцов (f"eur{N}.txt") — так разные размеры выборки в
    разных запусках не перезаписывают и не путают друг друга на диске
    (см. также _invalidate_stale_donor_cache(), которая чистит такие
    файлы через glob("eur*.txt"), а не по одному фиксированному имени).

    bin_dir: передаётся в download_file()/which_curl_ignoring_dir(), чтобы
        не подцепить потенциально битый curl.exe из папки бинарников htslib
        (см. v13-заметку в начале файла и core/network_utils.py).
    """
    label = "все доступные" if eur_sample_count is None else str(eur_sample_count)
    print(f"\n[1/3] Список европейских образцов (запрошено: {label})...")
    panel_file = output_dir / SAMPLES_FILENAME
    if not panel_file.exists():
        if not download_file(MIRRORS[0] + SAMPLES_FILENAME, panel_file, bin_dir=bin_dir):
            raise RuntimeError(
                f"Не удалось скачать список образцов 1000 Genomes "
                f"({SAMPLES_FILENAME}) ни с одного зеркала. Проверьте "
                f"подключение к интернету и повторите попытку."
            )

    european = []
    with panel_file.open('r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[2] == 'EUR':
                european.append(parts[0])
                if eur_sample_count is not None and len(european) >= eur_sample_count:
                    break

    if eur_sample_count is not None and len(european) < eur_sample_count:
        print(
            f"  ⚠ В файле панели найдено только {len(european)} образцов EUR "
            f"(запрошено {eur_sample_count}) — используются все найденные."
        )

    eur_file = output_dir / f"eur{len(european)}.txt"
    with eur_file.open('w', newline='\n') as f:
        for s in european:
            f.write(s + '\n')
    print(f"  ✓ Отобрано {len(european)} европейских образцов → {eur_file.name}")
    return european


def create_eur20_list(output_dir: Path, bin_dir: Optional[Path] = None) -> list[str]:
    """
    Устаревшая обёртка для обратной совместимости с любым внешним кодом,
    который мог импортировать create_eur20_list() напрямую (жёстко 20
    образцов) — новый код должен использовать create_eur_donor_list()
    с явным eur_sample_count.
    """
    return create_eur_donor_list(output_dir, bin_dir=bin_dir, eur_sample_count=20)


def _extract_signature_from_positions_json(positions_json: Path) -> Optional[str]:
    """
    save_position_cache() в adapters/*.py пишет файл как
    '<chip_signature>.positions.json'. Извлекаем сигнатуру из имени файла,
    а не пересчитываем её заново — так гарантированно совпадает с тем,
    что main.py/gui/app.py положат в chip_signature.txt для сравнения.
    """
    name = positions_json.name
    suffix = ".positions.json"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return None


def create_chip_positions_from_json(positions_json: Path, output_dir: Path) -> Set[Tuple[str, int]]:
    """
    Вместо повторного парсинга сырого CSV, читаем позиции чипа из JSON,
    сохранённого save_position_cache() тем же самым парсером
    (ftdna_v3/myheritage_v5/vcf_source), который уже отработал в основном
    пайплайне.
    """
    print(f"\n[2/3] Позиции чипа (из {positions_json.name})...")
    payload = json.loads(positions_json.read_text(encoding="utf-8"))
    positions: Set[Tuple[str, int]] = {(str(chrom), int(pos)) for chrom, pos in payload}

    pos_file = output_dir / "ftdna_pos.txt"
    with pos_file.open('w', newline='\n') as f:
        for c, p in sorted(positions):
            f.write(f"{c}\t{p}\n")
    print(f"  ✓ {len(positions)} уникальных позиций → {pos_file.name}")
    return positions


def create_chip_positions(csv_path: Path, output_dir: Path) -> Set[Tuple[str, int]]:
    print("\n[2/3] Позиции чипа...")
    positions: Set[Tuple[str, int]] = set()
    with csv_path.open('r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if len(row) != 4: continue
            chrom, pos_str, result = row[1].strip(), row[2].strip(), row[3].strip().upper()
            if result == '--': continue
            chrom = chrom.replace('chr', '')
            if chrom == 'XY': chrom = 'X'
            if chrom == 'M': chrom = 'MT'
            try:
                positions.add((chrom, int(pos_str)))
            except ValueError:
                continue

    pos_file = output_dir / "ftdna_pos.txt"
    with pos_file.open('w', newline='\n') as f:
        for c, p in sorted(positions):
            f.write(f"{c}\t{p}\n")
    print(f"  ✓ {len(positions)} уникальных позиций → {pos_file.name}")
    return positions


def process_chromosome(
    chrom: int,
    eur_samples: Set[str],
    chip_positions: Set[Tuple[str, int]],
    output_dir: Path,
    htslib: HtslibTools,
    working_suffix_by_mirror: dict[str, str],
    cancel_check: Optional[Callable[[], bool]] = None,
    raw_cache_dir: Optional[Path] = None,
    bcftools_threads: int = 1,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> bool:
    """
    genome_build (промт "TopMed/HRC", п.3): прокидывается в
    download_chromosome_vcf() (выбор зеркал/шаблона/суффиксов и правильной
    директории общего кэша) и используется здесь же для формирования
    CHROM-колонки в pos_chr{chrom}.txt через _chrom_label_for_build() —
    GRCh38-релиз 1000 Genomes обычно хранит CHROM с префиксом "chr", в
    отличие от GRCh37 (см. GENOME_BUILD_CHROM_PREFIX выше). Без этого
    bcftools view -R тихо не находил бы ни одной позиции на GRCh38.

    raw_cache_dir (v16, см. докстринг файла "Изменения v15 -> v16"):
    прокидывается напрямую в download_chromosome_vcf(), которое и решает,
    брать ли полную хромосому из общего кэша сырых хромосом вместо
    зеркал, и сохранять ли туда свежескачанный файл — см. докстринг
    download_chromosome_vcf() про то, почему логика кэша живёт именно
    там, а не здесь.

    bcftools_threads (v18, "ускорение фильтрации уже скачанных файлов"):
    прокидывается в `bcftools view --threads N` — де/компрессия BGZF
    хорошо распараллеливается, а именно она (не сеть) — узкое место при
    фильтрации многогигабайтного файла хромосомы, взятого из общего
    кэша. Значение считает download_donors_for_chip() с учётом того,
    сколько хромосом уже обрабатывается параллельно (max_parallel_chromosomes),
    чтобы не переподписать процессор.
    """
    out_file = output_dir / f"kgp_sub_{chrom}.vcf.gz"
    if out_file.exists() and out_file.stat().st_size > 0:
        print(f"  ✓ chr{chrom} уже готов")
        return True

    if cancel_check and cancel_check():
        raise DownloadCancelled("Скачивание отменено пользователем")

    tmp_vcf = output_dir / f"ALL.chr{chrom}.download.vcf.gz"
    print(f"\n  --- chr{chrom} ---")

    if not download_chromosome_vcf(
        chrom, tmp_vcf, htslib, working_suffix_by_mirror, cancel_check,
        raw_cache_dir=raw_cache_dir, genome_build=genome_build,
    ):
        return False

    if cancel_check and cancel_check():
        raise DownloadCancelled("Скачивание отменено пользователем")

    eur_sample_count = len(eur_samples)
    print(f"  Фильтрую по {eur_sample_count} образцам и позициям чипа (bcftools)...")

    eur_file = output_dir / f"eur{eur_sample_count}_chr{chrom}.txt"
    with eur_file.open('w', newline='\n') as f:
        for s in eur_samples:
            f.write(s + '\n')

    # v18 (патч "почему фильтрация из кэша всё равно долгая"): позиции
    # раньше писались в порядке итерации Python set() — то есть НЕ
    # отсортированы по номеру позиции. bcftools view -R (regions-file)
    # для индексированного входного VCF использует random access через
    # .tbi, но при этом ожидает регионы в отсортированном порядке —
    # иначе htslib вынужден постоянно прыгать по индексированному файлу
    # туда-сюда вместо одного последовательного прохода вперёд, что на
    # десятках-сотнях тысяч разбросанных позиций (как в вашем чипе)
    # превращает фильтрацию в минуты вместо секунд, даже когда сам VCF
    # уже лежит на диске (из общего кэша сырых хромосом) и сеть тут
    # вообще не участвует. Сортировка по позиции — единственное, что
    # меняется, сам список позиций тот же.
    chrom_positions = sorted(
        (p for c, p in chip_positions if c == str(chrom))
    )
    # ⚠ ИСПРАВЛЕНИЕ (найдено по факту: 0 записей на ВСЕХ 22 хромосомах,
    # одинаково с VPN и без него — то есть проблема не сетевая, а в
    # имени контига). GENOME_BUILD_CHROM_PREFIX["grch38"] = "chr"
    # предполагал, что GRCh38-релиз 1000 Genomes называет CHROM как
    # "chr1" (по аналогии с самим референсным .fa) — но реальный
    # заголовок скачанного файла показывает "##contig=<ID=1>": этот
    # конкретный датасет (20190312_biallelic_SNV_and_INDEL) хранит CHROM
    # БЕЗ префикса, как классический phase3. Из-за этого pos_chr{N}.txt
    # искал "chr21", когда в файле реально "21" — bcftools view -R не
    # находил НИ ОДНОГО совпадения, независимо от качества сети.
    #
    # По аналогии с уже проверенным фиксом той же природы в
    # template/assembler.py::load_imputed_genotypes() (см. его
    # комментарии про "обе формы хромосомы") — пишем ОБЕ формы имени
    # хромосомы для каждой позиции. bcftools -R не требует, чтобы каждая
    # строка регион-файла совпала с чем-то в индексе — несовпавшие формы
    # просто безвредно игнорируются, а совпавшая форма находится всегда,
    # независимо от того, как реально называет свои контиги конкретный
    # датасет доноров.
    pos_file = output_dir / f"pos_chr{chrom}.txt"
    with pos_file.open('w', newline='\n') as f:
        for p in chrom_positions:
            f.write(f"{chrom}\t{p}\n")
            f.write(f"chr{chrom}\t{p}\n")
    # общем кэше сырых хромосом фильтруем в ДВА шага — сначала (один раз
    # на genome_build, дальше из кэша) убираем 2504-20=2484 ненужных
    # образца, затем фильтруем по позициям УЖЕ МАЛЕНЬКИЙ файл (~0.8% от
    # исходного размера). Без raw_cache_dir двухшаговая схема не имеет
    # смысла — кэшировать нечего, а дополнительный проход по диску только
    # замедлил бы разовую фильтрацию, поэтому остаётся старая однокомандная
    # схема (-S и -R вместе, один проход по входному файлу).
    subset_source = tmp_vcf
    subset_has_samples_filtered = False
    if raw_cache_dir is not None:
        try:
            subset_source = _ensure_eur_subset(
                tmp_vcf, eur_file, chrom, htslib, output_dir,
                bcftools_threads=bcftools_threads, raw_cache_dir=raw_cache_dir,
                eur_sample_count=eur_sample_count,
            )
            subset_has_samples_filtered = True
        except RuntimeError as e:
            print(f"  ⚠ Подвыборка {eur_sample_count} образцов не удалась ({e}) — фильтрую напрямую из полной хромосомы.")
            subset_source = tmp_vcf
            subset_has_samples_filtered = False

    if subset_has_samples_filtered:
        # Образцы и биаллельные SNP уже отфильтрованы на этапе подвыборки
        # (см. _ensure_eur_subset) — здесь остаётся только -R по позициям.
        cmd = [
            htslib.bcftools_path, "view",
            "-R", str(pos_file),
            "--threads", str(max(1, bcftools_threads)),
            str(subset_source),
            "-Oz", "-o", str(out_file),
        ]
    else:
        cmd = [
            htslib.bcftools_path, "view",
            "-S", str(eur_file), "--force-samples",
            "-R", str(pos_file),
            "-m2", "-M2", "-v", "snps",
            "--threads", str(max(1, bcftools_threads)),
            str(subset_source),
            "-Oz", "-o", str(out_file),
        ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"   ✗ bcftools ошибка:\n     {result.stderr}")
            out_file.unlink(missing_ok=True)
            return False

        if htslib.has_tabix:
            subprocess.run([htslib.tabix_path, "-p", "vcf", str(out_file)], check=True, capture_output=True)
            print(f"  ✓ chr{chrom}: отфильтрован и проиндексирован")
        else:
            print(f"  ✓ chr{chrom}: отфильтрован (без индексации)")

        # Промт "проверка 'донор не пустой' после скачивания/фильтрации":
        # returncode==0 и ненулевой размер файла подтверждают только, что
        # bcftools записал корректный VCF-заголовок — этого недостаточно,
        # если тело файла оборвалось на середине (нестабильное сетевое
        # соединение при скачивании исходной хромосомы, VPN/прокси и
        # т.п. — см. докстринг _count_vcf_records()). Явно считаем число
        # записей и отбраковываем хромосому как неудавшуюся, если их 0.
        record_count = _count_vcf_records(htslib.bcftools_path, out_file)
        if record_count == 0:
            print(
                f"  ✗ chr{chrom}: донор получился пустым (0 записей после "
                f"фильтрации) — {_EMPTY_DONOR_HINT}"
            )
            out_file.unlink(missing_ok=True)
            out_file.with_suffix(out_file.suffix + ".tbi").unlink(missing_ok=True)
            return False
        elif record_count > 0:
            print(f"  ℹ chr{chrom}: проверено — {record_count} записей в доноре")
        # record_count == -1 — сама проверка не удалась (bcftools упал на
        # чтении только что построенного файла) — не трактуем это как
        # признак пустого файла, продолжаем как раньше.

        # Промт "выровнять CHROM доноров с sample.vcf.gz для merge" — см.
        # докстринг _ensure_donor_chrom_prefix(). Делается ПОСЛЕ проверки
        # "донор не пустой" (нет смысла переименовывать заведомо пустой
        # файл) и ДО чистки временных файлов, чтобы при неудаче
        # переименования cleanup всё равно отработал одинаково.
        if not _ensure_donor_chrom_prefix(htslib, out_file, chrom, genome_build):
            out_file.unlink(missing_ok=True)
            out_file.with_suffix(out_file.suffix + ".tbi").unlink(missing_ok=True)
            return False

        eur_file.unlink(missing_ok=True)
        pos_file.unlink(missing_ok=True)
        tmp_vcf.unlink(missing_ok=True)
        tmp_vcf.with_suffix(tmp_vcf.suffix + ".tbi").unlink(missing_ok=True)
        # Локальную временную копию подвыборки (не путать с записью в
        # raw_cache_dir, куда она уже сохранена хардлинком/копией внутри
        # _ensure_eur_subset) можно удалить — если она была взята прямо
        # из кэша, subset_source указывает НА файл в кэше, и удалять его
        # нельзя (это разные объекты только когда файл только что
        # построен в output_dir и ещё не был захардлинкан "на месте").
        if subset_has_samples_filtered and subset_source.parent == output_dir:
            subset_source.unlink(missing_ok=True)
            subset_source.with_suffix(subset_source.suffix + ".tbi").unlink(missing_ok=True)
        return True
    except Exception as e:
        print(f"   Ошибка обработки chr{chrom}: {e}")
        out_file.unlink(missing_ok=True)
        return False


# ---------------------------------------------------------------------------
# Часть 1.1: удалённая фильтрация без полного скачивания хромосомы
# ---------------------------------------------------------------------------
def _probe_bcftools_remote_support(
    htslib: HtslibTools,
    cancel_check: Optional[Callable[[], bool]] = None,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> bool:
    """
    genome_build (промт "TopMed/HRC", п.3): проба выполняется на зеркале/
    шаблоне ИМЕННО той сборки, для которой реально будет запущено
    скачивание доноров — иначе для panel="topmed" эта функция продолжала
    бы тестировать GRCh37-URL (REMOTE_PROBE_CHROM=21 у GRCh37) и могла бы
    дать ложноположительный/ложноотрицательный результат относительно
    единственного подтверждённого GRCh38-зеркала (см. GRCH38_MIRRORS).

    Проверяет, умеет ли установленный bcftools/htslib читать удалённые
    файлы по HTTP(S) через Range-запросы (собран с libcurl) — это нужно
    process_chromosome_remote(), которая фильтрует доноров прямо по URL
    без скачивания хромосомы целиком (см. докстринг модуля, "Часть 1.1").

    Пробный запрос — `bcftools view -h` (только заголовок VCF, без тела)
    на самой маленькой хромосоме (chr21) с первого зеркала. Это требует
    лишь нескольких Range-запросов к .tbi-индексу и началу файла, а не
    скачивания тела VCF — секунды даже на медленном канале.

    Возвращает False (никогда не бросает исключение, кроме отмены) при
    ЛЮБОЙ проблеме — отсутствии libcurl, сетевой ошибке, таймауте,
    отсутствии bcftools: вызывающий код просто откатывается на полное
    скачивание для всех хромосом, как и раньше в v11.
    """
    if not htslib.has_bcftools:
        return False
    if cancel_check and cancel_check():
        return False

    # v13-фикс: раньше проба использовала ТОЛЬКО VCF_SUFFIX_CANDIDATES[0]
    # ("v5b") — но живая проверка на реальном S3-зеркале 1000genomes
    # показала, что для chr21 там существует только суффикс "v5a" (см.
    # диагностику в докстринге модуля выше). Проба с одним лишь "v5b"
    # получала бы "No such file or directory" и ЛОЖНО делала вывод, что
    # bcftools не умеет читать удалённые файлы вовсе — хотя реальная
    # причина в отсутствующем конкретном файле, а не в отсутствии
    # libcurl. Перебираем оба кандидата суффикса на первом зеркале, как
    # это и так делает download_chromosome_vcf()/process_chromosome_remote()
    # для реальных хромосом — иначе проба и боевой путь расходятся в
    # логике и дают противоречивые результаты.
    last_stderr_tail = ""
    build_mirrors = _mirrors_for_build(genome_build)
    build_template = _vcf_template_for_build(genome_build)
    build_suffix_candidates = _vcf_suffix_candidates_for_build(genome_build)
    for suffix in build_suffix_candidates:
        if cancel_check and cancel_check():
            return False
        probe_url = build_mirrors[0] + build_template.format(
            chrom=REMOTE_PROBE_CHROM, suffix=suffix,
        )
        cmd = [htslib.bcftools_path, "view", "-h", probe_url]
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=REMOTE_PROBE_TIMEOUT,
            )
        except subprocess.TimeoutExpired:
            print(f"ℹ Проба удалённого доступа bcftools (суффикс {suffix}) не "
                  f"уложилась в {REMOTE_PROBE_TIMEOUT}с — пробую другой суффикс "
                  f"(если есть) или скачиваю хромосомы целиком.")
            continue
        except (FileNotFoundError, OSError) as e:
            print(f"ℹ Проба удалённого доступа bcftools не удалась ({e}) — "
                  f"буду скачивать хромосомы целиком.")
            return False

        if result.returncode == 0 and result.stdout.strip():
            print(
                f"✓ bcftools поддерживает удалённое чтение по HTTP(S) (libcurl, "
                f"суффикс {suffix}) — буду фильтровать доноров без полного "
                f"скачивания хромосом, где это возможно."
            )
            return True

        last_stderr_tail = (result.stderr or "").strip()[-300:]

    stderr_tail = last_stderr_tail
    print(
        f"ℹ bcftools не смог прочитать удалённый файл напрямую "
        f"(вероятно, сборка без libcurl): {stderr_tail or 'нет вывода stderr'} — "
        f"буду скачивать хромосомы целиком."
    )
    return False


# ---------------------------------------------------------------------------
# v13, Часть 1: диагностика РЕАЛЬНОЙ удалённой фильтрации (не только -h)
# ---------------------------------------------------------------------------
def diagnose_remote_filter(
    htslib: HtslibTools,
    output_dir: Path,
    test_chrom: int = REMOTE_PROBE_CHROM,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> dict:
    """
    genome_build (промт "TopMed/HRC", п.3, докрутка диагностики за компанию
    с боевым путём): выбирает зеркала/шаблон имени файла через
    _mirrors_for_build()/_vcf_template_for_build()/
    _vcf_suffix_candidates_for_build() — без этого параметра диагностика
    для panel="topmed" тестировала бы GRCh37-зеркала/файлы, что не имеет
    отношения к тому, что реально будет использовано при скачивании
    доноров под TopMed, и могла бы давать вводящий в заблуждение результат.

    ⚠ diag_positions ниже подобраны как валидные координаты в пределах
    контига chr21 GRCh37/hg19 (см. комментарий у diag_positions) — для
    genome_build="grch38" те же числовые координаты почти наверняка НЕ
    совпадут с реальными вариантами GRCh38-релиза (сборки используют
    разную нумерацию для одного и того же биологического локуса). Это не
    ломает диагностику МЕХАНИЗМА удалённого доступа (bcftools либо
    успешно подключится и отфильтрует 0 записей, либо провалится с той
    же сетевой/Range-ошибкой независимо от того, есть ли варианты на этих
    координатах) — но report["records"] == 0 для GRCh38 может означать
    просто "нет вариантов на этих координатах", а не "фильтрация не
    работает". Пользователю, ориентирующемуся на report["ok"], это не
    мешает (ok=False только при реальной ошибке subprocess/сети), но
    report["records"] для GRCh38 менее показателен, чем для GRCh37.

    Проверяет, что удалённая фильтрация (та же команда, что использует
    process_chromosome_remote() на каждой из 22 хромосом в бою) реально
    отдаёт непустой валидный результат — а не только то, что
    _probe_bcftools_remote_support() умеет открыть `-h` (заголовок без
    тела). Разница существенна: короткий запрос заголовка требует лишь
    нескольких Range-запросов к .tbi/началу файла, тогда как реальная
    фильтрация по `-S`/`-R` требует читать множество Range-блоков по
    всему телу файла — там могут вылезти проблемы, которых нет при
    простом `-h` (обрыв соединения на середине, битый удалённый .tbi,
    прокси/файрвол, обрезающий длинные keep-alive запросы, и т.п.).

    Строит минимальные eur20.txt (2-3 образца, не все 20 — тест не
    должен зависеть от полного списка) и pos.txt (несколько РЕАЛЬНЫХ
    позиций test_chrom — намеренно не заимствованы из чипа пользователя,
    чтобы диагностику можно было запускать до того, как позиции чипа
    вообще посчитаны/известны), выполняет РОВНО ТУ ЖЕ команду, что
    process_chromosome_remote(), во ВРЕМЕННЫЙ файл (никогда не пишет в
    kgp_sub_{chrom}.vcf.gz — чтобы не перепутать диагностику с реальным
    кэшем доноров и не блокировать последующий настоящий прогон), и
    считает записи через `bcftools view -H | подсчёт строк`.

    Возвращает словарь-отчёт:
      {
        "ok": bool,               — успех/неудача теста
        "records": int,           — сколько строк вариантов вернулось
        "stderr_tail": str,       — хвост stderr bcftools (для диагностики)
        "duration_sec": float,    — время выполнения теста
        "url": str,               — какой URL реально тестировался
      }

    НИКОГДА не бросает исключение — это диагностика, а не критический
    путь: любая ошибка (нет bcftools, нет сети, битый .tbi) отражается в
    отчёте ("ok": False, "stderr_tail": ...), а не прерывает вызывающий
    код.
    """
    import time as _time

    report = {
        "ok": False, "records": 0, "stderr_tail": "",
        "duration_sec": 0.0, "url": "",
    }

    if not htslib.has_bcftools:
        report["stderr_tail"] = "bcftools не найден (проверьте --bin-dir)"
        return report

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Минимальный список образцов — тест механизма фильтрации, а не
    # содержимого; полные 20 образцов не нужны и лишь замедлили бы тест.
    tmp_eur = output_dir / "_diag_eur.txt"
    tmp_pos = output_dir / "_diag_pos.txt"
    tmp_out = output_dir / f"_diag_chr{test_chrom}.vcf.gz"
    tmp_out_tbi = tmp_out.with_suffix(tmp_out.suffix + ".tbi")

    try:
        try:
            eur_samples = create_eur20_list(output_dir, bin_dir=htslib.bin_dir)
        except RuntimeError as e:
            report["stderr_tail"] = f"Не удалось получить список образцов: {e}"
            return report
        few_samples = eur_samples[:3] if len(eur_samples) >= 3 else eur_samples
        if not few_samples:
            report["stderr_tail"] = "Список европейских образцов пуст"
            return report
        with tmp_eur.open("w") as f:
            for s in few_samples:
                f.write(s + "\n")

        # Несколько заведомо валидных позиций в пределах контига chr21
        # (b37/hg19, длина ~48.1 Мб — см. ##contig в заголовке VCF) —
        # не обязаны совпадать с реальным чипом пользователя, здесь
        # проверяется сам МЕХАНИЗМ удалённой фильтрации, а не то, что
        # чип действительно содержит эти rsid. Для genome_build="grch38"
        # эти координаты не гарантированно валидны (см. докстринг выше) —
        # диагностика механизма всё равно корректна, показательность
        # report["records"]==0 менее надёжна.
        diag_positions = [9411239, 9411264, 9411267, 9411302, 9411313]
        # См. комментарий в process_chromosome() — реальный GRCh38-релиз
        # 1000 Genomes хранит CHROM без префикса "chr", вопреки
        # GENOME_BUILD_CHROM_PREFIX. Пишем обе формы, чтобы диагностика
        # не давала ложный "0 записей" из-за того же несовпадения имени
        # контига, которое ловилось в боевом пути.
        with tmp_pos.open("w") as f:
            for pos in diag_positions:
                f.write(f"{test_chrom}\t{pos}\n")
                f.write(f"chr{test_chrom}\t{pos}\n")

        started = _time.monotonic()
        found_ok = False
        last_stderr = ""
        used_url = ""

        build_mirrors = _mirrors_for_build(genome_build)
        build_template = _vcf_template_for_build(genome_build)
        build_suffix_candidates = _vcf_suffix_candidates_for_build(genome_build)

        for mirror in build_mirrors:
            mirror_name = mirror.split("/")[2]
            for suffix in build_suffix_candidates:
                url = mirror + build_template.format(chrom=test_chrom, suffix=suffix)
                cmd = [
                    htslib.bcftools_path, "view",
                    "-S", str(tmp_eur), "--force-samples",
                    "-R", str(tmp_pos),
                    url, "-Oz", "-o", str(tmp_out),
                ]
                try:
                    result = subprocess.run(
                        cmd, capture_output=True, text=True,
                        timeout=REMOTE_CHROM_TIMEOUT,
                    )
                except subprocess.TimeoutExpired:
                    last_stderr = f"Таймаут на {mirror_name} [{suffix}]"
                    continue

                if result.returncode == 0 and tmp_out.exists() and tmp_out.stat().st_size > 0:
                    found_ok = True
                    used_url = url
                    break
                last_stderr = (result.stderr or "").strip()[-500:] or f"пусто ({mirror_name}/{suffix})"
            if found_ok:
                break

        duration = _time.monotonic() - started
        report["duration_sec"] = round(duration, 2)
        report["url"] = used_url
        # Фикс: last_stderr — это хвост stderr от ПОСЛЕДНЕЙ НЕУДАЧНОЙ
        # попытки (например, v5b, которого нет на S3, — обычная работа
        # перебора суффиксов). Если в итоге found_ok=True благодаря
        # следующему суффиксу/зеркалу (например, v5a), не нужно показывать
        # пользователю этот "хвост" как будто это ошибка текущего успешного
        # результата — иначе отчёт с "Успех: ДА" вводит в заблуждение
        # строкой [E::hts_open_format]... из уже отброшенной попытки.
        report["stderr_tail"] = "" if found_ok else last_stderr

        if not found_ok:
            return report

        count_cmd = [htslib.bcftools_path, "view", "-H", str(tmp_out)]
        count_res = subprocess.run(count_cmd, capture_output=True, text=True)
        if count_res.returncode != 0:
            report["stderr_tail"] = (count_res.stderr or "").strip()[-500:]
            return report

        records = len([l for l in count_res.stdout.splitlines() if l.strip()])
        report["records"] = records
        report["ok"] = records > 0
        if records == 0:
            report["stderr_tail"] = (
                "Команда выполнилась успешно, но вернула 0 записей — "
                "проверьте позиции/образцы (для реальной диагностики это "
                "не проблема самого механизма, но обратите внимание)."
            )
        return report
    except Exception as e:  # noqa: BLE001 — диагностика не должна падать
        report["stderr_tail"] = f"Неожиданная ошибка диагностики: {e}"
        return report
    finally:
        tmp_eur.unlink(missing_ok=True)
        tmp_pos.unlink(missing_ok=True)
        tmp_out.unlink(missing_ok=True)
        tmp_out_tbi.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# v14: классификация ошибок удалённого Range-доступа (промт "точечный патч
# remote-фильтрации для крупных хромосом").
#
# Подтверждено реальным прогоном на 758 990 позициях чипа (chr1/chr2/chr3):
# ошибка "файла с этим суффиксом нет" (404 / "No such file or directory")
# принципиально отличается от ошибки "зеркало не может корректно отдать
# Range-блок для ЭТОГО файла вообще" — вторая не зависит от суффикса и
# повторный перебор суффиксов на том же зеркале почти наверняка провалится
# так же. Наблюдавшиеся в логе маркеры второго типа:
#   [E::bgzf_read_block] Failed to read BGZF header at offset ...
#   [E::hts_itr_next] Failed to seek to offset ...: Invalid seek
# (плюс общий "Error: BCF read error", который bcftools печатает вслед за
# обоими случаями выше — это НЕ самостоятельный маркер, а следствие).
_TRANSIENT_RANGE_ERROR_MARKERS = (
    "bgzf_read_block",
    "Failed to read BGZF header",
    "hts_itr_next",
    "Invalid seek",
)


def _is_transient_range_error(stderr: str) -> bool:
    """
    True, если stderr указывает на системную ошибку Range-доступа зеркала
    (битый/оборванный блок, рассинхрон индекса и тела файла) — а не на
    "файла с этим суффиксом нет" (404/No such file or directory).

    Используется в process_chromosome_remote(), чтобы при таком типе
    ошибки сразу переходить к следующему ЗЕРКАЛУ, не тратя время на
    оставшиеся суффиксы текущего — они статистически проваливаются так
    же (см. докстринг модуля, реальный лог прогона v13).
    """
    return any(marker in stderr for marker in _TRANSIENT_RANGE_ERROR_MARKERS)


def process_chromosome_remote(
    chrom: int,
    eur_samples: Set[str],
    chip_positions: Set[Tuple[str, int]],
    output_dir: Path,
    htslib: HtslibTools,
    working_suffix_by_mirror: dict[str, str],
    cancel_check: Optional[Callable[[], bool]] = None,
    bcftools_threads: int = 1,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> bool:
    """
    genome_build (промт "TopMed/HRC", п.3): та же роль, что и в
    process_chromosome() — выбор зеркал/шаблона/суффиксов через
    _mirrors_for_build()/_vcf_template_for_build()/
    _vcf_suffix_candidates_for_build() и правильный префикс CHROM в
    pos_chr{chrom}.txt через _chrom_label_for_build().

    Фильтрует доноров для ОДНОЙ хромосомы напрямую по удалённому URL —
    `bcftools view -S eur20.txt -R pos.txt "<URL>" -Oz -o kgp_sub_N.vcf.gz`
    — без промежуточного скачивания хромосомы целиком (см. "Часть 1.1").
    htslib сам делает HTTP Range-запросы через libcurl за нужными блоками,
    подтягивая удалённый .tbi по мере надобности.

    Возвращает False при ЛЮБОЙ неудаче (сеть, битое зеркало, таймаут,
    ошибка bcftools) — НЕ бросает исключение (кроме DownloadCancelled),
    чтобы process_chromosome_auto() мог тихо откатиться на
    process_chromosome() (полное скачивание) именно для этой хромосомы,
    не роняя обработку остальных.
    """
    out_file = output_dir / f"kgp_sub_{chrom}.vcf.gz"
    if out_file.exists() and out_file.stat().st_size > 0:
        print(f"  ✓ chr{chrom} уже готов")
        return True

    if cancel_check and cancel_check():
        raise DownloadCancelled("Скачивание отменено пользователем")

    eur_file = output_dir / f"eur{len(eur_samples)}_chr{chrom}.txt"
    with eur_file.open('w', newline='\n') as f:
        for s in eur_samples:
            f.write(s + '\n')

    # v18: та же сортировка позиций, что и в process_chromosome() — см.
    # докстринг там, почему это важно для bcftools view -R.
    chrom_positions = sorted(
        (p for c, p in chip_positions if c == str(chrom))
    )
    # См. подробный комментарий в process_chromosome() про причину этого
    # фикса — реальный GRCh38-релиз 1000 Genomes хранит CHROM без
    # префикса "chr" ("##contig=<ID=1>"), а не с ним, вопреки
    # GENOME_BUILD_CHROM_PREFIX. Пишем обе формы, чтобы не зависеть от
    # соглашения конкретного датасета.
    pos_file = output_dir / f"pos_chr{chrom}.txt"
    with pos_file.open('w', newline='\n') as f:
        for p in chrom_positions:
            f.write(f"{chrom}\t{p}\n")
            f.write(f"chr{chrom}\t{p}\n")

    ok = False
    # v14, Шаг 3 промта: суммарное время, реально потраченное на remote-
    # попытки для этой хромосомы (по всем зеркалам/суффиксам) — печатается
    # по завершении независимо от исхода. Нужно для будущей калибровки
    # порога remote_skip_large_chroms (Шаг 2): если для крупных хромосом
    # это время стабильно мало (как показал реальный прогон — единицы
    # секунд на bgzf_read_block/Invalid seek, а не таймауты), то основные
    # накладные расходы дают не сами remote-попытки, а последующий откат
    # на полное скачивание, а не перебор как таковой.
    remote_attempt_started = time.monotonic()
    build_mirrors = _mirrors_for_build(genome_build)
    build_template = _vcf_template_for_build(genome_build)
    build_suffix_candidates = _vcf_suffix_candidates_for_build(genome_build)
    try:
        for mirror in build_mirrors:
            if cancel_check and cancel_check():
                raise DownloadCancelled("Скачивание отменено пользователем")

            mirror_name = mirror.split('/')[2]
            known_suffix = working_suffix_by_mirror.get(mirror)
            suffixes_to_try = [known_suffix] if known_suffix else build_suffix_candidates

            for suffix in suffixes_to_try:
                url = mirror + build_template.format(chrom=chrom, suffix=suffix)
                print(f"\n  --- chr{chrom} (удалённо, {mirror_name}, суффикс {suffix}) ---")

                cmd = [
                    htslib.bcftools_path, "view",
                    "-S", str(eur_file), "--force-samples",
                    "-R", str(pos_file),
                    "-m2", "-M2", "-v", "snps",
                    "--threads", str(max(1, bcftools_threads)),
                    url,
                    "-Oz", "-o", str(out_file),
                ]
                try:
                    result = subprocess.run(
                        cmd, capture_output=True, text=True, timeout=REMOTE_CHROM_TIMEOUT,
                    )
                except subprocess.TimeoutExpired:
                    print(f"  ✗ Таймаут удалённого чтения chr{chrom} с {mirror_name}")
                    out_file.unlink(missing_ok=True)
                    continue

                if result.returncode == 0 and out_file.exists() and out_file.stat().st_size > 0:
                    if htslib.has_tabix:
                        tbi_res = subprocess.run(
                            [htslib.tabix_path, "-p", "vcf", str(out_file)],
                            capture_output=True,
                        )
                        if tbi_res.returncode != 0:
                            print(f"  ✗ Не удалось построить индекс для chr{chrom} — считаю неудачей")
                            out_file.unlink(missing_ok=True)
                            continue
                        print(f"  ✓ chr{chrom}: отфильтровано удалённо и проиндексировано "
                              f"(без полного скачивания)")
                    else:
                        print(f"  ✓ chr{chrom}: отфильтровано удалённо, без индексации (нет tabix)")

                    # Промт "проверка 'донор не пустой' после скачивания/
                    # фильтрации": returncode==0 у удалённой Range-
                    # фильтрации подтверждает только, что bcftools дошёл до
                    # конца команды — сам подтверждённый реальным прогоном
                    # баг именно в том, что при обрыве HTTPS-соединения на
                    # теле файла (нестабильная сеть/VPN/прокси) результат
                    # получается с корректным заголовком, но без единой
                    # записи, а returncode всё равно 0. Явно проверяем.
                    record_count = _count_vcf_records(htslib.bcftools_path, out_file)
                    if record_count == 0:
                        print(
                            f"  ✗ chr{chrom}: удалённая фильтрация вернула "
                            f"пустой результат (0 записей) с {mirror_name} "
                            f"[{suffix}] — {_EMPTY_DONOR_HINT} Пробую "
                            f"следующий вариант."
                        )
                        out_file.unlink(missing_ok=True)
                        out_file.with_suffix(out_file.suffix + ".tbi").unlink(missing_ok=True)
                        continue
                    elif record_count > 0:
                        print(f"  ℹ chr{chrom}: проверено — {record_count} записей")
                    # record_count == -1 — сама проверка не удалась, не
                    # трактуем как признак пустого файла (см. докстринг
                    # _count_vcf_records()).

                    # Промт "выровнять CHROM доноров с sample.vcf.gz для
                    # merge" — см. докстринг _ensure_donor_chrom_prefix().
                    if not _ensure_donor_chrom_prefix(htslib, out_file, chrom, genome_build):
                        out_file.unlink(missing_ok=True)
                        out_file.with_suffix(out_file.suffix + ".tbi").unlink(missing_ok=True)
                        continue

                    working_suffix_by_mirror[mirror] = suffix
                    ok = True
                    break
                else:
                    stderr_text = result.stderr or ""
                    stderr_tail = stderr_text.strip()[-300:]
                    print(f"  ✗ Удалённое чтение chr{chrom} с {mirror_name} [{suffix}] "
                          f"не удалось: {stderr_tail or 'нет вывода stderr'}")
                    out_file.unlink(missing_ok=True)

                    # v14, Шаг 1 промта: "файла с этим суффиксом нет"
                    # (404/No such file or directory) — нормальная ситуация,
                    # имеет смысл пробовать следующий суффикс на ЭТОМ ЖЕ
                    # зеркале, как и раньше. Но если ошибка указывает на
                    # системную проблему Range-доступа зеркала к ЭТОМУ
                    # файлу (bgzf_read_block/Invalid seek и т.п.) —
                    # подтверждено реальным прогоном, что второй суффикс
                    # на том же зеркале проваливается так же — не тратим
                    # время, сразу переходим к следующему зеркалу.
                    if _is_transient_range_error(stderr_text):
                        print(
                            f"  ℹ chr{chrom}: похоже на системную ошибку "
                            f"Range-доступа зеркала {mirror_name} (не "
                            f"'файла нет') — пропускаю оставшиеся суффиксы "
                            f"этого зеркала и перехожу к следующему."
                        )
                        break

            if ok:
                break
    finally:
        eur_file.unlink(missing_ok=True)
        pos_file.unlink(missing_ok=True)
        remote_elapsed = time.monotonic() - remote_attempt_started
        print(
            f"  ℹ chr{chrom}: на remote-попытки потрачено {remote_elapsed:.1f}с "
            f"(исход: {'успех' if ok else 'откат на полное скачивание'})"
        )

    return ok


def process_chromosome_auto(
    chrom: int,
    eur_samples: Set[str],
    chip_positions: Set[Tuple[str, int]],
    output_dir: Path,
    htslib: HtslibTools,
    working_suffix_by_mirror: dict[str, str],
    remote_capable: bool,
    cancel_check: Optional[Callable[[], bool]] = None,
    skip_remote: bool = False,
    raw_cache_dir: Optional[Path] = None,
    bcftools_threads: int = 1,
    genome_build: str = DEFAULT_GENOME_BUILD,
) -> bool:
    """
    genome_build (промт "TopMed/HRC", п.3): прокидывается и в
    process_chromosome_remote(), и в process_chromosome() (оба пути должны
    знать, с какой сборкой работают), и в _raw_cache_has_chrom() ниже —
    иначе проверка "уже есть в общем кэше сырых хромосом" всегда смотрела
    бы в GRCh37-часть кэша независимо от выбранной панели.

    raw_cache_dir (v16) — прокидывается только в process_chromosome()
    (полное скачивание); process_chromosome_remote() его не использует —
    общий кэш сырых хромосом актуален только для пути полного скачивания
    (см. докстринг файла).

    Единая точка выбора способа обработки одной хромосомы (Часть 1,
    цепочка приоритетов): удалённая фильтрация (1.1, быстрее и дешевле
    по трафику), если remote_capable — иначе/при неудаче сразу полное
    скачивание (process_chromosome(), которое само использует
    aria2c/curl/urllib внутри download_file(), см. 1.2/1.4).

    working_suffix_by_mirror — ОДИН общий словарь для обоих путей: URL
    (зеркало + суффикс) идентичен что для удалённой фильтрации, что для
    полного скачивания, так что найденный рабочий суффикс полезен в
    обоих направлениях и не нужно поддерживать раздельные кэши.

    skip_remote (v14, Шаг 2 промта "точечный патч remote-фильтрации для
        крупных хромосом"): если True, remote-путь для ЭТОЙ хромосомы не
        пробуется вовсе — сразу process_chromosome() (полное скачивание).
        Подтверждено реальным прогоном на 758 990 позициях чипа: для
        крупных хромосом (chr1/chr2/chr3 — счёт на ГБ) remote-путь
        систематически проваливается на всех зеркалах/суффиксах с
        ошибками bgzf_read_block/Invalid seek (см. _is_transient_range_error)
        и лишь откладывает неизбежный откат на полное скачивание. Каждая
        такая попытка сама по себе быстрая (секунды, не таймаут), но при
        параллельной обработке (max_parallel_chromosomes > 1) несколько
        таких "быстрых провалов" одновременно всё равно платят латентностью
        сети и создают лишнюю нагрузку — решение о конкретном списке/
        пороге хромосом принимает вызывающий код (download_donors_for_chip,
        параметр remote_skip_large_chroms), эта функция только его
        применяет.
    """
    if remote_capable and skip_remote:
        print(
            f"  ℹ chr{chrom}: remote-путь пропущен (в списке "
            f"remote_skip_large_chroms) — сразу полное скачивание."
        )

    # v17: если полная нефильтрованная хромосома уже лежит в общем кэше
    # сырых хромосом (скачана ранее для другого источника/чипа) — remote-
    # путь для неё смысла не имеет: process_chromosome() и так возьмёт её
    # оттуда почти мгновенно (см. download_chromosome_vcf() ->
    # _load_from_raw_cache()), тогда как remote-попытка может впустую
    # потратить до REMOTE_CHROM_TIMEOUT секунд на каждое зеркало.
    if remote_capable and not skip_remote and raw_cache_dir is not None:
        if _raw_cache_has_chrom(raw_cache_dir, chrom, genome_build=genome_build):
            print(
                f"  ℹ chr{chrom}: уже есть в общем кэше сырых хромосом "
                f"({raw_cache_dir}) — remote-путь пропущен, использую "
                f"локальный файл напрямую."
            )
            skip_remote = True

    if remote_capable and not skip_remote:
        try:
            if process_chromosome_remote(
                chrom, eur_samples, chip_positions, output_dir,
                htslib, working_suffix_by_mirror, cancel_check,
                bcftools_threads=bcftools_threads, genome_build=genome_build,
            ):
                return True
            print(
                f"  ℹ chr{chrom}: удалённая фильтрация не удалась ни на одном "
                f"зеркале — откатываюсь на полное скачивание."
            )
        except DownloadCancelled:
            raise
        except Exception as e:  # noqa: BLE001 — сбой удалённого пути не должен ронять хромосому целиком
            print(
                f"  ⚠ chr{chrom}: ошибка удалённой фильтрации ({e}) — "
                f"откатываюсь на полное скачивание."
            )

    return process_chromosome(
        chrom, eur_samples, chip_positions, output_dir,
        htslib, working_suffix_by_mirror, cancel_check,
        raw_cache_dir=raw_cache_dir,
        bcftools_threads=bcftools_threads,
        genome_build=genome_build,
    )


# ---------------------------------------------------------------------------
# Инвалидация устаревшего кэша (Задача 1, п.7)
# ---------------------------------------------------------------------------
def _invalidate_stale_donor_cache(
    output_dir: Path,
    expected_signature: str,
    progress_cb: Optional[Callable[[float, str], None]] = None,
) -> bool:
    """
    Если в output_dir уже лежит chip_signature.txt от ДРУГОГО чипа —
    удаляет все kgp_sub_*.vcf.gz(.tbi) и служебные файлы, привязанные к
    сигнатуре, чтобы process_chromosome() не принял их за "уже готовые"
    для нового чипа (см. докстринг файла). Возвращает True, если
    инвалидация действительно произошла.

    Если chip_signature.txt отсутствует — ничего не делает: либо это
    первый запуск, либо на диске уже лежат недостроенные/частичные файлы
    без сигнатуры, для которых process_chromosome() и так докачает
    недостающее по хромосомам (это нормальный сценарий докачки после
    сетевой ошибки, Задача 1, п.6).
    """
    sig_file = output_dir / "chip_signature.txt"
    if not sig_file.exists():
        return False

    cached = sig_file.read_text(encoding="utf-8").strip()
    if cached == expected_signature:
        return False

    msg = (
        f"⚠ Кэш доноров в {output_dir} принадлежит другому чипу "
        f"(на диске signature={cached}, ожидается {expected_signature}) — "
        f"удаляю устаревших доноров перед перекачкой под новый чип."
    )
    print(msg)
    if progress_cb:
        progress_cb(0.0, "Обнаружен кэш другого чипа — очистка...")

    for chrom in range(1, 23):
        (output_dir / f"kgp_sub_{chrom}.vcf.gz").unlink(missing_ok=True)
        (output_dir / f"kgp_sub_{chrom}.vcf.gz.tbi").unlink(missing_ok=True)
    for name in _SIGNATURE_SCOPED_STATIC_FILES:
        (output_dir / name).unlink(missing_ok=True)
    # Имя файла со списком EUR-образцов зависит от eur_sample_count
    # (eur20.txt/eur120.txt/...) — фиксированный список имён его не
    # покрывает, чистим через glob.
    for stale_eur_file in output_dir.glob("eur*.txt"):
        stale_eur_file.unlink(missing_ok=True)
    return True


# ---------------------------------------------------------------------------
# Задача 1: единая вызываемая точка входа (без argparse/sys.exit)
# ---------------------------------------------------------------------------
def download_donors_for_chip(
    positions_json: Path,
    source: str,
    output_dir: Path,
    htslib: HtslibTools,
    progress_cb: Optional[Callable[[float, str], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    remote_filter: bool = True,
    max_parallel_chromosomes: int = DEFAULT_PARALLEL_CHROMOSOMES,
    remote_skip_large_chroms: Optional[Set[int]] = DEFAULT_REMOTE_SKIP_LARGE_CHROMS,
    raw_cache_dir: Optional[Path] = None,
    bcftools_threads: Optional[int] = None,
    genome_build: str = DEFAULT_GENOME_BUILD,
    eur_sample_count: Optional[int] = DEFAULT_EUR_SAMPLE_COUNT,
) -> list[Path]:
    """
    eur_sample_count (промт "Monomorphic sites / настраиваемое количество
        EUR-доноров"): сколько EUR-образцов 1000 Genomes использовать как
        доноров. None (по умолчанию, EUR_SAMPLE_COUNT_ALL) — взять всю
        доступную EUR-подвыборку панели (не только фиксированные 20) —
        уменьшает долю Monomorphic sites на QC MIS ценой большего трафика
        и времени обработки каждой хромосомы. Можно передать явное число
        (например, 20 — прежнее поведение, 100 — компромисс по
        трафику/скорости). См. create_eur_donor_list().

    genome_build (промт "TopMed/HRC", п.3): "grch37" (по умолчанию, HRC) |
        "grch38" (TopMed) — определяет, с каких зеркал/по какому шаблону
        имени файла качаются доноры 1000 Genomes (см. MIRRORS_BY_BUILD/
        VCF_TEMPLATE_BY_BUILD выше) и какой префикс CHROM ("" | "chr")
        используется при фильтрации по позициям чипа. Вызывающий код
        (main.py/gui/app.py) должен передавать сюда
        REFERENCE_PANELS[panel]["genome_build"] — иначе для TopMed доноры
        по-прежнему тянулись бы из GRCh37-релиза 1000 Genomes, что дало бы
        рассинхронизацию координат с GRCh38-референсом/панелью. По
        умолчанию (не передан явно) — GRCh37, поведение идентично версии
        до этого промта.

    bcftools_threads (v18, "ускорение фильтрации уже скачанных файлов"):
        сколько потоков де/компрессии BGZF передавать в каждый вызов
        `bcftools view --threads N` (process_chromosome()/
        process_chromosome_remote()). По умолчанию (None) считается
        автоматически как (число ядер CPU - 1) / max_parallel_chromosomes,
        чтобы не переподписать процессор — параллельно уже обрабатывается
        до max_parallel_chromosomes хромосом, у каждой свой bcftools-
        процесс. Явное число полностью переопределяет автоподбор.

    raw_cache_dir (v16, промт "Доноры для VCF-источника: понятная отмена +
        общий кэш сырых хромосом") — необязательный путь к общему кэшу
        ЕЩЁ НЕ ОТФИЛЬТРОВАННЫХ полных хромосом 1000 Genomes, ОБЩЕМУ для
        всех источников/чипов одной референсной сборки (genome_build),
        обычно donors/_raw_chromosomes/<genome_build>/ (см.
        main.py::raw_chromosome_cache_dir()). По умолчанию None —
        поведение идентично v15 (временный файл полной хромосомы удаляется
        после локальной фильтрации, ничего не кэшируется). При заданном
        пути process_chromosome()/download_chromosome_vcf() сначала ищут
        полную хромосому там (пропуская зеркала целиком, если найдена), а
        после свежего скачивания сохраняют её туда для следующего
        источника/чипа. НЕ используется путём удалённой фильтрации
        (process_chromosome_remote(), Часть 1.1) — там полная хромосома
        никогда не скачивается целиком, кэшировать нечего. Экономит
        трафик и время ЦЕНОЙ постоянного места на диске (~десятки ГБ на
        22 хромосомы 1000 Genomes phase3) — включается по явному желанию
        пользователя (GUI-чекбокс/CLI-флаг), по умолчанию выключено.

    Скачивает и фильтрует доноров 1000 Genomes для ОДНОГО источника/чипа.

    Предназначена для вызова напрямую из приложения (gui/app.py) — без
    отдельного процесса и argparse. Не использует модульные глобалы:
    htslib передаётся параметром (обычно уже существующий pipeline.HTSLIB
    из GUI/CLI), working_suffix_by_mirror — локальная переменная внутри
    этого вызова.

    positions_json: путь к <signature>.positions.json из
        save_position_cache() — сигнатура извлекается из ИМЕНИ файла.
    source: 'ftdna' | 'myheritage' | 'vcf' — только для сообщений в лог,
        разделение по подпапкам делает вызывающий код через output_dir.
    output_dir: donors/<source>/ — сюда пишутся kgp_sub_{1..22}.vcf.gz и
        chip_signature.txt.
    progress_cb(frac, text): вызывается минимум один раз на каждую
        завершённую хромосому (frac в диапазоне 0.0..1.0), текст вида
        "Обработано хромосом: 5/22". Детали (мёртвые зеркала, повторные
        попытки, удалённо/полностью) идут только в print()/лог.
    cancel_check() -> bool: проверяется перед каждой хромосомой И внутри
        отдельных сетевых операций (см. _run_cancelable) — отмена
        прерывает уже идущее скачивание, а не только ждёт следующей
        хромосомы. При параллельной обработке (см. max_parallel_chromosomes)
        отмена останавливает приём новых хромосом в работу и дожидается
        уже запущенных (они сами быстро завершаются, проверяя cancel_check
        изнутри).

    remote_filter: bool = True (Часть 1.1) — пробовать сначала отфильтровать
        доноров прямо по URL без полного скачивания хромосомы (нужен
        bcftools, собранный с libcurl — проверяется автоматически один
        раз перед циклом, см. _probe_bcftools_remote_support()). При
        False или неудачной пробе — как раньше, полное скачивание для
        всех хромосом (само по себе ускоренное за счёт 1.2/1.4).
    max_parallel_chromosomes: int = DEFAULT_PARALLEL_CHROMOSOMES (Часть 1.3) —
        сколько хромосом обрабатывать одновременно. 1 — старое строго
        последовательное поведение v11.
    remote_skip_large_chroms: Optional[Set[int]] = DEFAULT_REMOTE_SKIP_LARGE_CHROMS
        (v14, Шаг 2 промта "точечный патч remote-фильтрации для крупных
        хромосом") — номера хромосом, для которых remote-путь пропускается
        вовсе и сразу используется полное скачивание. По умолчанию
        {1..8} — подтверждено (chr1-3) и экстраполировано (chr4-8) по
        реальному прогону, см. докстринг константы выше. Передайте
        None или set() для поведения как в v13 (пробовать remote для
        всех хромосом).

    Возвращает список путей kgp_sub_{1..22}.vcf.gz (22 штуки при полном
    успехе), ОТСОРТИРОВАННЫЙ по номеру хромосомы — порядок не зависит от
    того, в каком порядке параллельные воркеры завершились.

    ВАЖНО: эта функция сама пишет chip_signature.txt при полном успехе,
    но вызывающий код (gui/app.py) всё равно обязан повторно вызвать
    main.check_donor_cache() после этого вызова — это единая точка
    правды для финального списка путей и единственный надёжный способ
    убедиться, что сигнатура реально совпадает и файлы валидны, без
    дублирования этой проверки здесь.

    Бросает:
      RuntimeError — bcftools не найден, либо не удалось скачать/
        отфильтровать одну из хромосом после всех попыток (включая
        удалённый путь и fallback на полное скачивание).
      DownloadCancelled — cancel_check() вернул True.
    """
    if not htslib.has_bcftools:
        raise RuntimeError(
            "bcftools не найден — укажите --bin-dir/папку с бинарниками "
            "перед скачиванием доноров."
        )

    positions_json = Path(positions_json)
    if not positions_json.exists():
        raise RuntimeError(f"Файл позиций чипа не найден: {positions_json}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def _cb(frac: float, text: str) -> None:
        if progress_cb:
            progress_cb(max(0.0, min(1.0, frac)), text)

    # v13: гарантируем CA-сертификаты для libcurl и предупреждаем о
    # возможном конфликте bin/curl.exe ДО пробы удалённого доступа —
    # раньше это требовало ручной настройки $env:CURL_CA_BUNDLE в каждой
    # новой консоли (см. докстринг core/network_utils.py). Не бросает
    # исключение при неудаче (например, нет интернета) — в этом случае
    # _probe_bcftools_remote_support() ниже просто вернёт False, и
    # прогон продолжится через полное скачивание, как и раньше.
    ensure_network_ready(htslib.bin_dir)

    expected_signature = _extract_signature_from_positions_json(positions_json)
    if expected_signature:
        _invalidate_stale_donor_cache(output_dir, expected_signature, progress_cb)
    else:
        print(
            f"⚠ Не удалось извлечь сигнатуру чипа из имени файла "
            f"{positions_json.name} — проверка устаревшего кэша (Задача 1, "
            f"п.7) пропущена, chip_signature.txt не будет записан."
        )

    print("=" * 70)
    print(f"СКАЧИВАНИЕ ДОНОРОВ 1000 GENOMES (source={source})")
    print(f"Папка: {output_dir.absolute()}")
    print("=" * 70)

    _cb(0.0, "Список европейских образцов...")
    if cancel_check and cancel_check():
        raise DownloadCancelled("Скачивание отменено пользователем")
    eur_samples = create_eur_donor_list(
        output_dir, bin_dir=htslib.bin_dir, eur_sample_count=eur_sample_count,
    )
    eur_set = set(eur_samples)
    print(
        f"ℹ Доноров-образцов в подвыборке: {len(eur_samples)}"
        + (" (вся доступная EUR-панель)" if eur_sample_count is None
           else f" (запрошено: {eur_sample_count})")
    )

    _cb(0.02, "Позиции чипа...")
    if cancel_check and cancel_check():
        raise DownloadCancelled("Скачивание отменено пользователем")
    chip_positions = create_chip_positions_from_json(positions_json, output_dir)

    remote_capable = False
    if remote_filter:
        _cb(0.04, "Проверка удалённого доступа bcftools...")
        if cancel_check and cancel_check():
            raise DownloadCancelled("Скачивание отменено пользователем")
        remote_capable = _probe_bcftools_remote_support(
            htslib, cancel_check, genome_build=genome_build,
        )

    # working_suffix_by_mirror и progress-счётчик разделяются между
    # параллельными воркерами (Часть 1.3) — оба обновляются под общим
    # progress_lock. Гонка за working_suffix_by_mirror сама по себе не
    # критична (это просто подсказка "какой суффикс сработал в прошлый
    # раз", а не источник истины), но пишем под тем же лock'ом для
    # предсказуемости и чтобы не полагаться на детали GIL.
    working_suffix_by_mirror: dict[str, str] = {}
    outputs_by_chrom: dict[int, Path] = {}
    failed_chroms: list[int] = []
    completed_count = 0
    progress_lock = threading.Lock()
    cancelled_flag = threading.Event()

    max_workers = max(1, int(max_parallel_chromosomes))

    # v18: автоподбор --threads для bcftools view, если вызывающий код не
    # задал число явно — делим (ядра CPU - 1) на число параллельно
    # обрабатываемых хромосом, чтобы не переподписать процессор (у каждой
    # параллельной хромосомы свой bcftools-процесс с собственными
    # --threads декомпрессии BGZF).
    if bcftools_threads is None:
        cpu_count = os.cpu_count() or 2
        resolved_bcftools_threads = max(1, (cpu_count - 1) // max_workers)
    else:
        resolved_bcftools_threads = max(1, int(bcftools_threads))
    print(
        f"ℹ bcftools --threads на хромосому: {resolved_bcftools_threads} "
        f"(ядер CPU: {os.cpu_count() or '?'}, параллельно хромосом: {max_workers})"
    )

    skip_set: Set[int] = set(remote_skip_large_chroms) if remote_skip_large_chroms else set()
    if remote_capable and skip_set:
        print(
            f"ℹ Remote-путь по умолчанию пропущен для хромосом "
            f"{sorted(skip_set)} (v14, см. DEFAULT_REMOTE_SKIP_LARGE_CHROMS) — "
            f"для них сразу используется полное скачивание."
        )
    if raw_cache_dir is not None:
        raw_cache_dir = Path(raw_cache_dir)
        raw_cache_dir.mkdir(parents=True, exist_ok=True)
        print(
            f"ℹ Общий кэш сырых хромосом включён: {raw_cache_dir} — полные "
            f"хромосомы 1000 Genomes будут переиспользованы между "
            f"источниками/чипами этой референсной сборки."
        )
    print(
        f"\n[3/3] Обработка хромосом 1-22 (сборка: {genome_build}, "
        f"{'удалённая фильтрация + ' if remote_capable else 'полное скачивание, '}"
        f"до {max_workers} хромосом(ы) параллельно)..."
    )

    def _worker(chrom: int) -> bool:
        if cancelled_flag.is_set() or (cancel_check and cancel_check()):
            raise DownloadCancelled("Скачивание отменено пользователем")
        with progress_lock:
            suffix_map = working_suffix_by_mirror
        return process_chromosome_auto(
            chrom, eur_set, chip_positions, output_dir,
            htslib, suffix_map, remote_capable, cancel_check,
            skip_remote=(chrom in skip_set),
            raw_cache_dir=raw_cache_dir,
            bcftools_threads=resolved_bcftools_threads,
            genome_build=genome_build,
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_chrom = {executor.submit(_worker, c): c for c in range(1, 23)}
        try:
            for future in concurrent.futures.as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    ok = future.result()
                except DownloadCancelled:
                    cancelled_flag.set()
                    raise

                if ok:
                    outputs_by_chrom[chrom] = output_dir / f"kgp_sub_{chrom}.vcf.gz"
                else:
                    failed_chroms.append(chrom)

                with progress_lock:
                    completed_count += 1
                    n = completed_count
                _cb(0.05 + 0.9 * n / 22, f"Обработано хромосом: {n}/22")
        except DownloadCancelled:
            # Ещё не запущенные задачи снимаем из очереди; уже запущенные
            # сами быстро завершатся (проверяют cancelled_flag/cancel_check
            # изнутри, включая poll-цикл _run_cancelable) — `with executor`
            # дождётся их при выходе из блока.
            for f in future_to_chrom:
                f.cancel()
            raise DownloadCancelled("Скачивание отменено пользователем")

    outputs = [outputs_by_chrom[c] for c in sorted(outputs_by_chrom)]

    if failed_chroms:
        raise RuntimeError(
            f"Не удалось скачать/отфильтровать {len(failed_chroms)} из 22 "
            f"хромосом: {', '.join(str(c) for c in sorted(failed_chroms))}. "
            f"Уже готовые хромосомы сохранены — повторный вызов их не "
            f"перекачает (докачка с места обрыва)."
        )

    if expected_signature:
        sig_file = output_dir / "chip_signature.txt"
        sig_file.write_text(expected_signature, encoding="utf-8")
        print(f"✓ chip_signature.txt записан: signature={expected_signature} ({sig_file})")
        _cb(1.0, "Доноры готовы, сигнатура записана")
    else:
        print(
            f"⚠ Сигнатура чипа неизвестна — {output_dir / 'chip_signature.txt'} "
            f"НЕ создан. Повторный вызов check_donor_cache() не примет эти "
            f"доноры, пока сигнатура не появится."
        )
        _cb(1.0, "Доноры скачаны (без сигнатуры)")

    print("=" * 70)
    print(f"ГОТОВО: 22/22 хромосом обработано")
    print("=" * 70)
    return outputs


# ---------------------------------------------------------------------------
# CLI — сохранён для обратной совместимости
# ---------------------------------------------------------------------------
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Скачивание доноров 1000 Genomes (v11)")
    parser.add_argument("--csv", required=False, type=Path, default=None,
                         help="ftdna.csv (нужен, только если --positions-json не указан)")
    parser.add_argument("--positions-json", required=False, type=Path, default=None,
                         help="<signature>.positions.json из save_position_cache() — "
                              "предпочтительный способ, работает для любого "
                              "источника (ftdna/myheritage/vcf), в отличие от --csv. "
                              "При указании этого флага скачивание идёт через "
                              "download_donors_for_chip() с автоматической "
                              "инвалидацией устаревшего кэша при смене чипа.")
    parser.add_argument("--source", choices=["ftdna", "myheritage", "vcf"], default=None,
                         help="Источник данных — используется только для сообщений/подсказок")
    parser.add_argument("--donors-subdir", type=str, default=None,
                         help="Подпапка внутри --output-dir для раздельного хранения "
                              "доноров по источнику, напр. ftdna/myheritage/vcf. "
                              "Без неё поведение остаётся плоским (legacy)")
    parser.add_argument("--output-dir", type=Path, default=Path("donors"))
    parser.add_argument("--bin-dir", type=Path, default=None, help="Папка с bcftools.exe/tabix.exe")
    parser.add_argument(
        "--no-remote-filter", dest="remote_filter", action="store_false", default=True,
        help="Отключить удалённую фильтрацию по URL (Часть 1.1) и всегда качать "
             "хромосомы целиком, как в v11. По умолчанию удалённый доступ "
             "пробуется автоматически (и тихо отключается, если bcftools собран "
             "без libcurl).",
    )
    parser.add_argument(
        "--parallel-chromosomes", type=int, default=DEFAULT_PARALLEL_CHROMOSOMES,
        help=f"Сколько хромосом обрабатывать одновременно (Часть 1.3). "
             f"1 = строго последовательно, как в v11. По умолчанию: "
             f"{DEFAULT_PARALLEL_CHROMOSOMES}.",
    )
    parser.add_argument(
        "--remote-skip-chroms", type=str,
        default=",".join(str(c) for c in sorted(DEFAULT_REMOTE_SKIP_LARGE_CHROMS)),
        help="v14: список номеров хромосом через запятую, для которых "
             "remote-путь (Часть 1.1) пропускается сразу — используется "
             "полное скачивание без бесполезного перебора зеркал/суффиксов. "
             f"По умолчанию: {','.join(str(c) for c in sorted(DEFAULT_REMOTE_SKIP_LARGE_CHROMS))} "
             "(подтверждено/экстраполировано по реальному прогону, см. "
             "DEFAULT_REMOTE_SKIP_LARGE_CHROMS в коде). Передайте пустую "
             "строку '', чтобы пробовать remote-путь для всех хромосом, "
             "как в v13.",
    )
    parser.add_argument(
        "--diagnose-remote", action="store_true", default=False,
        help="v13, Часть 1: не качать доноров, а только проверить, что "
             "удалённая фильтрация bcftools РЕАЛЬНО работает (не только "
             "-h, а полноценный -S/-R запрос на тестовой хромосоме). "
             "Печатает отчёт и завершается с кодом 0 (успех) / 1 (неудача).",
    )
    parser.add_argument(
        "--raw-cache-dir", type=Path, default=None,
        help="v16 (промт 'Доноры для VCF-источника: понятная отмена + общий "
             "кэш сырых хромосом'): путь к общему кэшу ЕЩЁ НЕ "
             "отфильтрованных полных хромосом 1000 Genomes, переиспользуемому "
             "между разными источниками/чипами ОДНОЙ референсной сборки — "
             "обычно donors/_raw_chromosomes/<genome_build>/. По умолчанию "
             "не задан (кэш выключен, поведение как раньше). ⚠ Занимает "
             "~десятки ГБ на диске (22 полные хромосомы 1000 Genomes phase3).",
    )
    parser.add_argument(
        "--genome-build", choices=list(MIRRORS_BY_BUILD.keys()), default=DEFAULT_GENOME_BUILD,
        help="Промт 'TopMed/HRC', п.3: 'grch37' (по умолчанию, HRC) | 'grch38' "
             "(TopMed) — выбирает зеркала/шаблон имени файла GRCh37- или "
             "GRCh38-релиза 1000 Genomes и префикс CHROM ('' | 'chr') при "
             "фильтрации по позициям чипа. main.py/gui/app.py должны "
             "передавать сюда REFERENCE_PANELS[panel]['genome_build'].",
    )
    parser.add_argument(
        "--eur-sample-count", type=int, default=None,
        help="Промт 'Monomorphic sites / настраиваемое количество EUR-"
             "доноров': сколько EUR-образцов 1000 Genomes использовать как "
             "доноров. По умолчанию (не задан) — вся доступная EUR-"
             "подвыборка панели (обычно порядка 500 человек в phase3), что "
             "снижает долю Monomorphic sites на QC MIS по сравнению со "
             "старым фиксированным поведением (20). Укажите явное число "
             f"(например, 20), если трафик/скорость важнее (потолок для "
             f"справки: ~{MAX_EUR_SAMPLE_COUNT}, реальный предел зависит от "
             f"файла панели).",
    )
    args = parser.parse_args()

    bin_dir = args.bin_dir or (Path(os.environ["DONORS_BIN_DIR"]) if os.environ.get("DONORS_BIN_DIR") else None)
    htslib = HtslibTools(bin_dir)

    if args.diagnose_remote:
        # v13: настраиваем CA-сертификаты/предупреждаем о конфликте
        # curl.exe ДО диагностики — иначе диагностика будет тестировать
        # окружение, которое приложение само же ещё не подготовило.
        ensure_network_ready(bin_dir)
        check_dependencies(htslib)
        output_dir = args.output_dir / args.donors_subdir if args.donors_subdir else args.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        print("=" * 70)
        print(f"ДИАГНОСТИКА: реальная удалённая фильтрация bcftools (не только -h), сборка: {args.genome_build}")
        print("=" * 70)
        report = diagnose_remote_filter(htslib, output_dir, genome_build=args.genome_build)
        print(f"URL:              {report['url'] or '(не удалось найти рабочий)'}")
        print(f"Успех:            {'✓ ДА' if report['ok'] else '✗ НЕТ'}")
        print(f"Записей вернулось: {report['records']}")
        print(f"Время:            {report['duration_sec']:.2f}с")
        if report["stderr_tail"]:
            print(f"Диагностика:      {report['stderr_tail']}")
        print("=" * 70)
        sys.exit(0 if report["ok"] else 1)

    if not args.csv and not args.positions_json:
        sys.exit("ОШИБКА: укажите --csv или --positions-json (предпочтительно --positions-json).")

    # v13: единая настройка сетевого окружения (сертификаты + конфликт
    # curl.exe) ДО check_dependencies()/скачивания — см.
    # core/network_utils.py и заметку в начале файла.
    ensure_network_ready(bin_dir)
    check_dependencies(htslib)

    output_dir = args.output_dir / args.donors_subdir if args.donors_subdir else args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"bcftools: {'найден' if htslib.has_bcftools else 'НЕТ'}")
    print(f"tabix: {'найден' if htslib.has_tabix else 'НЕТ'}")

    if args.positions_json:
        # Предпочтительный путь: через новую функцию (Задача 1) — даёт
        # автоматическую инвалидацию устаревшего кэша "бесплатно" и для CLI.
        def _console_progress(frac: float, text: str) -> None:
            print(f"  [{frac * 100:5.1f}%] {text}")

        # v14: пустая строка -> set() (пробовать remote для всех хромосом,
        # поведение v13); непустая строка -> набор номеров хромосом.
        remote_skip_chroms = {
            int(x) for x in args.remote_skip_chroms.split(",") if x.strip()
        } if args.remote_skip_chroms.strip() else set()

        try:
            download_donors_for_chip(
                args.positions_json, args.source or "unknown", output_dir, htslib,
                progress_cb=_console_progress,
                remote_filter=args.remote_filter,
                max_parallel_chromosomes=args.parallel_chromosomes,
                remote_skip_large_chroms=remote_skip_chroms,
                raw_cache_dir=args.raw_cache_dir,
                genome_build=args.genome_build,
                eur_sample_count=args.eur_sample_count,
            )
        except (RuntimeError, DownloadCancelled) as e:
            sys.exit(f"ОШИБКА: {e}")
        return

    # Legacy-путь: чистый --csv без сигнатуры (обратная совместимость со
    # старыми вызовами v9 и раньше). Сигнатура здесь принципиально
    # неизвестна заранее (её ещё нужно посчитать так же, как это делает
    # парсер), поэтому проверка устаревшего кэша не выполняется — как и
    # раньше в v10 при отсутствии --positions-json.
    try:
        eur_samples = create_eur_donor_list(
            output_dir, bin_dir=bin_dir, eur_sample_count=args.eur_sample_count,
        )
    except RuntimeError as e:
        sys.exit(f"ОШИБКА: {e}")
    eur_set = set(eur_samples)
    chip_positions = create_chip_positions(args.csv, output_dir)

    # Легаси-путь (без --positions-json) тоже получает ускорения Части 1:
    # удалённую фильтрацию (если поддерживается) и ускоренную цепочку
    # полного скачивания (aria2c/curl/urllib) внутри process_chromosome_auto()/
    # download_file() — но остаётся строго последовательным (без
    # ThreadPoolExecutor), так как download_donors_for_chip() (и вся
    # инвалидация устаревшего кэша, Задача 1 п.7) для этого пути
    # недоступна — сигнатура чипа здесь принципиально не вычисляется
    # заранее (см. докстринг файла).
    remote_capable = (
        _probe_bcftools_remote_support(htslib, genome_build=args.genome_build)
        if args.remote_filter else False
    )
    working_suffix_by_mirror: dict[str, str] = {}
    print(f"\n[3/3] Обработка хромосом 1-22 (сборка: {args.genome_build})...")
    success = 0
    failed_chroms = []
    for chrom in range(1, 23):
        if process_chromosome_auto(chrom, eur_set, chip_positions, output_dir,
                                    htslib, working_suffix_by_mirror, remote_capable,
                                    genome_build=args.genome_build):
            success += 1
        else:
            failed_chroms.append(chrom)

    print("\n" + "=" * 70)
    print(f"ГОТОВО: {success}/22 хромосом обработано")
    if failed_chroms:
        print(f"✗ Не удалось: {', '.join(str(c) for c in failed_chroms)}")
        print("  Запустите скрипт ещё раз — готовые хромосомы будут пропущены.")
    else:
        print(
            f"⚠ Сигнатура чипа неизвестна (не указан --positions-json) — "
            f"chip_signature.txt НЕ создан. check_donor_cache() не примет эти "
            f"доноры, пока сигнатура не будет проставлена вручную или скрипт "
            f"не будет перезапущен с --positions-json."
        )
    print("=" * 70)


if __name__ == "__main__":
    main()
