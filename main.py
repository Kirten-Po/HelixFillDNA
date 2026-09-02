"""
main.py
Единая точка входа пайплайна конвертации FTDNA/MyHeritage/VCF -> 23andMe (Генотек).

Используется двояко:
  - как библиотека: gui/app.py делает `import main as pipeline` и обращается
    к pipeline.SOURCES, pipeline.HTSLIB, pipeline.ensure_reference_genome(...),
    pipeline.build_vcf(...), pipeline._concat_donors(...), pipeline.check_donor_cache(...)
    и т.д.;
  - как CLI: `python main.py --source ftdna --csv ftdna.csv --template template.txt`
    прогоняет все стадии подряд (кроме шага "загрузить на MIS руками").

ВАЖНО: до этого файла в проекте лежал текстовый документ-патч
(main_py_patch.py), сохранённый под именем main.py. Он не является валидным
кодом (см. IndentationError при запуске) — это были инструкции, что нужно
добавить в "настоящий" main.py. Этот файл — тот самый "настоящий" main.py,
куда фрагменты из патча (автозагрузка референса, санитайзер пароля,
надёжная распаковка ZIP) уже встроены как рабочий код.

=============================================================================
ИЗМЕНЕНИЯ В ЭТОЙ ВЕРСИИ (промт "Раздельное хранение доноров по источникам +
устранение Invalid alleles"):
=============================================================================

Задача A — устранена корневая причина Invalid alleles/950k SNP на MIS:
  раньше `chip_signature.txt` перезаписывался СРАЗУ после парсинга (шаг 1),
  ДО того как _check_donors/check_donor_cache успевал сравнить сигнатуру
  текущего чипа с тем, что реально лежит на диске в donors/. Сравнение
  "cached == chip_signature" в такой последовательности тривиально
  истинно всегда — кэш доноров от другого чипа тихо принимался как
  валидный, bcftools merge честно объединял sample с "не теми" донорами,
  получался union из двух чипов (~950k SNP), и MIS отбрасывал всё, чего
  нет в HRC-панели, как Invalid alleles.

  Исправление: сигнатура чипа теперь ВООБЩЕ не пишется в main.py/gui/app.py
  на этапе парсинга. Она пишется только в download_donors.py — и только
  сразу после успешной свежей загрузки/фильтрации всех 22 донорских
  хромосом. main.py/gui/app.py только ЧИТАЮТ и СРАВНИВАЮТ существующую
  сигнатуру через check_donor_cache() — без права её перезаписывать.

Задача B — раздельное хранение доноров по источникам:
  donors/ftdna/, donors/myheritage/, donors/vcf/, каждая со своим
  kgp_sub_1..22.vcf.gz + chip_signature.txt. Общая функция
  check_donor_cache() используется и CLI (main()), и GUI (_check_donors),
  чтобы логика проверки не дублировалась (по аналогии с core/archive_utils.py).
  Старые "плоские" donors/kgp_sub_*.vcf.gz не удаляются автоматически —
  только предупреждение в логе, если они обнаружены.

Задача C — post_merge_intersect как диагностический/защитный слой:
  после Задач A/B он должен быть no-op ("0 удалено"). Если удаляется
  больше 0 позиций — это сигнал регрессии, и в лог выводится явное
  предупреждение.

=============================================================================
ИЗМЕНЕНИЯ (промт "Детекция несоответствия источника и файла", Задача 2):
=============================================================================

detect_source_from_file(path) -> (source, confidence) — новая функция.
  Определяет вероятный формат файла (ftdna/myheritage/vcf) по его
  содержимому (VCF-заголовок, точный FTDNA-заголовок, '#'-комментарии/
  синонимы колонок MyHeritage, fallback по расширению) и не зависит от
  того, что выбрано в --source/GUI. Читает только "шапку" файла
  (до 40 строк), не парсит его целиком.

  Вызывается ДО "[0/7] Проверка референсного генома" — и в CLI (main(),
  под флагом --auto-detect-source, по умолчанию выключен), и в GUI
  (gui/app.py::App._run_stages_1_6(), всегда включено) — потому что
  ensure_reference_genome() может скачивать/проверять несколько ГБ, а
  источнику 'vcf' референс вообще не нужен: не тратим на это время до
  проверки соответствия источника файлу.

  В CLI при уверенности >= 0.8 и несовпадении с --source — предупреждение
  в stderr, запуск НЕ прерывается (никаких новых обязательных флагов).
  В GUI при том же условии показывается диалог с тремя вариантами
  ("Сменить источник" / "Продолжить с выбранным" / "Отмена") — см.
  gui/app.py::App._prompt_source_mismatch().

=============================================================================
ИЗМЕНЕНИЯ (промт "Поддержка выбора референсной панели HRC / TopMed r3",
Шаг 1 из плана реализации):
=============================================================================

Добавлена единая конфигурация REFERENCE_PANELS (HRC / TopMed) и явный
параметр `panel: str` во всех точках, где раньше подразумевалась только
GRCh37 (HRC): ensure_reference_genome(), _build_reference(),
_donor_source_dir(), _save_chip_signature(), check_donor_cache().

На этом шаге НЕ делается реальный лифтовер координат между GRCh37/GRCh38
(это отдельный следующий шаг плана) — цель шага 1 только в том, чтобы:
  - завести раздельные пути на диске для доноров/сигнатур под каждую
    панель (donors/<source>/<panel>/), чтобы кэш HRC и TopMed никогда не
    перепутались (та же логика защиты, что и в Задаче A/B, но по новой
    оси "панель", а не только "источник");
  - параметризовать скачивание референсного генома под нужную сборку;
  - прокинуть выбор панели через CLI (--panel) и подготовить точки
    расширения для GUI.

ВАЖНО: для source='vcf' референс по-прежнему не нужен (см. _needs_reference),
поэтому panel для 'vcf' влияет только на директорию доноров, а не на
скачивание reference-генома.

=============================================================================
ИЗМЕНЕНИЯ (промт "Поддержка выбора референсной панели HRC / TopMed r3",
следующий шаг — лифтовер координат GRCh37 -> GRCh38):
=============================================================================

Реальный лифтовер координат теперь реализован — core/liftover.py
(ChainLiftover, pure-Python, без внешнего liftOver/CrossMap/
`bcftools +liftover`, который не входит в основной репозиторий bcftools
и живёт отдельно у Giulio Genovese в freeseek/score).

Добавлено:
  - REFERENCE_PANELS["topmed"]["liftover_chain_url"/"liftover_chain_filename"]
    — UCSC chain-файл hg19ToHg38.over.chain.gz (тот же официальный
    источник, которым пользуется сам инструмент UCSC liftOver).
  - ensure_liftover_chain(project_root, panel, progress_cb) — аналог
    ensure_reference_genome(), но для chain-файла (десятки МБ, не
    гигабайты — упрощённая проверка целостности по размеру, без
    полноценного SHA-256-цикла). Возвращает None для панелей без
    "liftover_chain_url" в конфигурации (сейчас — "hrc": лифтовать
    GRCh37 в GRCh37 незачем).
  - _build_liftover(panel, progress_cb) -> Optional[ChainLiftover] —
    качает (если нужно) chain-файл и строит индекс ОДИН раз на запуск.
  - _supports_liftover(source) — по аналогии с _needs_reference():
    источники 'ftdna'/'myheritage' поддерживают параметр liftover в своих
    parser_fn (адаптеры применяют лифт ДО reference.base_at() — см.
    adapters/ftdna_v3.py/myheritage_v5.py); источник 'vcf' — ПОКА НЕТ,
    это отдельная, ещё не реализованная доработка (parse_vcf_source не
    имеет шага reference.base_at(), поэтому лифт для него должен
    применяться иначе — над самим VCF, а не внутри резолвинга ориентации,
    которого там нет). main()/gui/app.py явно предупреждают в лог, если
    выбраны panel="topmed" + source="vcf" одновременно.

main() (CLI) вызывает ensure_liftover_chain()/_build_liftover() сразу
после ensure_reference_genome() (шаг [0b/7]) и передаёт результат в
parser_fn(..., liftover=...) для источников, которые это поддерживают.
gui/app.py::_run_stages_1_6() делает то же самое через pipeline.*.

=============================================================================
ИЗМЕНЕНИЯ (промт "Именованные папки запуска"):
=============================================================================

Раньше ВСЕ файлы одного запуска (sample.vcf.gz, kgp_all.vcf.gz,
batch_merged*.vcf.gz, upload/, parse_result.pkl, rerun_results/,
genotek_23andme_v3*.txt) писались в одну общую output/ — при обработке
второго человека без ожидания писем MIS по первому это тихо перезаписывало
файлы первого запуска, а при перезапуске приложения между Этапом 1-6 и
Этапом 7 терялась привязка "какие файлы к какому человеку".

Теперь каждый запуск получает свою подпапку output/runs/<run_name>/ —
resolve_run_dir() (по умолчанию run_name — следующий свободный номер:
"1", "2", ...), validate_run_name() (запрет недопустимых для имени папки
Windows символов), list_runs()/load_run_info() (для списка "История
запусков" в GUI), save_run_info() (run_info.json с source/panel/
chip_signature/call_rate и т.д. в каждой папке запуска) и
attach_run_log_handler() (дублирование logger.*-сообщений в
<run_dir>/run.log; print()-сообщения дублирует сам вызывающий код — см.
_Tee в CLI main() и LogRedirector в gui/app.py).

donors/<source>/<panel>/ и референсные .fasta — БЕЗ ИЗМЕНЕНИЙ, остаются
общими на все запуски (см. докстринг промта: дорогой для перекачки кэш
не должен дублироваться на каждый запуск).

⚠ Breaking change для CLI: --output-dir теперь КОРЕНЬ (output/), а не
рабочая папка запуска — реальная рабочая папка:
<--output-dir>/runs/<--run-name>/. Новый флаг --run-name (по умолчанию —
автонумерация).

=============================================================================
ИЗМЕНЕНИЯ (фикс "Not a gzipped file" — повреждённый скачанный референс):
=============================================================================

_gunzip_file() раньше не обрабатывал случай, когда скачанный .gz-архив
оказывался повреждён (оборванная докачка, испорченный частичный файл с
прошлого запуска, сбой антивируса/диска и т.п.) — gzip.BadGzipFile
вылетал наружу сырым traceback'ом, а битый .gz оставался на диске.
Из-за этого следующий запуск пытался ДОКАЧАТЬ повреждённый файл (Range-
запрос от текущего размера файла на диске) поверх уже испорченных
данных — та же ошибка повторялась бесконечно, без способа
самостоятельно восстановиться без ручного вмешательства (найти и
удалить файл на диске).

Теперь _gunzip_file() перехватывает ошибки распаковки, удаляет и битый
.gz, и недописанный распакованный файл, и бросает понятный RuntimeError
с инструкцией "запустите ещё раз" — следующий запуск начнёт скачивание
заново с нуля, а не будет докачивать поверх мусора.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import pickle
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
import zlib
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

from adapters.ftdna_v3 import (
    parse_ftdna_v3, ReferenceGenome, FTDNAFormatError, StrandQualityError,
    save_position_cache as _save_position_cache_ftdna,
    save_position_cache_broad as _save_position_cache_broad_ftdna,
)
from adapters.myheritage_v5 import (
    parse_myheritage_v5, MyHeritageFormatError,
    save_position_cache as _save_position_cache_myheritage,
    save_position_cache_broad as _save_position_cache_broad_myheritage,
    COLUMN_SYNONYMS as _MH_COLUMN_SYNONYMS,
    MIN_MATCHED_COLUMNS as _MH_MIN_MATCHED_COLUMNS,
)
from adapters.vcf_source import (
    parse_vcf_source, VCFFormatError,
    save_position_cache as _save_position_cache_vcf,
    _open_text as _vcf_open_text,
)
from adapters.ancestry_v2 import (
    parse_ancestry_v2, AncestryFormatError,
    save_position_cache as _save_position_cache_ancestry,
    save_position_cache_broad as _save_position_cache_broad_ancestry,
    EXPECTED_HEADER as _ANCESTRY_HEADER_TOKENS,
)
from core.ancestry_convert import (
    prepare_ancestry_file, AncestryConvertError, CONVERTED_SUFFIX,
)
from core.pure_python_core import (
    build_vcf, split_autosomes, PureCoreError, _chrom_sort_key, UPLOAD_CHROMS,
    infer_male_from_variants,
)
from core.archive_utils import sanitize_password_text
from core.network_utils import ensure_network_ready
from core.liftover import ChainLiftover, LiftoverError
from template.skeleton import extract_skeleton, SkeletonError
from template.assembler import (
    load_imputed_genotypes, load_measured_genotypes, merge_dictionaries,
    assemble_final, validate_output, AssemblyError,
)
from mis_adapter import MISAdapter, MISAdapterError

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Пути / константы
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DONORS_DIR = PROJECT_ROOT / "donors"
# Промт "Отдельная папка reference/<panel>/" (продолжение промта TopMed/HRC):
# раньше ensure_reference_genome() клала .fasta прямо в корень проекта
# (PROJECT_ROOT / cfg["reference_filename"]) — для двух разных сборок генома
# (GRCh37 у HRC, GRCh38 у TopMed) это смешивало многогигабайтные файлы с
# остальным содержимым корня и никак не группировало их по панели, в
# отличие от donors/<source>/<panel>/. Теперь референсы хранятся в
# reference/<panel>/<reference_filename> — см. _reference_panel_dir() и
# _migrate_legacy_reference() ниже.
REFERENCE_ROOT = PROJECT_ROOT / "reference"
# Промт "HRC / TopMed", лифтовер координат: chain-файлы для лифтовера
# хранятся отдельно от самих референсных .fasta (reference/<panel>/) —
# они не привязаны к панели-получателю данных, а являются "мостом" между
# ДВУМЯ сборками (GRCh37 -> GRCh38), поэтому логичнее не дублировать их
# внутрь reference/topmed/ вместе с многогигабайтным .fasta, а держать в
# соседней подпапке reference/liftover/. См. ensure_liftover_chain().
LIFTOVER_ROOT = REFERENCE_ROOT / "liftover"
IS_WINDOWS = os.name == "nt"

# ---------------------------------------------------------------------------
# Референсные панели импутации (Шаг 1 промта "HRC / TopMed").
#
# Единая точка правды для всего, что отличается между панелями:
#   - display_name    — человекочитаемое имя (используется в GUI/подсказках);
#   - genome_build    — "grch37" | "grch38" — сборка генома этой панели;
#   - chrom_prefix    — префикс хромосомы в итоговом VCF ("" для HRC/GRCh37,
#                        "chr" для TopMed/GRCh38) — пока не применяется нигде,
#                        зарезервировано для следующего шага плана (сборка
#                        chr-файлов под TopMed);
#   - reference_filename / reference_url — какой .fasta скачивать и откуда;
#   - liftover_chain_filename / liftover_chain_url — UCSC chain-файл для
#                        переноса координат ИЗ GRCh37 (все адаптеры парсят
#                        физические координаты чипа именно в GRCh37) В
#                        genome_build этой панели. Отсутствует у "hrc"
#                        (genome_build той же GRCh37 — лифтовать нечего).
#   - donors_source   — просто метка для логов/диагностики;
#   - mis_panel_value — как называется панель на сайте MIS, чтобы можно
#                        было подставлять в текст инструкции пользователю.
#
# ⚠ Проверено при реализации (см. Промт_TopMed_HRC_v2.md, п.1): URL
# референса GRCh38 подтверждён по нескольким независимым источникам
# (документация 1000 Genomes/EBI, GATK, DRAGEN, публикации по
# alt-aware выравниванию на GRCh38) — это официальный FTP 1000 Genomes/
# EBI, откуда этот файл раздают и рекомендуют сами авторы датасета.
# ⚠ В отличие от human_g1k_v37.fasta.gz (HRC), этот файл раздаётся
# НЕСЖАТЫМ (.fa, не .fa.gz) — ensure_reference_genome() учитывает это
# явно (см. её докстринг/тело: скачивание через gzip-архив используется
# только если reference_url оканчивается на ".gz", иначе .fa качается
# напрямую в конечный путь, без промежуточной распаковки).
#
# ⚠ liftover_chain_url (hg19ToHg38.over.chain.gz) — официальный источник
# UCSC, тот же, которым пользуется сам инструмент liftOver. Лифтовер
# реализован (core/liftover.py::ChainLiftover, чистый Python, SNP-only —
# наши данные никогда не содержат инделов, см. adapters/*.py) — см. блок
# докстринга модуля "следующий шаг — лифтовер координат" выше.
# ---------------------------------------------------------------------------
REFERENCE_PANELS = {
    "hrc": {
        "display_name": "HRC r1.1 2016 (GRCh37/hg19)",
        "genome_build": "grch37",
        "chrom_prefix": "",
        "reference_filename": "human_g1k_v37.fasta",
        "reference_url": (
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
            "human_g1k_v37.fasta.gz"
        ),
        # Фикс "скачивание залипает на ~156 МБ из 851 МБ": равнозначные
        # зеркала ОДНОГО И ТОГО ЖЕ файла (побайтово идентичны, поэтому
        # докачка Range-запросом при переключении зеркала корректна —
        # см. _download_with_resume()). Зеркало на Amazon S3 —
        # официальное зеркало 1000 Genomes, уже используемое в этом же
        # проекте для скачивания донорских хромосом (download_donors.py,
        # см. также докстринг core/network_utils.py, где оно фигурирует
        # в ручной диагностике пользователя) — то есть его доступность с
        # машины пользователя уже подтверждена на практике. Порядок
        # важен: S3 идёт первым, потому что именно на EBI-зеркале
        # воспроизводилось залипание.
        # ⚠ Путь внутри S3-бакета 1000genomes НЕ содержит префикса
        # "vol1/ftp/" (в отличие от HTTP-зеркала EBI): корень бакета уже
        # соответствует каталогу ftp/. С "vol1/ftp/" зеркало отвечало
        # 404 Not Found — проверено живым прогоном.
        "reference_urls": [
            "https://1000genomes.s3.amazonaws.com/technical/reference/"
            "human_g1k_v37.fasta.gz",
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
            "human_g1k_v37.fasta.gz",
        ],
        "donors_source": "1000genomes_grch37",
        "mis_panel_value": "HRC r1.1 2016 (GRCh37/hg19)",
        # Промт "поправить ссылку для TopMed": HRC-панель по-прежнему
        # загружается на Michigan Imputation Server (там она и живёт) —
        # поведение не меняется, URL явный, а не жёстко вшит в
        # mis_adapter.py/gui/app.py/main.py, чтобы обе панели брали адрес
        # из одного и того же места конфигурации.
        "mis_upload_url": "https://imputationserver.sph.umich.edu",
        # Намеренно нет "liftover_chain_url" — исходные координаты чипа
        # уже в GRCh37, лифтовать в ту же сборку незачем.
    },
    "topmed": {
        "display_name": "TOPMed r3 (GRCh38/hg38)",
        "genome_build": "grch38",
        "chrom_prefix": "chr",
        # Официальный источник 1000 Genomes/EBI FTP — тот же файл, что
        # используют GATK/DRAGEN/GIAB и рекомендует сама документация
        # 1000 Genomes для alt-aware выравнивания на GRCh38 (см. комментарий
        # выше). Раздаётся НЕСЖАТЫМ — важно для ensure_reference_genome().
        "reference_filename": "GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "reference_url": (
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
            "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        ),
        "donors_source": "1000genomes_grch38",
        "mis_panel_value": "TOPMed r3",
        # Промт "поправить ссылку для TopMed": Michigan Imputation Server
        # (imputationserver.sph.umich.edu) НЕ предоставляет панель TOPMed
        # r3 — она доступна на отдельном сервисе того же семейства
        # (тот же движов Cloudgene/eMIS, тот же формат письма/curl-команды
        # для скачивания результатов, поэтому download_mis_results_smart()/
        # extract_all_results() не требуют изменений — отличается только
        # адрес страницы ЗАГРУЗКИ файлов). Подтверждено пользователем на
        # реальном прогоне: раньше сюда открывался Michigan, где TOPMed r3
        # отсутствует в списке Reference Panel.
        "mis_upload_url": "https://imputation.biodatacatalyst.nhlbi.nih.gov/",
        # Промт "HRC / TopMed", лифтовер координат: UCSC chain-файл
        # GRCh37 (hg19) -> GRCh38 (hg38). Официальный источник — тот же,
        # которым пользуется сам инструмент UCSC liftOver.
        "liftover_chain_filename": "hg19ToHg38.over.chain.gz",
        "liftover_chain_url": (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
            "hg19ToHg38.over.chain.gz"
        ),
        # Промт "встроить лифтовер HRC/TopMed в gui/app.py": обратный
        # chain-файл GRCh38 (hg38) -> GRCh37 (hg19) — нужен на Этапе 7,
        # чтобы перенести координаты результата Michigan Imputation Server
        # (который для TopMed приходит в GRCh38) обратно в GRCh37, где
        # живёт скелет трафарета (template/skeleton.py). Официальный
        # источник — тот же UCSC goldenPath, симметричный форвардному.
        "liftover_chain_reverse_filename": "hg38ToHg19.over.chain.gz",
        "liftover_chain_reverse_url": (
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
            "hg38ToHg19.over.chain.gz"
        ),
    },
}
DEFAULT_PANEL = "hrc"

# Задача 3: автозагрузка референсного генома (обратная совместимость —
# старые константы теперь являются алиасами на конфигурацию HRC-панели;
# новый код должен обращаться к REFERENCE_PANELS[panel], а не напрямую
# к этим именам).
REFERENCE_FILENAME = REFERENCE_PANELS[DEFAULT_PANEL]["reference_filename"]
REFERENCE_URL = REFERENCE_PANELS[DEFAULT_PANEL]["reference_url"]
REFERENCE_MIN_SIZE = int(2.5 * 1024 ** 3)  # 2.5 ГБ — минимальный ожидаемый размер .fasta
REFERENCE_SHA256_SIDECAR_SUFFIX = ".sha256"
# Промт "HRC / TopMed", лифтовер координат: chain-файлы UCSC (десятки МБ)
# на порядки меньше референсных .fasta (гигабайты) — полноценный SHA-256-
# цикл с sidecar-файлом здесь избыточен (см. ensure_liftover_chain()),
# минимальный размер используется только как защита от битой/HTML-страницы
# вместо настоящего файла.
LIFTOVER_CHAIN_MIN_SIZE = 100 * 1024  # 100 КБ

# Источники данных, которые понимает пайплайн.
# Ключ используется и в GUI (выпадающий список), и в CLI (--source), и как
# имя подпапки в donors/<source>/<panel>/ (Задача B + Шаг 1 промта TopMed).
#
# "save_position_cache" — Задача B, доп. пункт: вместо повторного парсинга
# сырого CSV в download_donors.py, сохраняем позиции чипа СРАЗУ после
# парсинга тем же самым парсером, который уже используется в пайплайне
# (устраняет целый класс багов для MyHeritage/VCF источников с их гибкими
# форматами).
SOURCES = {
    "ftdna": {
        "name": "FTDNA Family Finder (.csv)",
        "parser": parse_ftdna_v3,
        "save_position_cache": _save_position_cache_ftdna,
        "save_position_cache_broad": _save_position_cache_broad_ftdna,
    },
    "myheritage": {
        "name": "MyHeritage (.csv)",
        "parser": parse_myheritage_v5,
        "save_position_cache": _save_position_cache_myheritage,
        "save_position_cache_broad": _save_position_cache_broad_myheritage,
    },
    "ancestry": {
        "name": "AncestryDNA (.txt)",
        "parser": parse_ancestry_v2,
        "save_position_cache": _save_position_cache_ancestry,
        "save_position_cache_broad": _save_position_cache_broad_ancestry,
    },
    "vcf": {
        "name": "Готовый VCF (свой файл / WGS)",
        "parser": parse_vcf_source,
        # У VCF-источника нет person-specific QC-отбраковки (REF/ALT/GT
        # уже даны в файле) — chip_signature и chip_signature_broad были
        # бы идентичны. Отдельная broad-функция не нужна, см.
        # _resolve_chip_signature_mode() ниже.
        "save_position_cache": _save_position_cache_vcf,
        "save_position_cache_broad": None,
    },
}


def _needs_reference(source: str) -> bool:
    """VCF-источнику референс не нужен — REF/ALT/GT там уже разрешены."""
    return source in ("ftdna", "myheritage", "ancestry")


#: Источники, которым перед Этапом 1 нужен отдельный шаг приведения
#: файла к оформлению 23andMe v3 (Этап 0). Пока такой один — AncestryDNA:
#: его сырой экспорт отличается от 23andMe не содержанием, а оформлением
#: (5 колонок, коды хромосом 23-26, пропуск как аллель '0'), и отдельный
#: шаг оставляет на диске промежуточный файл, который можно проверить
#: глазами и залить в Генотек как есть — см. core/ancestry_convert.py.
_SOURCES_NEEDING_CONVERSION = ("ancestry",)


def _default_conversion_template(template_path: Optional[Path] = None) -> Optional[Path]:
    """
    Откуда взять '#'-шапку для конвертированного файла.

    Целевое оформление Этапа 0 — всегда v3, независимо от того, какой
    трафарет выбран для СБОРКИ ИТОГОВОГО файла на Этапе 7 (пользователь
    может выбрать v5, и это не должно менять промежуточный файл).
    Поэтому по умолчанию берётся samples/template_v3.txt, и только если
    его нет — переданный трафарет, а если нет и его, конвертер подставит
    встроенную шапку.
    """
    bundled_v3 = PROJECT_ROOT / "samples" / "template_v3.txt"
    if bundled_v3.is_file():
        return bundled_v3
    if template_path is not None and Path(template_path).is_file():
        return Path(template_path)
    return None


def prepare_source_file(
    source: str,
    csv_path: Path,
    output_dir: Path,
    template_path: Optional[Path] = None,
) -> tuple[Path, Optional[object]]:
    """
    Этап 0. Возвращает (файл_для_парсинга, статистика_конвертации).

    Для источников не из _SOURCES_NEEDING_CONVERSION возвращает исходный
    путь и None — то есть для FTDNA/MyHeritage/VCF поведение пайплайна не
    меняется вообще, вызов этой функции для них бесплатный.

    Для 'ancestry' конвертирует сырой экспорт в
    output_dir/<имя>_23andme_v3.txt и возвращает путь к нему; дальше все
    обычные этапы идут уже по конвертированному файлу. Если на вход
    подсунули УЖЕ конвертированный файл (или любой другой в формате
    23andMe), конвертация пропускается — см. prepare_ancestry_file().
    """
    csv_path = Path(csv_path)
    if source not in _SOURCES_NEEDING_CONVERSION:
        return csv_path, None

    stats = prepare_ancestry_file(
        csv_path, Path(output_dir),
        template_path=_default_conversion_template(template_path),
    )
    return Path(stats.out_path), stats


def _supports_liftover(source: str) -> bool:
    """
    Промт "HRC / TopMed", лифтовер координат: только 'ftdna'/'myheritage'
    адаптеры принимают параметр liftover в своей parser_fn (см.
    adapters/ftdna_v3.py::parse_ftdna_v3()/
    adapters/myheritage_v5.py::parse_myheritage_v5() — лифт применяется
    там ВНУТРИ, до reference.base_at()).

    'vcf' — ПОКА НЕТ. parse_vcf_source() не имеет шага reference.base_at()
    (REF/ALT уже даны в самом VCF), поэтому лифтовер для него должен
    применяться иначе — как отдельный проход над уже распарсенными
    (chrom, pos), а не внутри резолвинга ориентации, которого там просто
    нет. Это отдельная, ещё не реализованная доработка — main()/
    gui/app.py явно предупреждают в лог, если source='vcf' выбран вместе
    с panel="topmed" (координаты в этом случае НЕ переносятся).

    'ancestry' — ДА: parse_ancestry_v2() принимает тот же параметр
    liftover и применяет его в том же месте (сразу после нормализации
    chrom/pos, до broad_key и до reference.base_at()).
    """
    return source in ("ftdna", "myheritage", "ancestry")


def _panel_config(panel: Optional[str]) -> dict:
    """
    Возвращает конфигурацию панели из REFERENCE_PANELS, с явной ошибкой
    при неизвестном ключе (вместо тихого KeyError где-то в середине
    пайплайна) и с откатом на DEFAULT_PANEL при panel=None — для
    обратной совместимости со старым кодом/вызовами, которые ещё не
    прокидывают panel явно.
    """
    key = panel or DEFAULT_PANEL
    if key not in REFERENCE_PANELS:
        raise ValueError(
            f"Неизвестная референсная панель: {key!r}. "
            f"Доступные: {', '.join(REFERENCE_PANELS.keys())}"
        )
    return REFERENCE_PANELS[key]


def _reference_panel_dir(
    panel: Optional[str],
    reference_root: Path = REFERENCE_ROOT,
) -> Path:
    """
    reference_root/<panel>/ — раздельное хранение референсных геномов по
    панели (hrc/topmed), по аналогии с _donor_source_dir() для доноров.

    panel=None откатывается на DEFAULT_PANEL — та же обратная совместимость,
    что и в _panel_config()/_donor_source_dir().
    """
    panel_key = panel or DEFAULT_PANEL
    return Path(reference_root) / panel_key


def _migrate_legacy_reference(legacy_root: Path, panel_dir: Path, cfg: dict) -> Optional[Path]:
    """
    Одноразовая миграция референса со старого пути (PROJECT_ROOT/
    <reference_filename>, до появления папки reference/<panel>/) на новый
    (panel_dir/<reference_filename>).

    В отличие от _warn_if_legacy_flat_donors()/_warn_if_legacy_flat_output()
    (которые только предупреждают и НЕ переносят файлы — там неясно,
    к какому запуску/чипу привязаны старые файлы, трогать их автоматически
    рискованно), здесь миграция БЕЗОПАСНА для автоматического переноса:
    имя референсного файла целиком определяется конфигурацией панели
    (REFERENCE_PANELS), файл ровно один на панель, и он не привязан ни к
    какому конкретному запуску/человеку — переносить нечего перепутать.

    Переносит сам .fasta и его .sha256-sidecar (если есть), пишет
    logger.info о том, что и куда перенесено. Ничего не делает (не бросает
    исключение), если:
      - старого файла на legacy_root нет — нечего мигрировать;
      - новый файл уже существует по новому пути — миграция уже
        выполнялась раньше, повторный вызов идемпотентен.

    Возвращает новый путь к референсу, если миграция произошла, иначе None.
    """
    legacy_path = Path(legacy_root) / cfg["reference_filename"]
    if not legacy_path.exists():
        return None

    new_path = Path(panel_dir) / cfg["reference_filename"]
    if new_path.exists():
        # Новый файл уже на месте (например, миграция уже выполнялась
        # в прошлый запуск) — старый "хвост" в корне проекта не трогаем,
        # пусть пользователь сам решит, удалять ли его вручную.
        return None

    panel_dir = Path(panel_dir)
    panel_dir.mkdir(parents=True, exist_ok=True)

    legacy_sidecar = legacy_path.with_suffix(
        legacy_path.suffix + REFERENCE_SHA256_SIDECAR_SUFFIX
    )
    new_sidecar = new_path.with_suffix(new_path.suffix + REFERENCE_SHA256_SIDECAR_SUFFIX)

    logger.info(
        "ℹ Обнаружен референс старого формата (%s) — переношу в %s "
        "(новая структура хранения reference/<панель>/, см. docstring "
        "_migrate_legacy_reference).", legacy_path, new_path,
    )
    try:
        legacy_path.rename(new_path)
    except OSError:
        # Разные разделы диска/файловые системы — rename() не работает
        # между ними, откатываемся на copy+delete (аналог shutil.move()).
        shutil.move(str(legacy_path), str(new_path))

    if legacy_sidecar.exists():
        try:
            legacy_sidecar.rename(new_sidecar)
        except OSError:
            shutil.move(str(legacy_sidecar), str(new_sidecar))

    logger.info("✓ Референс перенесён: %s", new_path)
    return new_path


# ---------------------------------------------------------------------------
# Задача 2: детекция несоответствия источника (--source / выпадающий
# список GUI) и реального формата файла — до тяжёлых операций.
# ---------------------------------------------------------------------------
_FTDNA_HEADER_TOKENS = ("RSID", "CHROMOSOME", "POSITION", "RESULT")
_MYHERITAGE_MIN_COMMENT_LINES = 10
# Не читаем файл целиком — только "шапку". Референсный геном (до нескольких
# ГБ) и сам парсинг ещё не запущены на этом этапе, поэтому важно не тратить
# время/память на файлы, которые могут быть очень большими (WGS VCF).
_DETECT_MAX_SCAN_LINES = 40


def _detect_clean_header_token(token: str) -> str:
    """
    Нормализация токена заголовка для сравнения с COLUMN_SYNONYMS:
    убрать кавычки/пробелы по краям, привести к нижнему регистру.
    Логика намеренно продублирована в миниатюре (а не импортирована как
    приватная _clean_header_token() из adapters/myheritage_v5.py) —
    единственное, что здесь реально нужно от того модуля, это публичная
    константа COLUMN_SYNONYMS, тянуть привязку к приватному API соседнего
    адаптера ради одной строки нормализации не хочется.
    """
    v = token.strip()
    if len(v) >= 2 and v[0] == v[-1] and v[0] in ('"', "'"):
        v = v[1:-1]
    return v.strip().lower()


def _looks_like_ftdna_header(line: str) -> bool:
    """Точное совпадение с 'RSID,CHROMOSOME,POSITION,RESULT' (как того
    требует adapters/ftdna_v3.py::_validate_header) — без синонимов и
    гибкости, это как раз то, что отличает FTDNA от MyHeritage."""
    tokens = tuple(part.strip() for part in line.strip().split(","))
    return tokens == _FTDNA_HEADER_TOKENS


def _myheritage_header_synonym_score(tokens: list[str]) -> int:
    """Сколько из 4 канонических колонок (RSID/CHROMOSOME/POSITION/RESULT)
    нашлось среди токенов заголовка через adapters.myheritage_v5.COLUMN_SYNONYMS."""
    cleaned = [_detect_clean_header_token(t) for t in tokens]
    matched = 0
    for _canonical, synonyms in _MH_COLUMN_SYNONYMS.items():
        if any(tok in synonyms for tok in cleaned):
            matched += 1
    return matched


def detect_source_from_file(path: Path) -> tuple[Optional[str], float]:
    """
    Определяет вероятный источник
    ('ftdna' | 'myheritage' | 'ancestry' | 'vcf' | None)
    по СОДЕРЖИМОМУ файла, независимо от того, что выбрано в GUI/--source.

    Вызывается ДО тяжёлых операций пайплайна (в первую очередь до
    ensure_reference_genome(), которая может качать/проверять несколько
    ГБ и вообще не нужна источнику 'vcf') — поэтому читает только "шапку"
    файла (до _DETECT_MAX_SCAN_LINES строк), а не весь файл целиком.

    Возвращает (source, confidence), confidence в диапазоне 0.0..1.0.

    Правила проверяются по порядку, первое совпадение побеждает:
      1. VCF: первая непустая строка начинается с '##fileformat=VCF'
         -> ('vcf', 1.0). Для '*.vcf.gz' файл открывается ТЕМ ЖЕ
         способом, что и adapters.vcf_source._open_text() (gzip.open в
         текстовом режиме при расширении '.gz') — иначе для сжатых VCF
         детекция всегда бы молчала.
      2. FTDNA: САМАЯ первая непустая строка файла (то есть перед ней
         нет ни одной '#'-строки-комментария) точно равна
         'RSID,CHROMOSOME,POSITION,RESULT' (после strip каждого поля)
         -> ('ftdna', 0.95).
      3. AncestryDNA: первая непустая строка начинается с
         '#AncestryDNA' -> ('ancestry', 0.98); либо первая
         НЕ-комментарийная строка точно равна
         'rsid\tchromosome\tposition\tallele1\tallele2'
         -> ('ancestry', 0.95).
         ⚠ Проверяется СТРОГО ДО правила MyHeritage: у файлов Ancestry
         18 ведущих '#'-строк, поэтому под правило 4(а) ('10+ строк
         комментариев подряд') они попадают тоже, и без этой проверки
         ЛЮБОЙ Ancestry-файл определялся бы как 'myheritage' с
         уверенностью 0.9.
      4. MyHeritage: в первых _DETECT_MAX_SCAN_LINES строках
           (а) 10+ строк-комментариев '#' подряд в начале файла, ИЛИ
           (б) в первой не-комментарийной непустой строке табов больше,
               чем запятых (TSV), ИЛИ
           (в) среди токенов этой строки через COLUMN_SYNONYMS находится
               минимум MIN_MATCHED_COLUMNS (3) из 4 канонических колонок
         -> ('myheritage', 0.9).
      5. Fallback по расширению файла (если ни одна content-эвристика
         выше не сработала):
           '.vcf' / '.vcf.gz' -> ('vcf', 0.5)
           '.csv'              -> ('ftdna', 0.3) — самый частый CSV-
               источник в проекте, но заведомо более слабый сигнал, чем
               эвристики выше: '.csv' с равным успехом может оказаться и
               MyHeritage, поэтому уверенность специально ниже, чем у
               content-эвристик (используется чисто для info-сообщения в
               лог — при confidence < 0.8 предупреждение/диалог не
               показываются, см. вызывающий код).
         Прочие расширения -> (None, 0.0).

    Если файл нечитаем (проблема кодировки, прав доступа, файл пуст и
    т.п.) — исключение НЕ пробрасывается наружу: пишется предупреждение
    в лог, возвращается (None, 0.0).
    """
    path = Path(path)
    try:
        if not path.exists() or path.stat().st_size == 0:
            logger.warning("Автодетект источника: файл не найден или пуст: %s", path)
            return None, 0.0

        lines: list[str] = []
        with _vcf_open_text(path) as f:
            for _ in range(_DETECT_MAX_SCAN_LINES):
                raw_line = f.readline()
                if not raw_line:
                    break
                lines.append(raw_line.rstrip("\r\n"))

        non_empty = [l for l in lines if l.strip() != ""]
        if not non_empty:
            logger.warning(
                "Автодетект источника: не удалось прочитать содержимое файла: %s", path
            )
            return None, 0.0

        # --- 1. VCF -----------------------------------------------------
        if non_empty[0].lstrip("\ufeff").startswith("##fileformat=VCF"):
            return "vcf", 1.0

        # --- 2. FTDNA (только если это САМАЯ первая строка файла, то
        #        есть перед заголовком нет '#'-комментариев) ------------
        first_raw = lines[0].lstrip("\ufeff") if lines else ""
        if first_raw.strip() and not first_raw.lstrip().startswith("#"):
            if _looks_like_ftdna_header(first_raw):
                return "ftdna", 0.95

        # --- 3. AncestryDNA (СТРОГО до MyHeritage, см. докстринг) -----
        if non_empty[0].lstrip("\ufeff").lower().startswith("#ancestrydna"):
            return "ancestry", 0.98
        for line in non_empty:
            if line.lstrip().startswith("#"):
                continue
            # Первая не-комментарийная строка — либо заголовок Ancestry,
            # либо это другой формат; в обоих случаях дальше не смотрим.
            tokens = tuple(
                t.strip().strip('"').lower() for t in line.split("\t")
            )
            if tokens == _ANCESTRY_HEADER_TOKENS:
                return "ancestry", 0.95
            break

        # --- 4. MyHeritage -------------------------------------------------
        leading_comment_lines = 0
        header_candidate: Optional[str] = None
        for line in non_empty:
            if line.lstrip().startswith("#"):
                leading_comment_lines += 1
                continue
            header_candidate = line
            break

        if leading_comment_lines >= _MYHERITAGE_MIN_COMMENT_LINES:
            return "myheritage", 0.9

        if header_candidate is not None:
            tabs = header_candidate.count("\t")
            commas = header_candidate.count(",")
            if tabs > commas:
                return "myheritage", 0.9

            delimiter = "," if commas >= tabs else "\t"
            tokens = header_candidate.split(delimiter)
            if _myheritage_header_synonym_score(tokens) >= _MH_MIN_MATCHED_COLUMNS:
                return "myheritage", 0.9

        # --- 5. Fallback по расширению --------------------------------------
        name_lower = path.name.lower()
        if name_lower.endswith(".vcf.gz") or name_lower.endswith(".vcf"):
            return "vcf", 0.5
        if name_lower.endswith(".csv"):
            return "ftdna", 0.3

        return None, 0.0

    except Exception as e:  # noqa: BLE001 — детекция не должна ронять пайплайн
        logger.warning("Автодетект источника не удался для %s: %s", path, e)
        return None, 0.0


# ---------------------------------------------------------------------------
# Обёртка над htslib-бинарниками (bcftools/tabix/bgzip)
# ---------------------------------------------------------------------------
class HtslibTools:
    def __init__(self, bin_dir: Optional[Path]):
        self.bin_dir = Path(bin_dir) if bin_dir else None
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
    def has_bgzip(self) -> bool:
        return self.bgzip_path is not None

    @property
    def has_tabix(self) -> bool:
        return self.tabix_path is not None

    @property
    def has_bcftools(self) -> bool:
        return self.bcftools_path is not None


# Глобальный инстанс — GUI и CLI переопределяют его после выбора папки bin/
# (см. gui/app.py: `pipeline.HTSLIB = pipeline.HtslibTools(bd)`).
HTSLIB = HtslibTools(None)


# ---------------------------------------------------------------------------
# Прогресс-бар / спиннер для консоли (используются автозагрузкой референса)
# ---------------------------------------------------------------------------
class ProgressBar:
    def __init__(self, total: int, label: str = "", width: int = 40):
        self.total, self.current, self.label, self.width = total, 0, label, width

    def update(self, n: int = 1):
        self.current += n
        pct = self.current / self.total * 100 if self.total > 0 else 0
        filled = int(self.width * self.current / self.total) if self.total > 0 else self.width
        bar = "█" * filled + "░" * (self.width - filled)
        print(f"\r{self.label}: [{bar}] {pct:.1f}% "
              f"({self.current / 1024**2:.1f}/{self.total / 1024**2:.1f} МБ)",
              end="", flush=True)
        if self.current >= self.total:
            print()

    def finish(self):
        if self.current < self.total:
            self.current = self.total
            self.update(0)


class Spinner:
    """Простой спиннер-контекстменеджер для долгих операций без прогресса (%)."""
    def __init__(self, label: str = ""):
        self.label = label

    def __enter__(self):
        print(f"{self.label}...", flush=True)
        return self

    def __exit__(self, *exc):
        print(f"  ✓ {self.label} — готово")
        return False


# ---------------------------------------------------------------------------
# Задача 3: автозагрузка референсного генома
# ---------------------------------------------------------------------------
class _IncompleteDownloadError(Exception):
    """Внутреннее исключение: соединение оборвалось раньше, чем скачался весь файл."""


class _MirrorNotFoundError(Exception):
    """
    Внутреннее исключение: конкретное зеркало вернуло 404 — файла по этому
    адресу нет. В отличие от обрыва связи, повторять попытки к нему
    бессмысленно: зеркало исключается из ротации, чтобы не тратить на него
    оставшиеся попытки (см. _download_with_resume()).
    """


def _download_with_resume(url: str, dest: Path, max_retries: int = 5) -> None:
    """
    Скачивает файл с поддержкой докачки (Range-запрос), с автоматическими
    повторами при обрыве соединения.

    Фикс "тихого обрыва соединения": раньше, если сервер/прокси/нестабильный
    интернет обрывал TCP-соединение раньше времени, response.read() просто
    возвращал пустой b"" — это неотличимо от штатного конца потока, поэтому
    скачивание считалось УСПЕШНЫМ, хотя реально скачалась только часть файла
    (например, 156 из 851 МБ). Прогресс-бар при этом даже рисовал обманчивые
    "100%" (ProgressBar.finish() принудительно дорисовывает 100%, если конец
    так и не был достигнут по счётчику). На следующей стадии распаковка
    такого обрубленного .gz падала с "EOFError: Compressed file ended before
    the end-of-stream marker was reached".

    Теперь после каждой попытки итоговый размер файла сверяется с ожидаемым
    (Content-Length из заголовков) — если файл короче, это считается обрывом
    связи, а не успехом: функция сама докачивает (Range-запрос с текущего
    места) ещё раз, до max_retries попыток с растущей паузой между ними, и
    только если все попытки исчерпаны — поднимает понятную ошибку.

    Фикс "скачивание намертво замедляется на одной и той же отметке":
    url может быть не одной строкой, а СПИСКОМ равнозначных зеркал одного
    и того же файла (см. _reference_urls()/REFERENCE_PANELS[...]["reference_urls"]).
    Подтверждено живым прогоном у пользователя: скачивание с ftp.1000genomes.ebi.ac.uk
    стабильно "залипало" примерно на 156 МБ из 851 МБ (причём реконнект на
    100 МБ не помогал — значит, дело в абсолютной отметке файла/маршруте до
    конкретного зеркала, а не в длительности одного соединения), при этом
    тот же самый файл нормально качался браузером. На каждую следующую
    попытку берётся СЛЕДУЮЩЕЕ зеркало по кругу — уже скачанная часть файла
    при этом не теряется, докачка идёт обычным Range-запросом с той же
    позиции (файл на всех зеркалах побайтово одинаков).
    """
    urls = [url] if isinstance(url, str) else list(url)
    if not urls:
        raise ValueError("_download_with_resume: не передано ни одного URL")

    last_error: Optional[Exception] = None
    dead: set[str] = set()  # зеркала, ответившие 404 — больше не пробуем
    for attempt in range(1, max_retries + 1):
        alive = [u for u in urls if u not in dead] or urls
        current_url = alive[(attempt - 1) % len(alive)]
        try:
            _download_attempt(current_url, dest)
            return
        except _MirrorNotFoundError as e:
            last_error = e
            dead.add(current_url)
            print(f"  Зеркало недоступно (404): {current_url} — пробую следующее.")
            continue
        except (_IncompleteDownloadError, RuntimeError) as e:
            last_error = e
            if attempt < max_retries:
                wait = min(5 * attempt, 30)
                still_alive = [u for u in urls if u not in dead] or urls
                next_url = still_alive[attempt % len(still_alive)]
                print(f"  Загрузка прервалась ({e}). Попытка {attempt}/{max_retries}, "
                      f"докачиваю через {wait} с"
                      + (f" с другого зеркала: {next_url}" if len(still_alive) > 1 else "")
                      + "...")
                time.sleep(wait)
                continue
    raise RuntimeError(
        f"Не удалось скачать файл за {max_retries} попыток: {last_error}\n"
        f"Проверьте подключение к интернету и запустите приложение ещё раз — "
        f"скачивание продолжится с места обрыва."
    ) from last_error


def _reference_urls(cfg: dict) -> list[str]:
    """
    Список равнозначных зеркал для скачивания референса выбранной панели.

    Берётся из cfg["reference_urls"], если он задан, иначе — одиночный
    cfg["reference_url"] (обратная совместимость: старый код/тесты, а также
    панели, для которых зеркала ещё не перечислены, продолжают работать
    без изменений).
    """
    urls = cfg.get("reference_urls")
    if urls:
        return list(urls)
    return [cfg["reference_url"]]


# Промт "виснет на ~18%, хотя тот же файл нормально качается в браузере
# без VPN": подтверждено пользователем — сервер/IP не заблокированы (иначе
# и браузер бы не смог), значит дело не в адресе, а в том, что что-то на
# пути соединения (DPI-оборудование провайдера и т.п.) регулирует именно
# ДОЛГИЙ непрерывный поток на одном TCP/TLS-соединении. Браузеры обычно
# качают большие файлы несколькими последовательными/параллельными Range-
# запросами, а не одним сплошным потоком до конца файла — раньше этот код
# делал именно так (одно соединение до самого конца). Разбивка на куски
# по _DOWNLOAD_CHUNK_BYTES с НОВЫМ HTTP-соединением на каждый кусок — самый
# простой способ вести себя похоже на браузер, не завязываясь на то, что
# именно регулирует поток на стороне провайдера/сети.
_DOWNLOAD_CHUNK_BYTES = 100 * 1024 * 1024  # ~100 МБ на одно HTTP-соединение

# Детект троттлинга: если в течение _SPEED_WINDOW_SECONDS средняя скорость
# держится ниже _MIN_ACCEPTABLE_SPEED_BYTES — считаем, что зеркало/маршрут
# "придушили", и переключаемся на следующее зеркало вместо того, чтобы
# сутками тянуть остаток файла. Порог намеренно низкий (200 КБ/с): цель —
# отсечь именно залипание почти в ноль, а не наказывать медленный, но
# рабочий канал.
_SPEED_WINDOW_SECONDS = 60
_MIN_ACCEPTABLE_SPEED_BYTES = 200 * 1024  # 200 КБ/с


def _download_attempt(url: str, dest: Path) -> None:
    """
    Одна попытка скачивания/докачки файла — используется
    _download_with_resume(). Качает файл последовательными кусками по
    _DOWNLOAD_CHUNK_BYTES, каждый раз заново открывая HTTP-соединение
    (обычный Range-запрос "bytes=<offset>-"), а не одним непрерывным
    потоком до конца файла — см. комментарий у _DOWNLOAD_CHUNK_BYTES.
    """
    existing_size = dest.stat().st_size if dest.exists() else 0
    if existing_size > 0:
        print(f"  Найден частично скачанный файл ({existing_size / 1024**2:.1f} МБ), докачиваю...")

    total: Optional[int] = None
    bar: Optional[ProgressBar] = None

    try:
        while True:
            headers = {"User-Agent": "Mozilla/5.0"}
            if existing_size > 0:
                headers["Range"] = f"bytes={existing_size}-"

            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=45) as response:
                resumed = getattr(response, "status", 200) == 206
                content_length = int(response.headers.get("Content-Length", 0))

                # ⚠ КРИТИЧНО (фикс "скачалось 100%, но архив побит лишними
                # данными в конце"): мы попросили Range, а сервер/зеркало
                # ответило 200 и отдаёт ФАЙЛ ЦЕЛИКОМ с нулевой позиции.
                # Раньше в этом случае полное тело всё равно дописывалось
                # в КОНЕЦ уже скачанного куска ("ab") — получался файл
                # вида <кусок><весь файл>: gzip успешно читал первый член
                # и падал на мусоре после него ("Not a gzipped file"
                # изнутри _read_gzip_header, вызванного из read(), то есть
                # НЕ на первых байтах файла, а при попытке прочитать
                # следующий член). По размеру такой файл выглядел
                # завершённым, поэтому проверка "final_size < total"
                # его пропускала. Теперь при 200-ответе на Range-запрос
                # докачка невозможна в принципе — начинаем файл заново.
                if existing_size > 0 and not resumed:
                    print(
                        "  ⚠ Сервер не поддержал докачку (ответил 200 вместо 206) — "
                        "начинаю файл заново, чтобы не испортить архив."
                    )
                    existing_size = 0
                    total = content_length
                    bar = ProgressBar(total, "  Скачивание файла")
                    bar.current = 0

                if total is None:
                    total = existing_size + content_length if resumed else content_length
                    bar = ProgressBar(total, "  Скачивание файла")
                    bar.current = existing_size if resumed else 0
                mode = "ab" if existing_size > 0 else "wb"
                read_this_connection = 0
                # Детект "залипшей" скорости: соединение формально живо
                # (read() не падает по таймауту, данные идут), но скорость
                # падает до единиц КБ/с — так выглядит троттлинг на стороне
                # сети/зеркала. Ждать в таком режиме бессмысленно: 700 МБ
                # на скорости 5 КБ/с — это больше суток. Прерываем попытку
                # и отдаём управление внешнему ретраю, который возьмёт
                # СЛЕДУЮЩЕЕ зеркало (см. _download_with_resume()).
                window_start = time.monotonic()
                window_bytes = 0
                with open(dest, mode) as f:
                    while read_this_connection < _DOWNLOAD_CHUNK_BYTES:
                        try:
                            chunk = response.read(1024 * 1024)
                        except (OSError, TimeoutError) as e:
                            raise _IncompleteDownloadError(
                                f"обрыв соединения на {bar.current / 1024**2:.1f} МБ ({e})"
                            ) from e
                        if not chunk:
                            break
                        f.write(chunk)
                        bar.update(len(chunk))
                        read_this_connection += len(chunk)

                        window_bytes += len(chunk)
                        elapsed = time.monotonic() - window_start
                        if elapsed >= _SPEED_WINDOW_SECONDS:
                            speed = window_bytes / elapsed
                            if speed < _MIN_ACCEPTABLE_SPEED_BYTES:
                                raise _IncompleteDownloadError(
                                    f"скорость упала до {speed / 1024:.0f} КБ/с на "
                                    f"{bar.current / 1024**2:.1f} МБ — похоже на "
                                    f"троттлинг этого зеркала"
                                )
                            window_start = time.monotonic()
                            window_bytes = 0

            existing_size = dest.stat().st_size
            if total and total > 0 and existing_size >= total:
                break
            if read_this_connection == 0:
                # Соединение ничего не отдало, но по размеру файл ещё не
                # полный — не зацикливаемся бесконечно, отдаём решение
                # внешнему ретраю (_download_with_resume).
                break

        # Сверяем итоговый размер с ожидаемым: тихий обрыв соединения
        # (сервер/DPI закрыл поток раньше времени) выглядит как штатный
        # конец очередного куска, но итоговый файл при этом короче total.
        final_size = dest.stat().st_size
        if total and total > 0 and final_size != total:
            # Строго "!=", а не "<": файл БОЛЬШЕ ожидаемого — такой же
            # признак порчи, как и меньше (лишние данные дописаны в конец,
            # см. комментарий про 200-ответ на Range выше). Раньше проверка
            # "<" такой случай пропускала, и битый архив уходил в распаковку.
            raise _IncompleteDownloadError(
                f"размер не совпал: на диске {final_size / 1024**2:.1f} МБ, "
                f"ожидалось {total / 1024**2:.1f} МБ"
            )
        if bar:
            bar.finish()
    except urllib.error.HTTPError as e:
        if e.code == 416:
            print("  Файл уже скачан полностью")
            return
        if e.code == 404:
            raise _MirrorNotFoundError(
                f"404 Not Found: {url}"
            ) from e
        raise RuntimeError(f"Ошибка HTTP при скачивании: {e.code} {e.reason}") from e
    except (_IncompleteDownloadError, _MirrorNotFoundError):
        raise
    except Exception as e:
        raise RuntimeError(
            f"Не удалось скачать файл: {e}\n"
            f"Проверьте подключение к интернету — запустите ещё раз, "
            f"скачивание продолжится с места обрыва."
        ) from e


def _gunzip_stream(fin, fout) -> tuple[int, int, bool]:
    """
    Распаковывает gzip-поток fin в fout вручную, через zlib, вместо
    gzip.open() + shutil.copyfileobj().

    Зачем не gzip.open (фикс "Not a gzipped file (b'\\x01w')" на файле,
    который скачался ПОЛНОСТЬЮ и правильного размера):

    gzip-файл может состоять из нескольких members подряд (так устроен,
    например, BGZF), поэтому gzip.open() после конца каждого члена
    пытается прочитать заголовок следующего. Если сразу после ПОСЛЕДНЕГО,
    полностью корректного члена в файле лежат посторонние байты, gzip.open
    поднимает BadGzipFile — хотя все полезные данные к этому моменту уже
    прочитаны и записаны. Диагностический признак ровно такой, как в
    отчёте пользователя: исключение приходит из _read_gzip_header,
    вызванного из read() (то есть НЕ на первых байтах файла), сообщение
    показывает первые 2 байта мусора, а размер файла точно совпадает с
    Content-Length сервера.

    Эта функция читает members по очереди и, если очередной member не
    начинается корректно, ЗАВЕРШАЕТСЯ УСПЕШНО (сигнализируя trailing=True),
    но только при условии, что хотя бы один member уже был прочитан
    полностью и что-то было записано. Обрыв ПОСРЕДИ сжатого потока
    (настоящая порча/недокачка) по-прежнему приводит к EOFError — такой
    файл принимать нельзя.

    Возвращает (сколько_байт_записано, сколько_members, был_ли_мусор_в_хвосте).
    """
    buf_size = 1024 * 1024
    written = 0
    members = 0
    pending = b""

    while True:
        if not pending:
            pending = fin.read(buf_size)
            if not pending:
                break  # штатный конец файла

        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            while True:
                out = dec.decompress(pending)
                if out:
                    fout.write(out)
                    written += len(out)
                if dec.eof:
                    pending = dec.unused_data
                    members += 1
                    break
                pending = fin.read(buf_size)
                if not pending:
                    raise EOFError(
                        "архив обрывается посреди сжатого потока "
                        "(файл скачан не полностью или повреждён)"
                    )
        except zlib.error as e:
            if members > 0 and written > 0:
                # Полезные данные уже прочитаны целиком — то, что идёт
                # дальше, к содержимому архива отношения не имеет.
                return written, members, True
            raise EOFError(f"не удалось распаковать первый gzip-member: {e}") from e

    return written, members, False


def _gunzip_file(gz_path: Path, out_path: Path) -> None:
    """
    Распаковывает gz_path в out_path.

    Фикс "Not a gzipped file"/gzip.BadGzipFile: если скачанный .gz-архив
    повреждён (оборванная докачка, испорченный частичный файл с прошлого
    запуска, сбой антивируса/диска и т.п.), раньше исключение вылетало
    наружу сырым traceback'ом, а битый .gz оставался на диске — следующий
    запуск пытался ДОКАЧАТЬ его (Range-запрос от текущего размера файла)
    поверх уже испорченных данных, и ошибка повторялась бесконечно без
    возможности самостоятельно восстановиться.

    Теперь при ошибке распаковки удаляются и битый .gz, и недописанный
    распакованный файл — следующий запуск начнёт скачивание заново с
    нуля (докачка от повреждённого файла всё равно ни к чему хорошему не
    привела бы), и бросается понятный RuntimeError вместо сырого
    traceback'а.
    """
    with Spinner(f"  Распаковка {gz_path.name} (может занять несколько минут)"):
        try:
            with gz_path.open("rb") as fin, out_path.open("wb") as fout:
                written, members, trailing = _gunzip_stream(fin, fout)
            if trailing:
                # Не ошибка: сам gzip-поток прочитан ДО КОНЦА и корректно
                # завершён (см. докстринг _gunzip_stream) — данные целы.
                logger.warning(
                    "⚠ После корректно завершённого gzip-потока в %s обнаружены "
                    "лишние байты — они проигнорированы. Распаковано %.2f ГБ. "
                    "Целостность распакованного файла дополнительно проверяется "
                    "по размеру и SHA-256 (см. ensure_reference_genome).",
                    gz_path.name, written / 1024 ** 3,
                )
        except (gzip.BadGzipFile, EOFError, OSError, zlib.error) as e:
            gz_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)
            raise RuntimeError(
                f"Скачанный архив {gz_path.name} повреждён и не распаковывается "
                f"({e}). Битый файл удалён — запустите ещё раз, чтобы скачать "
                f"его заново с нуля (докачка поверх повреждённого файла привела "
                f"бы к той же ошибке)."
            ) from e


def _check_gzip_magic(gz_path: Path) -> None:
    """
    Быстрая проверка, что скачанный файл вообще начинается как gzip
    (сигнатура 0x1f 0x8b), ДО того как тратить минуты на распаковку 851 МБ
    в ~3 ГБ. Если вместо архива скачалась HTML-страница с ошибкой/заглушка
    провайдера/обрезанный мусор — узнаём об этом сразу.
    """
    with gz_path.open("rb") as f:
        magic = f.read(2)
    if magic != b"\x1f\x8b":
        raise RuntimeError(
            f"Скачанный файл {gz_path.name} не является gzip-архивом "
            f"(первые байты: {magic!r}). Обычно это значит, что вместо файла "
            f"скачалась страница с ошибкой или заглушка."
        )


def _download_and_gunzip_with_retries(
    url, gz_path: Path, out_path: Path, max_retries: int = 3,
) -> None:
    """
    Скачивает gz_path и распаковывает в out_path, автоматически повторяя
    ВЕСЬ цикл (скачивание с нуля + распаковка), если распаковка не удалась.

    Зачем это отдельно от ретраев внутри _download_with_resume(): та
    функция ловит только НЕПОЛНУЮ закачку (оборванное соединение, итоговый
    размер меньше Content-Length). Но бывает и другой случай — файл
    скачивается ПОЛНОСТЬЮ (100%, размер совпадает с ожидаемым), однако его
    содержимое всё равно испорчено с самого начала (например,
    gzip.BadGzipFile: "Not a gzipped file" на первых байтах, а не EOFError
    в конце) — по размеру такая закачка выглядит успешной, но байты не те.

    Частая причина на практике — антивирус или корпоративный прокси с SSL-
    инспекцией, который "на лету" фильтрует/пересобирает большие бинарные
    файлы и портит их, оставляя корректный итоговый размер. Требовать от
    пользователя вручную перезапускать приложение на каждый такой случай
    неудобно — вместо этого сама функция удаляет битый .gz (это уже делает
    _gunzip_file при ошибке) и качает заново с нуля, до max_retries раз.
    """
    last_error: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        # Гарантированно начинаем с чистого листа: остаток от предыдущей
        # неудачной попытки не должен участвовать в докачке (иначе битые
        # байты переживут любое число повторов).
        gz_path.unlink(missing_ok=True)
        _download_with_resume(url, gz_path)
        try:
            _check_gzip_magic(gz_path)
            _gunzip_file(gz_path, out_path)
            return
        except RuntimeError as e:
            last_error = e
            if attempt < max_retries:
                print(
                    f"  ⚠ Распаковка не удалась (попытка {attempt}/{max_retries}) — "
                    f"качаю архив заново с нуля..."
                )
                continue
    raise RuntimeError(
        f"Не удалось скачать и корректно распаковать {gz_path.name} за "
        f"{max_retries} попыток: {last_error}\n"
        f"Файл несколько раз скачивался полностью, но оказывался испорчен — "
        f"похоже, что-то портит его при скачивании (антивирус или "
        f"корпоративный прокси с SSL-инспекцией — попробуйте временно "
        f"отключить проверку HTTPS-трафика для этого приложения или "
        f"скачать через другую сеть)."
    ) from last_error


def _sha256_of_file(path: Path, progress_cb: Optional[Callable[[float, str], None]] = None) -> str:
    """Считает SHA-256 файла потоково (без загрузки целиком в память)."""
    h = hashlib.sha256()
    total = path.stat().st_size or 1
    read = 0
    with path.open("rb") as f:
        while True:
            chunk = f.read(4 * 1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
            read += len(chunk)
            if progress_cb:
                progress_cb(0.95 + 0.05 * min(read / total, 1.0), "Проверка целостности (SHA-256)...")
    return h.hexdigest()


def ensure_reference_genome(
    project_root: Optional[Path] = None,
    progress_cb: Optional[Callable[[float, str], None]] = None,
    panel: str = DEFAULT_PANEL,
) -> Path:
    """
    Гарантирует наличие референсного генома для выбранной панели импутации
    в reference/<panel>/ (см. REFERENCE_ROOT/_reference_panel_dir() и
    докстринг промта "Отдельная папка reference/<panel>/"). Используется и
    из GUI, и из CLI — единая точка правды.

    panel: ключ REFERENCE_PANELS ("hrc" | "topmed"). Имя файла и URL
        скачивания берутся из REFERENCE_PANELS[panel] — так что для HRC
        и TopMed это РАЗНЫЕ файлы на диске (не конфликтуют между собой,
        не надо перекачивать один при переключении на другой и обратно).

    project_root: определяет и корень для legacy-пути (где раньше мог
        лежать референс — PROJECT_ROOT/<reference_filename>, до появления
        reference/<panel>/), и корень, от которого строится
        reference_root = project_root/"reference". Если не передан —
        берётся PROJECT_ROOT (и, соответственно, модульная константа
        REFERENCE_ROOT). Параметр оставлен ради переопределения в тестах
        (изолированная временная директория) без изменения глобального
        состояния модуля — при вызове без аргумента поведение идентично
        обращению к REFERENCE_ROOT напрямую.

    ⚠ Скачивание автоматически определяет, сжат ли источник (см. ниже,
    ветвление по cfg["reference_url"].endswith(".gz")) — HRC раздаётся как
    .fasta.gz (нужна распаковка), а подтверждённый источник GRCh38 для
    TopMed раздаётся уже несжатым .fa (распаковка не требуется).
    """
    cfg = _panel_config(panel)

    def notify(p: float, text: str) -> None:
        if progress_cb:
            progress_cb(p, text)

    legacy_root = Path(project_root) if project_root else PROJECT_ROOT
    reference_root = legacy_root / "reference"
    panel_dir = _reference_panel_dir(panel, reference_root)
    panel_dir.mkdir(parents=True, exist_ok=True)

    # Одноразовая миграция референса со старого "плоского" пути
    # (legacy_root/<reference_filename>) в panel_dir — см. докстринг
    # _migrate_legacy_reference(). Делается ДО проверки существования
    # файла по новому пути, чтобы уже скачанный референс не перекачивался
    # заново только из-за смены структуры хранения.
    _migrate_legacy_reference(legacy_root, panel_dir, cfg)

    ref_path = panel_dir / cfg["reference_filename"]
    sidecar_path = ref_path.with_suffix(ref_path.suffix + REFERENCE_SHA256_SIDECAR_SUFFIX)

    notify(0.0, f"Проверка наличия референса ({cfg['display_name']})...")

    if ref_path.exists() and ref_path.stat().st_size >= REFERENCE_MIN_SIZE:
        print(f"✓ Референс найден: {ref_path.name} ({ref_path.stat().st_size / 1024**3:.2f} ГБ)")
        notify(0.9, "Проверка целостности (SHA-256)...")
        actual_hash = _sha256_of_file(ref_path, progress_cb=notify)
        if sidecar_path.exists():
            expected_hash = sidecar_path.read_text(encoding="utf-8").strip()
            if actual_hash != expected_hash:
                raise RuntimeError(
                    f"Референс {ref_path.name} не прошёл проверку SHA-256: "
                    f"файл изменился или повреждён с прошлого запуска "
                    f"(было {expected_hash[:16]}…, стало {actual_hash[:16]}…). "
                    f"Удалите {ref_path.name} и {sidecar_path.name} и запустите заново, "
                    f"чтобы перекачать референс."
                )
            print(f"✓ SHA-256 совпадает с сохранённым: {actual_hash}")
        else:
            sidecar_path.write_text(actual_hash, encoding="utf-8")
            print(f"ℹ SHA-256 референса зафиксирован для последующих проверок: {actual_hash}")
        notify(1.0, "Референс готов")
        return ref_path

    if ref_path.exists():
        print(
            f"⚠ Референс {ref_path.name} повреждён или неполон "
            f"({ref_path.stat().st_size / 1024**3:.2f} ГБ) — перекачиваю"
        )
        ref_path.unlink(missing_ok=True)
        sidecar_path.unlink(missing_ok=True)

    is_gzipped_source = cfg["reference_url"].endswith(".gz")
    print(
        f"ℹ Референс не найден ({cfg['display_name']}), начинаю скачивание "
        f"(размер архива/распакованного файла зависит от сборки генома)"
    )
    print(f"  Источник: {cfg['reference_url']}")
    notify(0.1, "Скачивание...")

    if is_gzipped_source:
        # HRC (human_g1k_v37.fasta.gz) — источник сжат gzip'ом, качаем во
        # временный .gz и распаковываем в ref_path.
        gz_path = panel_dir / (cfg["reference_filename"] + ".gz")
        notify(0.8, "Скачивание и распаковка...")
        print("Распаковываю референс (gzip)...")
        _download_and_gunzip_with_retries(_reference_urls(cfg), gz_path, ref_path)
        gz_path.unlink(missing_ok=True)
    else:
        # TopMed (GRCh38_full_analysis_set_plus_decoy_hla.fa) — официальный
        # источник 1000 Genomes/EBI раздаёт файл НЕСЖАТЫМ (проверено при
        # реализации, см. комментарий у REFERENCE_PANELS["topmed"]) —
        # качаем сразу в ref_path, распаковка не нужна и не выполняется.
        _download_with_resume(_reference_urls(cfg), ref_path)
        notify(0.8, "Скачивание завершено (несжатый источник, распаковка не требуется)")

    if not ref_path.exists() or ref_path.stat().st_size < REFERENCE_MIN_SIZE:
        size_gb = ref_path.stat().st_size / 1024 ** 3 if ref_path.exists() else 0
        raise RuntimeError(
            f"После распаковки размер референса подозрительно мал "
            f"({size_gb:.2f} ГБ, ожидалось >= 2.5 ГБ). Запустите ещё раз — "
            f"скачивание докачается с места обрыва."
        )

    notify(0.95, "Проверка целостности (SHA-256)...")
    actual_hash = _sha256_of_file(ref_path, progress_cb=notify)
    sidecar_path.write_text(actual_hash, encoding="utf-8")
    print(f"ℹ SHA-256 референса зафиксирован для последующих проверок: {actual_hash}")

    print(f"✓ Референс готов: {ref_path} ({ref_path.stat().st_size / 1024**3:.2f} ГБ)")
    notify(1.0, "Референс готов")
    return ref_path


# ---------------------------------------------------------------------------
# Промт "HRC / TopMed", лифтовер координат: автозагрузка UCSC chain-файла.
# ---------------------------------------------------------------------------
def ensure_liftover_chain(
    project_root: Optional[Path] = None,
    panel: str = DEFAULT_PANEL,
    progress_cb: Optional[Callable[[float, str], None]] = None,
    direction: str = "forward",
) -> Optional[Path]:
    """
    Аналог ensure_reference_genome(), но для chain-файла лифтовера
    (reference/liftover/<liftover_chain_filename>).

    direction (промт "встроить лифтовер HRC/TopMed в gui/app.py"):
      "forward" (по умолчанию) — GRCh37 -> сборка панели, читает
          "liftover_chain_url"/"liftover_chain_filename" (как раньше,
          используется при парсинге исходного файла на Этапе 1).
      "reverse" — сборка панели -> GRCh37, читает
          "liftover_chain_reverse_url"/"liftover_chain_reverse_filename"
          (используется на Этапе 7, чтобы перенести результат Michigan
          Imputation Server обратно в GRCh37, где живёт скелет трафарета).

    Возвращает None для панелей/направлений, у конфигурации которых нет
    соответствующего URL (сейчас — "hrc" вообще не имеет ни одного из
    двух: исходные координаты чипа уже в GRCh37, лифтовать в ту же
    сборку незачем — см. REFERENCE_PANELS). Вызывающий код
    (_build_liftover() ниже) трактует None как "лифтовер для этой
    панели/направления не нужен", а не как ошибку.

    В отличие от ensure_reference_genome() (файл на 2.5+ ГБ, нужна
    докачка по факту частично скачанного размера, обязательная SHA-256-
    проверка при каждом запуске) chain-файл — десятки МБ, качается за
    секунды даже на медленном канале. Здесь используется упрощённая
    проверка: если файл уже есть и его размер не подозрительно мал
    (>= LIFTOVER_CHAIN_MIN_SIZE) — считаем его валидным без повторного
    SHA-256 на каждый запуск (тот же принцип, что и в
    core/network_utils.py::ensure_ca_bundle() для cacert.pem — файл
    маленький, риск незамеченной порчи на диске уже скачанного файла
    несопоставимо ниже риска "забыли докачать/перекачали не до конца",
    что как раз ловит проверка размера).
    """
    if direction not in ("forward", "reverse"):
        raise ValueError(f"Неизвестное направление лифтовера: {direction!r} (ожидается 'forward'/'reverse')")

    cfg = _panel_config(panel)
    url_key = "liftover_chain_url" if direction == "forward" else "liftover_chain_reverse_url"
    filename_key = "liftover_chain_filename" if direction == "forward" else "liftover_chain_reverse_filename"
    url = cfg.get(url_key)
    if not url:
        return None

    def notify(p: float, text: str) -> None:
        if progress_cb:
            progress_cb(p, text)

    root = Path(project_root) if project_root else PROJECT_ROOT
    chain_dir = root / "reference" / "liftover"
    chain_dir.mkdir(parents=True, exist_ok=True)
    chain_path = chain_dir / cfg[filename_key]

    if chain_path.exists() and chain_path.stat().st_size >= LIFTOVER_CHAIN_MIN_SIZE:
        print(f"✓ Chain-файл лифтовера найден: {chain_path.name} "
              f"({chain_path.stat().st_size / 1024**2:.1f} МБ)")
        notify(1.0, "Chain-файл готов")
        return chain_path

    if chain_path.exists():
        print(
            f"⚠ Chain-файл {chain_path.name} подозрительно мал "
            f"({chain_path.stat().st_size} байт) — перекачиваю"
        )
        chain_path.unlink(missing_ok=True)

    print(f"ℹ Chain-файл лифтовера не найден ({cfg['display_name']}), скачиваю")
    print(f"  Источник: {url}")
    notify(0.0, "Скачивание chain-файла лифтовера...")
    _download_with_resume(url, chain_path)

    if not chain_path.exists() or chain_path.stat().st_size < LIFTOVER_CHAIN_MIN_SIZE:
        size = chain_path.stat().st_size if chain_path.exists() else 0
        raise RuntimeError(
            f"Chain-файл {chain_path.name} скачался подозрительно маленьким "
            f"({size} байт, ожидалось >= {LIFTOVER_CHAIN_MIN_SIZE} байт) — "
            f"похоже, вместо файла скачалась HTML-страница с ошибкой. "
            f"Проверьте подключение к интернету и запустите ещё раз."
        )

    print(f"✓ Chain-файл готов: {chain_path} ({chain_path.stat().st_size / 1024**2:.1f} МБ)")
    notify(1.0, "Chain-файл готов")
    return chain_path


def _build_liftover(
    panel: str,
    project_root: Optional[Path] = None,
    progress_cb: Optional[Callable[[float, str], None]] = None,
    direction: str = "forward",
) -> Optional[ChainLiftover]:
    """
    None для панелей/направлений без соответствующего chain URL (сейчас —
    'hrc' целиком, для 'topmed' оба направления заданы). Для 'topmed' —
    качает (если нужно) chain-файл через ensure_liftover_chain() и строит
    ChainLiftover ОДИН раз на запуск: сам ChainLiftover парсит файл и
    строит bisect-индекс в конструкторе (см. core/liftover.py) —
    повторного парсинга на каждую позицию/каждый вызов парсера не
    происходит, объект переиспользуется на все ~700-900 тыс. позиций
    типичного чипа за один вызов.

    direction (промт "встроить лифтовер HRC/TopMed в gui/app.py"):
      "forward" (по умолчанию, поведение не изменилось) — GRCh37 -> сборка
          панели, используется при парсинге исходного файла (Этап 1).
      "reverse" — сборка панели -> GRCh37, используется на Этапе 7 для
          переноса результата Michigan Imputation Server обратно в GRCh37.

    Бросает LiftoverError, если chain-файл скачался, но не парсится
    (нестандартный формат — например, tStrand='-', см. докстринг
    core/liftover.py) — это фатальная ошибка конфигурации, а не что-то,
    от чего можно тихо откатиться (в отличие от ensure_network_ready(),
    где неудача просто означает "удалённая фильтрация недоступна, работаем
    как раньше" — здесь без лифтовера panel="topmed" в принципе не может
    дать корректный результат, поэтому лучше упасть явно и сразу).
    """
    chain_path = ensure_liftover_chain(
        project_root=project_root, panel=panel, progress_cb=progress_cb, direction=direction,
    )
    if chain_path is None:
        return None
    return ChainLiftover(chain_path)


def liftover_positions_forward(
    positions: list[tuple[str, int]],
    liftover: ChainLiftover,
) -> list[tuple[str, int]]:
    """
    Промт "встроить лифтовер HRC/TopMed в gui/app.py": лифтует список
    позиций ВПЕРЁД (GRCh37 -> сборка панели, тем же chain-файлом, что и
    parser_fn на Этапе 1) — используется в gui/app.py::_run_stage_7() для
    переноса координат скелета трафарета (всегда GRCh37,
    template/skeleton.py) в сборку панели ПЕРЕД фильтрацией результата
    MIS через template/assembler.py::load_imputed_genotypes(..., panel_pos=...),
    которая применяет panel_pos как прямой позиционный фильтр (-R) к VCF
    результата импутации — а этот VCF для не-HRC панелей приходит уже в
    сборке панели, не в GRCh37.

    Позиции, которые лифтовер не смог перенести (нет chain-блока для
    хромосомы, gap между блоками выравнивания — см.
    core/liftover.py::ChainLiftover.lift()), молча отбрасываются: они
    физически не могут совпасть ни с одной позицией результата MIS в
    сборке панели, включать их в фильтр бессмысленно.

    Хромосома в возвращаемых кортежах — КАНОНИЧЕСКАЯ (без префикса
    "chr"), как её всегда возвращает ChainLiftover.lift(). Если формат
    CHROM в самом VCF результата MIS требует префикса (для TopMed это
    обычно так, см. REFERENCE_PANELS["topmed"]["chrom_prefix"]) —
    сопоставление префикса остаётся заботой вызывающего кода/самого
    load_imputed_genotypes(), эта функция сознательно не знает и не
    должна знать формат CHROM в конкретном VCF результата.
    """
    lifted_positions: list[tuple[str, int]] = []
    dropped = 0
    for chrom, pos in positions:
        result = liftover.lift(str(chrom), int(pos))
        if result is not None:
            lifted_positions.append(result)
        else:
            dropped += 1
    if dropped:
        logger.info(
            "liftover_positions_forward: %d из %d позиций скелета не "
            "перенесены вперёд на сборку панели — исключены из фильтра "
            "результата MIS.",
            dropped, len(positions),
        )
    return lifted_positions


def liftback_imputed_genotypes(
    genotypes: dict[str, str],
    liftover: ChainLiftover,
) -> tuple[dict[str, str], int]:
    """
    Промт "встроить лифтовер HRC/TopMed в gui/app.py": переносит словарь
    импутированных генотипов от
    template/assembler.py::load_imputed_genotypes() (ключ
    f"{chrom}_{pos}", координаты сборки панели) ОБРАТНО в GRCh37 —
    используется в gui/app.py::_run_stage_7() ПЕРЕД
    template/assembler.py::merge_dictionaries()/assemble_final(), которые
    ищут генотип по ключу f"{chrom}_{pos}" в координатах скелета
    трафарета (всегда GRCh37).

    Хромосома в ключах входного словаря принимается в любом виде ("1"
    или "chr1" — префикс, если есть, снимается перед вызовом
    liftover.lift(), т.к. core/liftover.py::ChainLiftover работает с
    каноническими именами хромосом без префикса). Ключи результата —
    КАНОНИЧЕСКИЕ (без "chr"), совпадающие с форматом, который использует
    HRC/GRCh37 (REFERENCE_PANELS["hrc"]["chrom_prefix"] == "") и,
    соответственно, скелет трафарета/template/assembler.py.

    Возвращает (новый_словарь, число_отброшенных) — позиции, которые
    лифтовер не смог перенести (нет chain-блока, gap между блоками
    выравнивания), не попадают в новый словарь и учитываются во втором
    элементе кортежа. Никогда не бросает исключение сама по себе —
    единственный источник исключений здесь — сам liftover.lift(), у
    которого нет исключений в штатном пути (см. докстринг
    core/liftover.py::ChainLiftover.lift() — при неудаче возвращает None,
    не бросает).

    Если на одну и ту же целевую (chrom, pos) в GRCh37 отобразилось
    несколько исходных позиций (маловероятно для точечных SNP, но
    теоретически возможно на стыке chain-блоков) — последнее совпадение
    побеждает, коллизии логируются предупреждением, но не прерывают
    работу.
    """
    result: dict[str, str] = {}
    dropped = 0
    collisions = 0
    for key, value in genotypes.items():
        try:
            chrom_raw, pos_str = key.rsplit("_", 1)
            pos = int(pos_str)
        except (ValueError, TypeError):
            logger.warning(
                "liftback_imputed_genotypes: не удалось разобрать ключ %r "
                "(ожидался формат '<chrom>_<pos>') — пропускаю.", key,
            )
            dropped += 1
            continue

        chrom = chrom_raw[3:] if chrom_raw.lower().startswith("chr") else chrom_raw
        lifted = liftover.lift(chrom, pos)
        if lifted is None:
            dropped += 1
            continue

        new_chrom, new_pos = lifted
        new_key = f"{new_chrom}_{new_pos}"
        if new_key in result:
            collisions += 1
        result[new_key] = value

    if dropped:
        logger.warning(
            "liftback_imputed_genotypes: %d из %d позиций не перенесены "
            "обратно в GRCh37 (нет chain-блока/gap между блоками "
            "выравнивания) — исключены из результата.",
            dropped, len(genotypes),
        )
    if collisions:
        logger.warning(
            "liftback_imputed_genotypes: %d коллизий целевых позиций "
            "GRCh37 (несколько исходных позиций отобразились в одну и ту "
            "же) — сохранено последнее встреченное значение.",
            collisions,
        )

    return result, dropped


def _build_reference(
    args,
    source: str,
    panel: str = DEFAULT_PANEL,
    progress_cb: Optional[Callable[[float, str], None]] = None,
):
    """Создаёт ReferenceGenome только для источников, которым он нужен."""
    if not _needs_reference(source):
        if getattr(args, "reference", None) is not None:
            logger.info("Источник '%s' не использует референс — --reference игнорируется", source)
        return None
    ref_path = getattr(args, "reference", None)
    if ref_path is None:
        ref_path = ensure_reference_genome(progress_cb=progress_cb, panel=panel)
    return ReferenceGenome(Path(ref_path))


# ---------------------------------------------------------------------------
# Задача 4/6: надёжный пароль/распаковка результатов MIS.
# ---------------------------------------------------------------------------
def _sanitize_password(raw: str) -> str:
    return sanitize_password_text(raw)


def download_mis_results_smart(
    curl_command: str,
    results_dir: Path,
    password: str,
    on_file_error: Optional[Callable[[str, str], bool]] = None,
) -> None:
    """
    Скачивает и распаковывает результаты Michigan Imputation Server.
    При ошибке распаковки автоматически повторяет попытку с "очищенным"
    от невидимых unicode-символов паролем (частая причина "неверный пароль"
    при копировании из HTML-письма).

    on_file_error (промт "проверять уже скачанные файлы + предлагать
    повтор при ошибке"): прокидывается без изменений в
    MISAdapter.download_results() — см. её докстринг. None (по
    умолчанию, обратная совместимость для CLI/старых вызовов) — при
    ошибке скачивания файла повтор не предлагается, файл просто попадает
    в список неудавшихся и в конце бросается MISAdapterError со списком
    всех проблемных файлов; уже скачанные файлы при этом остаются на
    диске, и повторный запуск их не перекачивает.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # ⚠ Фикс: раньше sevenzip_path не передавался вовсе, из-за чего
    # MISAdapter искал 7z.exe только в PATH и захардкоженных путях
    # (C:\Program Files\7-Zip\...), полностью игнорируя папку бинарников
    # (--bin-dir/HTSLIB.bin_dir), куда пользователь обычно кладёт 7z.exe
    # вместе с bcftools.exe/tabix.exe. В этом случае распаковка молча
    # откатывалась на встроенный zipfile, который НЕ поддерживает
    # AES-256 (см. core/archive_utils.py) — архивы результатов MIS
    # зашифрованы именно им, поэтому распаковка могла падать или
    # создавать битые файлы. find_7z() внутри MISAdapter всё равно
    # проверяет explicit_path через .is_file() и откатывается на PATH/
    # захардкоженные пути сам, если файла там нет — так что передача
    # кандидата из bin_dir безопасна даже если 7z.exe там не лежит.
    sevenzip_candidate = (
        str(HTSLIB.bin_dir / ("7z.exe" if IS_WINDOWS else "7z"))
        if HTSLIB.bin_dir else None
    )

    adapter = MISAdapter(
        upload_dir=PROJECT_ROOT / "output" / "upload",
        results_dir=results_dir,
        bcftools_path=HTSLIB.bcftools_path,
        sevenzip_path=sevenzip_candidate,
    )

    zip_paths = adapter.download_results(curl_command, on_file_error=on_file_error)
    print(f"✓ Скачано архивов: {len(zip_paths)}")

    try:
        adapter.extract_all_results(zip_paths, password)
        return
    except MISAdapterError as e:
        sanitized = _sanitize_password(password)
        if sanitized == password:
            raise
        logger.warning(
            "Распаковка с исходным паролем не удалась (%s). "
            "Пробую ещё раз с очищенным от пробелов/невидимых символов паролем.", e,
        )
        adapter.extract_all_results(zip_paths, sanitized)


# ---------------------------------------------------------------------------
# Вспомогательные bcftools-обёртки для стадий 1-6
# ---------------------------------------------------------------------------
def _run_bcftools(args: list[str]) -> None:
    if not HTSLIB.has_bcftools:
        raise RuntimeError("bcftools не найден. Укажите папку с бинарниками (--bin-dir).")
    result = subprocess.run(
        [HTSLIB.bcftools_path, *args], capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bcftools {' '.join(args)} завершился с ошибкой:\n{result.stderr}")


def _index_vcf(vcf_path: Path) -> None:
    """tabix -p vcf <файл> — индексация BGZF-сжатого VCF."""
    if not HTSLIB.has_tabix:
        raise RuntimeError("tabix не найден. Укажите папку с бинарниками (--bin-dir).")
    subprocess.run(
        [HTSLIB.tabix_path, "-p", "vcf", "-f", str(vcf_path)],
        check=True, capture_output=True,
    )
    logger.info("Индекс построен: %s.tbi", vcf_path)


# ---------------------------------------------------------------------------
# Промт "проверка 'донор не пустой' после скачивания/фильтрации".
#
# Подтверждено реальным прогоном: kgp_sub_*.vcf.gz для panel="topmed"
# могут получиться практически пустыми (корректный заголовок + список
# образцов, 0 строк с вариантами) из-за нестабильного HTTPS-соединения
# при скачивании/удалённой фильтрации в download_donors.py (например,
# активный VPN/прокси мешает стабильности соединения). Сам
# download_donors.py теперь отбраковывает такие хромосомы как неудавшиеся
# (см. _count_vcf_records() там же) — но это защищает только СВЕЖЕЕ
# скачивание. Кэш, который уже был один раз ошибочно записан на диск ДО
# этого фикса (chip_signature.txt при этом успешно совпадает), иначе
# продолжал бы молча приниматься как валидный check_donor_cache() на
# каждом следующем запуске.
# ---------------------------------------------------------------------------
def _count_vcf_records(vcf_path: Path) -> int:
    """
    Считает число строк с вариантами в VCF (`bcftools view -H | количество
    строк`). Возвращает -1, если проверку выполнить не удалось (bcftools
    не найден/упал) — это трактуется как "проверка неубедительна", а не
    как "файл пуст", чтобы не отбраковывать иначе валидный кэш только
    из-за временной недоступности bcftools в момент самой проверки.
    """
    if not HTSLIB.has_bcftools:
        return -1
    try:
        result = subprocess.run(
            [HTSLIB.bcftools_path, "view", "-H", str(vcf_path)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            return -1
        return len([l for l in result.stdout.splitlines() if l.strip()])
    except Exception:
        return -1


def _donor_sample_order(vcf_path: Path) -> list[str]:
    """Имена образцов из строки #CHROM донорского VCF, в порядке колонок."""
    # errors="replace" — bcftools дописывает в заголовок команду вызова с
    # полным путём к файлу; на Windows при не-ASCII символах в пути эта
    # строка может оказаться не в UTF-8 (та же причина, что и в
    # core/pure_python_core.py). Имена образцов всегда чистый ASCII.
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\r\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


def _align_donor_sample_order(donor_vcfs: list[Path]) -> None:
    """
    Приводит порядок колонок образцов во всех донорских файлах к порядку
    первого из них.

    Зачем (промт "Покрытие X-хромосомы", найдено живым прогоном):
    `bcftools concat` требует не только одинакового НАБОРА образцов во
    всех файлах, но и одинакового ПОРЯДКА колонок, иначе падает с
    "Different sample names in <файл>. Perhaps bcftools merge is what you
    are looking for?". Релиз 1000 Genomes phase3 перечисляет образцы в
    chrX в другом порядке, чем в аутосомах, и этот порядок доживает до
    kgp_sub_X.vcf.gz — набор образцов при этом идентичен (проверено:
    те же 20 из 20), различается только их расположение по колонкам.

    Переупорядочивание идемпотентно и делается ОДИН РАЗ: результат
    записывается обратно в кэш доноров, поэтому на следующих запусках
    порядок уже совпадает и функция ничего не делает. Файл с ДРУГИМ
    набором образцов (а не только порядком) не трогается — это признак
    испорченного/чужого кэша, и пусть об этом честно скажет сам
    bcftools concat, а не тихо "починит" эта функция.
    """
    if len(donor_vcfs) < 2:
        return
    reference_order = _donor_sample_order(donor_vcfs[0])
    if not reference_order:
        return
    reference_set = set(reference_order)

    for donor in donor_vcfs[1:]:
        order = _donor_sample_order(donor)
        if order == reference_order or set(order) != reference_set:
            continue
        logger.info(
            "Донор %s перечисляет образцов в другом порядке, чем %s — "
            "переупорядочиваю (иначе bcftools concat откажется объединять)",
            donor.name, donor_vcfs[0].name,
        )
        order_file = donor.parent / f"_sample_order_{donor.stem}.txt"
        reordered = donor.parent / f"{donor.name}.reordered.tmp.vcf.gz"
        try:
            with order_file.open("w", encoding="utf-8", newline="\n") as f:
                for name in reference_order:
                    f.write(name + "\n")
            _run_bcftools([
                "view", "-S", str(order_file), str(donor),
                "-Oz", "-o", str(reordered),
            ])
            reordered.replace(donor)
            _index_vcf(donor)
            logger.info("Порядок образцов в %s приведён к общему", donor.name)
        finally:
            order_file.unlink(missing_ok=True)
            reordered.unlink(missing_ok=True)


def _concat_donors(donor_vcfs: list[Path], output_vcf: Path) -> Path:
    """bcftools concat kgp_sub_{1..22,X}.vcf.gz -Oz -o kgp_all.vcf.gz + индекс."""
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    _align_donor_sample_order(donor_vcfs)
    _run_bcftools(["concat", *[str(p) for p in donor_vcfs], "-Oz", "-o", str(output_vcf)])
    _index_vcf(output_vcf)
    logger.info("Доноры объединены: %d файлов -> %s", len(donor_vcfs), output_vcf)
    return output_vcf


def _merge_with_donors_bcftools(sample_vcf: Path, kgp_all: Path, output_vcf: Path) -> Path:
    """
    bcftools merge sample.vcf.gz kgp_all.vcf.gz -Oz -o batch_merged.vcf.gz + индекс.

    Задача D: флаг -0/--missing-to-ref убран безусловно. В строгом режиме
    (сигнатура по QC-варианты конкретного человека) позиции донора и
    сэмпла и так совпадают 1-в-1, так что удаление -0 ничего не меняет.
    В широком режиме (chip_signature_broad) донорская панель шире, чем
    sample_vcf, и -0 сфабриковал бы ложные 0/0 на позициях, отброшенных
    у этого человека по QC (self-complementary), которые почти всегда
    на самом деле гетерозиготны. Без -0 такие позиции у сэмпла честно
    остаются './.' (missing) в объединённом VCF.
    """
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run_bcftools([
        "merge", "--force-samples",
        str(sample_vcf), str(kgp_all),
        "-Oz", "-o", str(output_vcf),
    ])
    _index_vcf(output_vcf)
    logger.info("Merge sample+доноры готов: %s", output_vcf)
    return output_vcf


# ---------------------------------------------------------------------------
# Задача A/B + Шаг 1 промта TopMed: единая проверка и (только на этапе
# скачивания!) запись сигнатуры чипа доноров — общая логика для CLI и GUI,
# теперь раздельно ещё и по референсной панели.
# ---------------------------------------------------------------------------
def _donor_source_dir(
    source: str,
    donors_root: Path = DONORS_DIR,
    panel: str = DEFAULT_PANEL,
) -> Path:
    """
    donors_root/<source>/<panel>/ — раздельное хранение доноров и по
    источнику (ftdna/myheritage/vcf, Задача B), и по референсной панели
    (hrc/topmed, Шаг 1 промта TopMed).

    ⚠ Изменение путей на диске по сравнению с предыдущей версией: раньше
    доноры лежали в donors/<source>/, теперь — в donors/<source>/hrc/ по
    умолчанию (panel="hrc" везде, где вызывающий код ещё не прокидывает
    panel явно). Если у вас уже есть скачанные доноры по старому пути
    donors/<source>/kgp_sub_*.vcf.gz и chip_signature.txt — перенесите их
    вручную в donors/<source>/hrc/, иначе они будут не видны
    check_donor_cache() и попросят перекачку. Это сознательный компромисс
    ради того, чтобы кэш HRC и TopMed никогда не могли перепутаться
    (по аналогии с тем, как раньше была устранена путаница кэша между
    разными source — см. докстринг Задачи A/B выше).
    """
    panel_key = panel or DEFAULT_PANEL
    return Path(donors_root) / source / panel_key


# ---------------------------------------------------------------------------
# Промт "Доноры для VCF-источника: понятная отмена + общий кэш сырых
# хромосом", Шаг 4/main.py-часть: единая точка правды для пути к общему
# кэшу ЕЩЁ НЕ ОТФИЛЬТРОВАННЫХ полных хромосом 1000 Genomes.
#
# В отличие от _donor_source_dir() (раздельно по source И panel), этот
# кэш раздельный ТОЛЬКО по genome_build ("grch37"/"grch38") — потому что
# полная, нефильтрованная хромосома 1000 Genomes совершенно одинакова
# независимо от того, какой источник (ftdna/myheritage/vcf) или какой
# конкретно чип её потом фильтрует локальным bcftools; отличаются только
# ПОЗИЦИИ фильтрации, а не сами скачиваемые данные. HRC и TopMed
# используют разные сборки генома (GRCh37 vs GRCh38 — разные релизы
# 1000 Genomes), поэтому раздельность по genome_build обязательна, а вот
# раздельность по source/panel-ключу (а не по самой сборке) была бы
# неоправданной — именно она и вызывала повторное скачивание одних и тех
# же ГБ при переключении источника, которое устраняет этот кэш.
# ---------------------------------------------------------------------------
def raw_chromosome_cache_dir(donors_root: Path = DONORS_DIR, panel: str = DEFAULT_PANEL) -> Path:
    """
    donors_root/_raw_chromosomes/<genome_build>/ — общий на все
    источники/чипы этой референсной сборки. По умолчанию НЕ используется
    нигде автоматически: вызывающий код (GUI-чекбокс "Хранить сырые
    хромосомы...", CLI-флаг --raw-chromosome-cache) сам решает, включать
    ли его, и передаёт результат в
    download_donors.download_donors_for_chip(..., raw_cache_dir=...).

    ⚠ Место на диске: 22 полные хромосомы 1000 Genomes phase3 весят
    суммарно порядка 15-20 ГБ (самая крупная, chr1, — около 1.1 ГБ в
    сжатом .vcf.gz виде; в сумме по всем 22 — по разным оценкам от ~13 до
    ~20+ ГБ в зависимости от того, какой конкретно релиз/суффикс
    зеркала). Это ДОПОЛНИТЕЛЬНО к месту, которое уже занимают
    отфильтрованные doнорские файлы каждого source/panel
    (donors/<source>/<panel>/, обычно десятки-сотни МБ, они на порядки
    меньше, потому что содержат только позиции конкретного чипа).
    """
    return Path(donors_root) / "_raw_chromosomes" / _panel_config(panel)["genome_build"]


def _save_chip_signature(
    chip_signature: str,
    source: str,
    donors_root: Path = DONORS_DIR,
    panel: str = DEFAULT_PANEL,
) -> Path:
    """
    Сохраняет сигнатуру чипа в donors_root/<source>/<panel>/chip_signature.txt.

    ВАЖНО (Задача A): эта функция больше НЕ вызывается из main()/GUI на
    этапе парсинга — только из download_donors.py, и только сразу после
    успешной свежей загрузки/фильтрации всех 22 донорских хромосом.
    Вызов на этапе парсинга делал сравнение сигнатур бессмысленным
    (сигнатура сравнивалась сама с собой уже ПОСЛЕ перезаписи) — именно
    так тихо принимался кэш доноров от другого чипа.
    """
    donors_dir = _donor_source_dir(source, donors_root, panel)
    donors_dir.mkdir(parents=True, exist_ok=True)
    sig_file = donors_dir / "chip_signature.txt"
    sig_file.write_text(chip_signature, encoding="utf-8")
    return sig_file


# ---------------------------------------------------------------------------
# Задача D: выбор режима сигнатуры (строгий по умолчанию / широкий —
# переиспользование доноров между разными людьми на одном чипе).
# ---------------------------------------------------------------------------
def _resolve_chip_signature_mode(
    result, source: str, reuse_donors_across_people: bool = False,
) -> tuple[str, Callable[[Path, object], Path]]:
    """
    Возвращает (signature, save_position_cache_fn) для дальнейшего
    использования в save_pos_fn(...)/check_donor_cache(...)/
    download_donors_for_chip(...).

    reuse_donors_across_people=False (по умолчанию): строгий режим —
        result.chip_signature (как раньше, поведение не меняется).
    reuse_donors_across_people=True: широкий режим —
        result.chip_signature_broad, если источник его поддерживает
        (ftdna/myheritage). Для source='vcf' широкого режима нет
        (chip_signature_broad там не считается) — тихо откатываемся на
        обычную сигнатуру с предупреждением в лог, запуск не падает.
    """
    save_broad_fn = SOURCES[source].get("save_position_cache_broad")
    if not reuse_donors_across_people or save_broad_fn is None:
        if reuse_donors_across_people and save_broad_fn is None:
            logger.warning(
                "⚠ Переиспользование доноров между людьми запрошено, но источник "
                "'%s' его не поддерживает (нет person-specific QC, широкая и "
                "строгая сигнатуры совпадают) — использую обычный режим.",
                source,
            )
        return result.chip_signature, SOURCES[source]["save_position_cache"]

    if not getattr(result, "chip_signature_broad", ""):
        logger.warning(
            "⚠ Переиспользование доноров запрошено, но chip_signature_broad "
            "пуста для источника '%s' — использую обычный режим.", source,
        )
        return result.chip_signature, SOURCES[source]["save_position_cache"]

    logger.info(
        "ℹ Режим переиспользования доноров включён: signature_broad=%s",
        result.chip_signature_broad,
    )
    return result.chip_signature_broad, save_broad_fn


def _warn_if_legacy_flat_donors(donors_root: Path = DONORS_DIR) -> None:
    """Задача B (обратная совместимость): если в donors/ напрямую лежат
    старые "плоские" kgp_sub_*.vcf.gz — не трогаем их, только предупреждаем."""
    flat_marker = Path(donors_root) / "kgp_sub_1.vcf.gz"
    if flat_marker.exists():
        logger.info(
            "ℹ Обнаружены доноры старого (плоского) формата в %s — они оставлены "
            "как есть, новые запуски используют %s/<source>/<panel>/.",
            donors_root, donors_root,
        )


def check_donor_cache(
    chip_signature: str,
    source: str,
    donors_root: Path = DONORS_DIR,
    panel: str = DEFAULT_PANEL,
) -> list[Path]:
    """
    Единая точка правды для CLI (main()) и GUI (_check_donors): проверяет,
    что кэш доноров под donors_root/<source>/<panel>/ полон (23 файла: 1-22 + X) и
    что его chip_signature.txt совпадает с chip_signature текущего запуска.

    Возвращает список путей kgp_sub_{1..22,X}.vcf.gz при успехе.
    Бросает RuntimeError с понятной инструкцией (--source, --panel,
    --donors-subdir) при отсутствии файлов или несовпадении сигнатуры —
    что и требуется по сценарию "сначала FTDNA, потом MyHeritage без
    явного download_donors" (а теперь ещё и "сначала HRC, потом TopMed
    для того же файла"): раньше он молча ломался (union из двух чипов/
    панелей), теперь падает явно.
    """
    _warn_if_legacy_flat_donors(donors_root)

    donors_dir = _donor_source_dir(source, donors_root, panel)
    sig_file = donors_dir / "chip_signature.txt"

    donor_vcfs: list[Path] = []
    missing: list[Path] = []
    # UPLOAD_CHROMS = 1..22 + X: тот же список, что уходит на сервер
    # импутации (core/pure_python_core.py) и что качает download_donors.py
    # (DONOR_CHROMS). Здесь берём его из core, потому что main.py
    # сознательно не импортирует download_donors.
    for chrom in UPLOAD_CHROMS:
        f = donors_dir / f"kgp_sub_{chrom}.vcf.gz"
        if f.exists():
            donor_vcfs.append(f)
        else:
            missing.append(f)

    panel_key = panel or DEFAULT_PANEL
    # ⚠ Фикс попутно с промтом "TopMed/HRC" п.3: download_donors.py CLI
    # никогда не имел флага --panel (только --donors-subdir/--genome-build) —
    # эта подсказка годами предлагала пользователю невалидную команду.
    # Теперь используется --genome-build (см. download_donors.py::main()),
    # значение берётся из REFERENCE_PANELS[panel_key]["genome_build"] — без
    # него доноры для TopMed молча тянулись бы из GRCh37-релиза 1000
    # Genomes (см. download_donors.py::DEFAULT_GENOME_BUILD).
    genome_build = _panel_config(panel_key)["genome_build"]
    download_cmd = (
        f"  python download_donors.py --source {source} "
        f"--genome-build {genome_build} "
        f"--donors-subdir {source}/{panel_key} "
        f"--positions-json \"<путь к donors/{source}/{panel_key}/<signature>.positions.json>\" "
        f"--output-dir donors --bin-dir <папка с bcftools>"
        f"\n  (необязательно: добавьте --raw-cache-dir "
        f"{raw_chromosome_cache_dir(donors_root, panel_key)} — общий кэш "
        f"нефильтрованных хромосом 1000 Genomes, переиспользуемый между "
        f"источниками/чипами этой сборки, экономит трафик ценой места на "
        f"диске)"
    )

    if missing or not sig_file.exists():
        missing_names = ", ".join(f.name for f in missing[:5])
        if len(missing) > 5:
            missing_names += f" и ещё {len(missing) - 5}"
        # Отдельная подсказка для самого частого случая после добавления
        # X в пайплайн: кэш доноров, скачанный прежней версией, полон по
        # аутосомам, и не хватает ровно kgp_sub_X.vcf.gz. Перекачивать
        # всё заново не нужно — download_donors пропускает уже готовые
        # хромосомы и дотянет только недостающие.
        only_x_missing = (
            len(missing) == 1 and missing[0].name == "kgp_sub_X.vcf.gz"
            and sig_file.exists()
        )
        hint = (
            "\nЭто ожидаемо для кэша доноров, скачанного до появления "
            "поддержки X-хромосомы: не хватает только её. Уже скачанные "
            "хромосомы повторно не качаются — дотянется одна X."
            if only_x_missing else ""
        )
        what_missing = missing_names or "нет файла сигнатуры"
        raise RuntimeError(
            f"Донорские файлы для источника '{source}' (панель '{panel_key}') "
            f"отсутствуют или ещё не скачаны ({len(missing)} из "
            f"{len(UPLOAD_CHROMS)} отсутствует: {what_missing}, "
            f"папка: {donors_dir}).{hint}\nСначала запустите:\n{download_cmd}"
        )

    cached = sig_file.read_text(encoding="utf-8").strip()
    if cached != chip_signature:
        logger.warning(
            "⚠ Кэш доноров НЕ соответствует текущему чипу (source=%s, panel=%s): "
            "ожидалась сигнатура %s, на диске %s.",
            source, panel_key, chip_signature, cached,
        )
        raise RuntimeError(
            f"⚠ Кэш доноров НЕ соответствует текущему чипу (source={source}, "
            f"panel={panel_key}): ожидалась сигнатура {chip_signature}, на диске "
            f"{cached}.\nНужна повторная загрузка доноров. Запустите:\n{download_cmd}"
        )

    # Промт "проверка 'донор не пустой' после скачивания/фильтрации":
    # chip_signature.txt мог быть записан ДО того, как в
    # download_donors.py появилась защита от пустых доноров (см.
    # _count_vcf_records() там же) — то есть сигнатура совпадает, но сами
    # kgp_sub_*.vcf.gz могут быть пустыми (корректный заголовок, 0
    # записей) из-за оборванного в своё время сетевого соединения
    # (VPN/прокси и т.п.). Полная проверка всех файлов на каждый
    # запуск была бы слишком дорогой — проверяем только пару файлов
    # (первый и последний в списке, обычно chr1 и chr22), этого
    # достаточно, чтобы поймать типичный случай "весь кэш скачался
    # пустым при массовом сетевом сбое", не превращая обычную проверку
    # кэша в полный повторный подсчёт по всем хромосомам.
    sample_files = donor_vcfs[:1] + donor_vcfs[-1:] if len(donor_vcfs) > 1 else donor_vcfs[:1]
    empty_samples: list[Path] = []
    for f in sample_files:
        count = _count_vcf_records(f)
        if count == 0:
            empty_samples.append(f)
    if empty_samples:
        names = ", ".join(p.name for p in empty_samples)
        raise RuntimeError(
            f"⚠ Кэш доноров (source={source}, panel={panel_key}) в {donors_dir} "
            f"числится актуальным по сигнатуре, но следующие файлы пусты "
            f"(0 записей после фильтрации): {names}. Вероятно, при прошлом "
            f"скачивании оборвалось сетевое соединение (активный VPN/прокси "
            f"нередко тому виной). Удалите содержимое {donors_dir} и "
            f"перекачайте доноров заново:\n{download_cmd}"
        )

    logger.info(
        "✓ Кэш доноров (%s/%s) актуален: signature=%s", source, panel_key, chip_signature,
    )
    return donor_vcfs


# ---------------------------------------------------------------------------
# Именованные/нумерованные папки запуска (промт "Именованные папки
# запуска").
#
# Раньше все запуски писали в одну общую output/ — второй человек, начатый
# без ожидания писем MIS по первому, тихо перезаписывал файлы первого
# (sample.vcf.gz, upload/, parse_result.pkl и т.д.), а перезапуск
# приложения между Этапом 1-6 и Этапом 7 (письмо может прийти через
# несколько часов) терял привязку "какие файлы к какому человеку".
#
# Теперь каждый запуск получает СВОЮ подпапку output/runs/<run_name>/ —
# все файлы этого конкретного запуска живут внутри неё. donors/<source>/
# <panel>/ и референсные .fasta остаются ОБЩИМИ на все запуски (см.
# докстринг промта) — они дорогие для перекачки и по конструкции
# переиспользуемы (chip_signature.txt, Задачи A/B/D), поэтому дублировать
# их на каждый запуск не нужно.
# ---------------------------------------------------------------------------
RUNS_SUBDIR_NAME = "runs"
_INVALID_RUN_NAME_CHARS = re.compile(r'[\\/:*?"<>|]')


def validate_run_name(name: str) -> str:
    """
    Проверяет и нормализует имя запуска (пользовательский ввод в GUI или
    --run-name в CLI) — должно быть валидным именем папки на Windows:
    не пустое, без ведущих/конечных пробелов, без символов
    \\ / : * ? " < > | и не "." / "..".
    """
    cleaned = (name or "").strip()
    if not cleaned:
        raise ValueError("Имя запуска не может быть пустым")
    if _INVALID_RUN_NAME_CHARS.search(cleaned):
        raise ValueError(
            f"Имя запуска {cleaned!r} содержит недопустимые для имени "
            f"папки символы (\\ / : * ? \" < > |) — используйте буквы, "
            f"цифры, пробелы, дефисы и подчёркивания."
        )
    if cleaned in (".", ".."):
        raise ValueError(f"Имя запуска не может быть {cleaned!r}")
    return cleaned


def _next_run_name(runs_root: Path) -> str:
    """
    Следующий свободный числовой run_name внутри runs_root — по
    умолчанию "1", если папок ещё нет, иначе (макс. существующее
    числовое имя) + 1. Нечисловые имена запусков (пользователь мог
    переименовать запуск) в подсчёте max игнорируются, а не ломают его.
    """
    runs_root = Path(runs_root)
    runs_root.mkdir(parents=True, exist_ok=True)
    existing_nums = [
        int(p.name) for p in runs_root.iterdir()
        if p.is_dir() and p.name.isdigit()
    ]
    return str(max(existing_nums, default=0) + 1)


def list_runs(output_root: Path) -> list[Path]:
    """
    Список существующих папок запуска (output_root/runs/*), новые
    первыми (по времени последнего изменения) — используется GUI для
    списка "История запусков".
    """
    runs_root = Path(output_root) / RUNS_SUBDIR_NAME
    if not runs_root.exists():
        return []
    dirs = [p for p in runs_root.iterdir() if p.is_dir()]
    return sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)


def resolve_run_dir(
    output_root: Path,
    run_name: Optional[str] = None,
    must_exist: bool = False,
) -> tuple[Path, str]:
    """
    Определяет (и, если нужно, создаёт) папку конкретного запуска
    output_root/runs/<run_name>/.

    run_name=None -> автонумерация (_next_run_name()).

    must_exist=False (по умолчанию, соответствует Этапу 1 — началу
        нового запуска): папка НЕ должна уже существовать — коллизия
        имён теперь явная ошибка, а не тихая перезапись чужого запуска,
        как было раньше с общей output/.
    must_exist=True (Этап 7/"Продолжить" — привязка к уже выполненному
        Этапу 1-6, в т.ч. после перезапуска приложения): папка ДОЛЖНА
        уже существовать, иначе явная ошибка вместо создания пустой
        папки с тем же именем.
    """
    runs_root = Path(output_root) / RUNS_SUBDIR_NAME
    runs_root.mkdir(parents=True, exist_ok=True)

    if run_name is None:
        run_name = _next_run_name(runs_root)
    else:
        run_name = validate_run_name(run_name)

    run_dir = runs_root / run_name

    if must_exist:
        if not run_dir.is_dir():
            raise RuntimeError(
                f"Папка запуска '{run_name}' не найдена в {runs_root} — "
                f"проверьте имя запуска (Этап 2/'Продолжить' должен "
                f"указывать на уже выполненный Этап 1-6)."
            )
    else:
        if run_dir.exists():
            raise RuntimeError(
                f"Папка запуска '{run_name}' уже существует ({run_dir}) — "
                f"выберите другое имя запуска, чтобы не затереть "
                f"предыдущий запуск с тем же именем."
            )
        run_dir.mkdir(parents=True, exist_ok=True)

    return run_dir, run_name


def _warn_if_legacy_flat_output(output_root: Path) -> None:
    """
    Промт "Именованные папки запуска", Шаг 6 (обратная совместимость):
    если в output/ напрямую (не в output/runs/<name>/) лежат файлы
    старого формата (sample.vcf.gz и т.п. — как писал пайплайн до этого
    промта) — НЕ трогаем и НЕ переносим их автоматически (риск задеть
    недообработанный запуск, для которого ещё не был выполнен Этап 7),
    только предупреждаем один раз в лог/консоль.
    """
    legacy_marker = Path(output_root) / "sample.vcf.gz"
    if legacy_marker.exists():
        logger.info(
            "ℹ Обнаружены файлы старого формата прямо в %s (sample.vcf.gz и "
            "т.п.) — новая версия их не использует и не трогает; если Этап 7 "
            "по ним уже выполнен (или не нужен), можно удалить вручную. "
            "Новые запуски используют %s/%s/<имя_запуска>/.",
            output_root, output_root, RUNS_SUBDIR_NAME,
        )


def save_run_info(run_dir: Path, **fields) -> Path:
    """
    Записывает/дополняет run_dir/run_info.json — метаданные запуска для
    списка "История запусков" в GUI (например: "3 — ftdna, HRC, call
    rate 96.4%, 22.08.2026 16:40" вместо голого имени папки).

    Обновление частичное и накопительное: существующие поля в файле
    сохраняются, новые/переданные — накладываются поверх. Так Этап 1-6
    пишет source/panel/chip_signature при старте, а Этап 7 того же
    запуска позже дописывает call_rate/finished_at, не затирая
    остальное. Поля со значением None игнорируются (не перезаписывают
    уже сохранённое значение).
    """
    run_dir = Path(run_dir)
    info_path = run_dir / "run_info.json"
    data: dict = {}
    if info_path.exists():
        try:
            data = json.loads(info_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            data = {}
    data.update({k: v for k, v in fields.items() if v is not None})
    info_path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    return info_path


def load_run_info(run_dir: Path) -> dict:
    """Читает run_dir/run_info.json, {} если файла нет или он повреждён —
    не бросает исключение (используется для отображения в GUI-списке)."""
    info_path = Path(run_dir) / "run_info.json"
    if not info_path.exists():
        return {}
    try:
        return json.loads(info_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def format_run_label(run_dir: Path) -> str:
    """
    Человекочитаемая подпись запуска для списка "История запусков",
    собранная из run_info.json: "<имя> — <источник>, <панель>[, call
    rate X%][, дата]". Если метаданных нет (старый/повреждённый запуск)
    — просто имя папки, без падения.
    """
    run_dir = Path(run_dir)
    info = load_run_info(run_dir)
    parts = [run_dir.name]
    if info.get("source"):
        parts.append(str(info["source"]))
    panel_cfg = REFERENCE_PANELS.get(info.get("panel", ""))
    if panel_cfg:
        parts.append(panel_cfg["display_name"])
    elif info.get("panel"):
        parts.append(str(info["panel"]))
    if info.get("call_rate") is not None:
        try:
            parts.append(f"call rate {float(info['call_rate']):.1f}%")
        except (TypeError, ValueError):
            pass
    ts = info.get("finished_at") or info.get("started_at")
    if ts:
        parts.append(str(ts))
    return " — ".join(parts)


def attach_run_log_handler(run_dir: Path) -> logging.Handler:
    """
    Добавляет к корневому логгеру FileHandler, дублирующий все
    logger.info/warning/error-сообщения (из main.py и остальных модулей
    проекта — они используют logging.getLogger(__name__) с
    propagate=True по умолчанию, поэтому долетают до root) в
    <run_dir>/run.log.

    ⚠ print()-сообщения logging НЕ перехватывает (это отдельный поток
    вывода) — их дублирование в run.log обеспечивает вызывающий код
    отдельно: _Tee в CLI main() и file_path у LogRedirector в
    gui/app.py.

    Возвращает handler, чтобы вызывающий код мог снять его
    (logging.getLogger().removeHandler(handler)) при завершении запуска
    или переключении на другой запуск в рамках одного долгоживущего
    GUI-процесса — иначе сообщения следующего запуска продолжали бы
    дублироваться и в лог предыдущего.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(run_dir / "run.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logging.getLogger().addHandler(handler)
    return handler


class _Tee(io.TextIOBase):
    """
    Дублирует запись в исходный поток (обычно sys.stdout/sys.stderr) и в
    файл — используется CLI main() для сохранения текстового вывода
    print() в <run_dir>/run.log (промт "Именованные папки запуска").
    logger.*-сообщения сюда не попадают (см. attach_run_log_handler) —
    Tee нужен только для дублирования print().
    """
    def __init__(self, original, file_obj):
        self._original = original
        self._file = file_obj

    def write(self, s):
        self._original.write(s)
        self._file.write(s)
        return len(s)

    def flush(self):
        self._original.flush()
        self._file.flush()


# ---------------------------------------------------------------------------
# Задача C: post-merge intersect — диагностический/защитный слой.
# ---------------------------------------------------------------------------
def _post_merge_intersect(
    merged_vcf: Path,
    donor_vcfs: list[Path],
    output_vcf: Path,
    bcftools_path: Optional[str] = None,
    kgp_all_vcf: Optional[Path] = None,
) -> tuple[Path, int, int]:
    """
    Оставляет в merged VCF только позиции, реально присутствующие в донорах
    текущего запуска (объединение позиций всех kgp_sub_*.vcf.gz).

    ⚠ ЭТО НЕ ДИАГНОСТИКА, А ОБЯЗАТЕЛЬНЫЙ ШАГ для не-HRC панелей (найдено
    живым прогоном на panel="topmed"): без него в отправку на MIS/BioData
    Catalyst попадают десятки/сотни тысяч позиций, физически отсутствующих
    в донорской подвыборке 1000 Genomes (сайт не полиморфен в этой
    подвыборке). Сервер режет такие позиции как "Invalid Alleles" (ALT="."
    без вариантной альтернативы, см. core/pure_python_core.py::_write_vcf_line)
    и как "SNPs call rate < 90%" (донорам ставится "./." — bcftools merge
    без -0/--missing-to-ref их не подставляет фиктивным 0/0, см. Задачу D,
    и это правильное поведение) — после чего НИ ОДИН chunk не проходит QC.
    Раньше эта функция считалась чисто диагностическим слоем ("должно быть
    no-op после Задач A/B") — предположение оказалось неверным для панелей,
    где донорская подвыборка не покрывает позиции чипа так же плотно, как
    HRC. Вызывающий код (main()/gui/app.py) теперь ПРЕРЫВАЕТ запуск при
    сбое этой функции, а не продолжает молча — см. их комментарии.

    kgp_all_vcf (фикс "Failed to read the regions"): если передан и
    существует — используется НАПРЯМУЮ как аргумент -R вместо ручной сборки
    текстового common_pos.txt. bcftools/htslib официально поддерживают
    tabix-индексированный VCF/BCF как источник regions — это убирает
    хрупкий текстовый парсер regidx, который на Windows дважды падал с
    "Failed to read the regions" (сначала из-за CRLF, потом повторно уже
    после фикса CRLF — по неустановленной до конца причине, вероятно
    связанной с масштабом текстового файла на конкретной сборке htslib).
    kgp_all_vcf — обычно output_dir/"kgp_all.vcf.gz", уже готовый
    BGZF+tabix-индексированный результат _concat_donors(), несущий ровно ту
    же информацию о позициях, что раньше собиралась вручную через
    bcftools query по каждому donor_vcf.

    Если kgp_all_vcf не передан/не существует, либо путь через него не
    сработал — используется старый текстовый путь (common_pos.txt) как
    fallback, для обратной совместимости и подстраховки. Логика сортировки
    текстового fallback (_sort_key/_chrom_sort_key) и newline="\n" при
    записи common_pos.txt сохранены без изменений — это уже подтверждённые
    фиксы (см. историю правок в комментариях старой версии функции), просто
    теперь это запасной, а не основной путь.
    """
    bcftools_path = bcftools_path or HTSLIB.bcftools_path
    if not bcftools_path:
        raise RuntimeError("bcftools не найден для post-merge intersect")

    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    def _sort_key(line: str):
        chrom, pos = line.split("\t")
        # Каноническая сортировка (_chrom_sort_key), а не лексикографическая
        # по строке chrom — для TopMed ("chr1".."chr22") лексикографический
        # порядок не совпадает с физическим порядком в tabix-индексе
        # merged_vcf и мог приводить к аварийному завершению bcftools.
        return (_chrom_sort_key(chrom), int(pos))

    def _build_text_regions() -> Path:
        common_pos = output_vcf.parent / "common_pos.txt"
        positions: set[str] = set()
        for donor_vcf in donor_vcfs:
            res = subprocess.run(
                [bcftools_path, "query", "-f", "%CHROM\t%POS\n", str(donor_vcf)],
                capture_output=True, text=True, check=True,
            )
            positions.update(line for line in res.stdout.splitlines() if line)
        # newline="\n" явно — без этого питоновский текстовый режим на
        # Windows транслирует "\n" в "\r\n" при записи, что htslib-парсер
        # regions-файла на некоторых сборках не переваривает.
        with common_pos.open("w", encoding="utf-8", newline="\n") as f:
            for line in sorted(positions, key=_sort_key):
                f.write(line + "\n")
        return common_pos

    def _count_records(vcf_path: Path) -> int:
        res = subprocess.run(
            [bcftools_path, "view", "-H", str(vcf_path)],
            capture_output=True, text=True,
        )
        if res.returncode != 0:
            raise PureCoreError(
                f"bcftools view -H ({vcf_path.name}) завершился с кодом "
                f"{res.returncode}:\n"
                f"{res.stderr.strip() or '(stderr пуст — вероятно, аварийное завершение процесса)'}"
            )
        return len(res.stdout.splitlines())

    before = _count_records(merged_vcf)

    # Удаляем возможный "хвост" от прерванного прошлого запуска —
    # залоченный или недописанный output_vcf с прошлой попытки на Windows
    # иногда даёт неинформативные коды возврата у процесса, который
    # пытается его перезаписать (антивирус держит файл открытым и т.п.).
    output_vcf.unlink(missing_ok=True)

    def _run_intersect(regions_arg: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [bcftools_path, "view", "-R", regions_arg, str(merged_vcf),
             "-Oz", "-o", str(output_vcf)],
            capture_output=True, text=True,
        )

    used_kgp_all = kgp_all_vcf is not None and Path(kgp_all_vcf).exists()
    res = None
    if used_kgp_all:
        res = _run_intersect(str(Path(kgp_all_vcf)))
        if res.returncode != 0:
            logger.warning(
                "⚠ Post-merge intersect через kgp_all.vcf.gz не удался "
                "(код %s) — пробую fallback через текстовый common_pos.txt.",
                res.returncode,
            )
            output_vcf.unlink(missing_ok=True)
            used_kgp_all = False

    if not used_kgp_all:
        common_pos = _build_text_regions()
        res = _run_intersect(str(common_pos))

    if res.returncode != 0:
        # ⚠ Раньше был check=True без сохранения stderr — реальная причина
        # от bcftools терялась, наружу долетал только бессмысленный код
        # возврата (напр. 0xFFFFFFFF на Windows).
        raise PureCoreError(
            f"bcftools view -R (post-merge intersect) завершился с кодом "
            f"{res.returncode}:\n"
            f"{res.stderr.strip() or '(stderr пуст — вероятно, аварийное завершение процесса, не обычная ошибка bcftools; проверьте антивирус/залоченные файлы)'}"
        )
    _index_vcf(output_vcf)

    after = _count_records(output_vcf)
    removed = before - after
    if removed > 0:
        logger.info(
            "✓ Post-merge intersect: %d → %d позиций (%d удалено — позиции, "
            "отсутствующие в донорской подвыборке, корректно отфильтрованы).",
            before, after, removed,
        )
    else:
        logger.warning(
            "⚠ Post-merge intersect: 0 удалено (%d → %d). Для panel != 'hrc' "
            "это может означать, что кэш доноров/regions-фильтр не сработал "
            "как ожидалось — стоит перепроверить итоговый call rate после "
            "Этапа 7.",
            before, after,
        )
    return output_vcf, before, after


def _normalize_vcf(vcf_path: Path, reference_fasta: Path, output_vcf: Path,
                    bcftools_path: Optional[str] = None) -> Path:
    """
    Опциональная нормализация multiallelic-сайтов (Задача C, необязательный
    чекбокс в GUI, по умолчанию выключен). НЕ входит в критический путь
    фикса Invalid alleles — отдельная оптимизация, а не фикс основного
    бага. `bcftools norm -f <reference> -m-both`.
    """
    bcftools_path = bcftools_path or HTSLIB.bcftools_path
    if not bcftools_path:
        raise RuntimeError("bcftools не найден для нормализации (bcftools norm)")
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [bcftools_path, "norm", "-f", str(reference_fasta), "-m-both",
         str(vcf_path), "-Oz", "-o", str(output_vcf)],
        check=True, capture_output=True,
    )
    _index_vcf(output_vcf)
    logger.info("Нормализация (bcftools norm -m-both) выполнена: %s", output_vcf)
    return output_vcf


# ---------------------------------------------------------------------------
# CLI: прогон полного пайплайна из терминала
# ---------------------------------------------------------------------------
def _parse_args():
    parser = argparse.ArgumentParser(
        description="Конвертер ДНК-файлов (FTDNA/MyHeritage/VCF) -> 23andMe формат для Генотек",
    )
    parser.add_argument("--source", choices=list(SOURCES.keys()), default="ftdna")
    parser.add_argument(
        "--panel", choices=list(REFERENCE_PANELS.keys()), default=DEFAULT_PANEL,
        help=(
            "Референсная панель импутации: 'hrc' — HRC r1.1 2016 (GRCh37/hg19, "
            "поведение по умолчанию, без изменений), 'topmed' — TOPMed r3 "
            "(GRCh38/hg38, включает лифтовер координат GRCh37 -> GRCh38 через "
            "core/liftover.py для источников ftdna/myheritage; source='vcf' "
            "лифтовер пока не поддерживает, см. _supports_liftover())."
        ),
    )
    parser.add_argument("--csv", required=True, type=Path, help="Исходный файл (ftdna.csv и т.п.)")
    parser.add_argument("--template", required=True, type=Path, help="Трафарет (template.txt)")
    parser.add_argument("--reference", type=Path, default=None,
                         help="Путь к .fasta референса выбранной панели (по умолчанию — автозагрузка)")
    parser.add_argument("--bin-dir", type=Path, default=None, help="Папка с bcftools/tabix/bgzip")
    parser.add_argument("--donors-dir", type=Path, default=DONORS_DIR,
                         help="Корень хранилища доноров (внутри — подпапки по --source и --panel)")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "output",
                         help="Корень для папок запуска (промт 'Именованные папки запуска'). "
                              "⚠ Breaking change: раньше файлы писались прямо сюда, теперь "
                              "реальная рабочая папка запуска — <--output-dir>/runs/<--run-name>/.")
    parser.add_argument("--run-name", type=str, default=None,
                         help="Имя папки запуска внутри <--output-dir>/runs/. По умолчанию — "
                              "следующий свободный номер (1, 2, 3, ...). На Этапе 7 укажите "
                              "ТО ЖЕ имя, что было напечатано на Этапе 1-6 этого человека — "
                              "иначе Этап 7 не найдёт parse_result.pkl/upload/ от нужного запуска.")
    parser.add_argument("--format", choices=["v3", "v5"], default="v3")
    parser.add_argument("--rsq-threshold", type=float, default=0.30)
    parser.add_argument("--post-merge-intersect", action="store_true", default=True,
                         help="Диагностический post-merge intersect (Задача C, включён по умолчанию)")
    parser.add_argument("--normalize", action="store_true", default=False,
                         help="Опциональная bcftools norm -m-both (Задача C, выключено по умолчанию)")
    parser.add_argument("--auto-detect-source", action="store_true", default=False,
                         help="Сверить --source с реальным форматом файла (Задача 2). "
                              "По умолчанию выключено. При несовпадении с уверенностью "
                              ">= 0.8 — предупреждение в stderr, запуск НЕ прерывается.")
    parser.add_argument("--reuse-donors-across-people", action="store_true", default=False,
                         help="Опционально (Задача D, выключено по умолчанию): считать "
                              "сигнатуру по ВСЕМ измеренным позициям чипа (до QC-отбраковки "
                              "конкретного человека), чтобы кэш доноров переиспользовался "
                              "между разными людьми на одном и том же чипе. Требует, чтобы "
                              "_merge_with_donors_bcftools() не использовал -0/--missing-to-ref "
                              "(это обеспечено безусловно в этой версии).")
    parser.add_argument("--raw-chromosome-cache", action="store_true", default=False,
                         help="Опционально (промт 'Доноры для VCF-источника: понятная "
                              "отмена + общий кэш сырых хромосом', выключено по умолчанию): "
                              "подсказывает путь к общему кэшу ЕЩЁ НЕ отфильтрованных "
                              "полных хромосом 1000 Genomes (donors/_raw_chromosomes/"
                              "<genome_build>/, см. raw_chromosome_cache_dir()) в "
                              "инструкции для download_donors.py --raw-cache-dir, чтобы "
                              "второй/третий source/чип не перекачивал заново те же самые "
                              "многогигабайтные файлы зеркал 1000 Genomes. ⚠ Занимает "
                              "дополнительно ~десятки ГБ на диске. Это main.py CLI сам не "
                              "скачивает доноров (см. --run-name/[3/7] в докстринге файла) "
                              "— флаг только влияет на текст подсказки и на run_info.json.")
    return parser.parse_args()


def main() -> None:
    global HTSLIB
    args = _parse_args()

    HTSLIB = HtslibTools(args.bin_dir)
    if not HTSLIB.has_bcftools:
        sys.exit("ОШИБКА: bcftools не найден. Укажите --bin-dir или добавьте в PATH.")

    # v13 (промт "Диагностика + устойчивая настройка CA-сертификатов"):
    # гарантируем CA-сертификаты для libcurl (bcftools) и предупреждаем о
    # возможном конфликте bin/curl.exe ДО того, как что-либо в пайплайне
    # (в первую очередь удалённая фильтрация доноров, Этап 3) попытается
    # обратиться к HTTPS через bcftools. Раньше это требовало ручной
    # настройки $env:CURL_CA_BUNDLE в каждой новой консоли — см. докстринг
    # core/network_utils.py. Не бросает исключение при неудаче (например,
    # нет интернета) — в этом случае удалённая фильтрация просто не
    # включится, и пайплайн продолжит работу через полное скачивание.
    ensure_network_ready(args.bin_dir)

    output_root = Path(args.output_dir)
    _warn_if_legacy_flat_output(output_root)
    try:
        output_dir, run_name = resolve_run_dir(output_root, args.run_name, must_exist=False)
    except RuntimeError as e:
        sys.exit(str(e))

    # промт "Именованные папки запуска": дублируем весь текстовый вывод
    # этого запуска в <run_dir>/run.log — logger.*-сообщения через
    # attach_run_log_handler() (root logger, ловит main.py и остальные
    # модули), print()-сообщения через простой _Tee поверх sys.stdout/
    # sys.stderr. CLI — одноразовый процесс, поэтому переключение
    # sys.stdout/stderr не восстанавливается явно (процесс завершится
    # сразу после main()); долгоживущий GUI делает то же самое иначе,
    # см. gui/app.py::LogRedirector.
    attach_run_log_handler(output_dir)
    _run_log_file = (output_dir / "run.log").open("a", encoding="utf-8")
    sys.stdout = _Tee(sys.stdout, _run_log_file)
    sys.stderr = _Tee(sys.stderr, _run_log_file)

    print(f"ℹ Папка запуска: {output_dir} (имя запуска: {run_name!r})")
    save_run_info(
        output_dir,
        run_name=run_name,
        started_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        source=args.source,
        panel=args.panel,
        csv_filename=Path(args.csv).name,
        format=args.format,
        rsq_threshold=args.rsq_threshold,
        normalize=args.normalize,
        reuse_donors_across_people=args.reuse_donors_across_people,
        raw_chromosome_cache=args.raw_chromosome_cache,
    )

    if not Path(args.csv).exists():
        sys.exit(f"ОШИБКА: файл с данными не найден: {args.csv}")

    panel_cfg = _panel_config(args.panel)
    print(f"ℹ Референсная панель: {panel_cfg['display_name']}")

    # Задача 2: детекция несоответствия источника (--source) и реального
    # формата файла. Вызывается ДО ensure_reference_genome() (может качать
    # несколько ГБ и не нужна источнику 'vcf'), поэтому вставлена прямо
    # здесь — сразу после проверки, что --csv указывает на существующий
    # файл, и до [0/7]. Проверка выполняется, только если пользователь
    # явно включил --auto-detect-source — CLI-поведение по умолчанию не
    # меняется.
    if args.auto_detect_source:
        detected_source, confidence = detect_source_from_file(Path(args.csv))
        if detected_source:
            logger.info(
                "ℹ Автодетект: файл похож на %s (уверенность %.2f), выбрано %s",
                detected_source, confidence, args.source,
            )
            if confidence >= 0.8 and detected_source != args.source:
                print(
                    f"⚠ Похоже, файл '{args.csv}' на самом деле в формате "
                    f"'{detected_source}' (уверенность {confidence:.2f}), а "
                    f"выбран --source {args.source}. Запуск ПРОДОЛЖАЕТСЯ с "
                    f"выбранным источником — если это ошибка, перезапустите "
                    f"с --source {detected_source}.",
                    file=sys.stderr,
                )
        else:
            logger.info("ℹ Автодетект: не удалось определить формат файла %s", args.csv)

    # --- Этап 0: приведение файла к оформлению 23andMe v3 --------------
    # Только для источников из _SOURCES_NEEDING_CONVERSION; для остальных
    # это no-op, возвращающий тот же путь. Стоит ДО проверки референса
    # (которая может качать гигабайты): если файл не того формата, узнать
    # об этом лучше сразу.
    csv_for_parsing = Path(args.csv)
    if args.source in _SOURCES_NEEDING_CONVERSION:
        print("[0a/7] Приведение исходного файла к оформлению 23andMe v3")
        try:
            # output_dir — папка ЭТОГО запуска (output/runs/<имя>), а не
            # общий корень output/: конвертированный файл принадлежит
            # конкретному запуску и лежит рядом с sample.vcf.gz и run.log.
            csv_for_parsing, conversion_stats = prepare_source_file(
                args.source, Path(args.csv), output_dir, args.template,
            )
        except AncestryConvertError as e:
            sys.exit(f"ОШИБКА: {e}")
        print(f"  {conversion_stats.summary()}")
        if not conversion_stats.skipped:
            save_run_info(output_dir,
                          converted_file=Path(conversion_stats.out_path).name)

    print("[0/7] Проверка референсного генома")
    reference = _build_reference(args, args.source, panel=args.panel)

    # Промт "HRC / TopMed", лифтовер координат: строится СРАЗУ после
    # референса (тот же этап "подготовка входных данных, до парсинга") —
    # _build_liftover() возвращает None для panel="hrc" (быстрый no-op) и
    # для panel="topmed" + source="vcf" (лифтовер для VCF-источника пока
    # не реализован, см. _supports_liftover()); в последнем случае явно
    # предупреждаем в лог, а не молча теряем корректность координат.
    print("[0b/7] Проверка chain-файла лифтовера (если панель этого требует)")
    liftover = _build_liftover(args.panel) if _supports_liftover(args.source) else None
    if args.panel != DEFAULT_PANEL and not _supports_liftover(args.source):
        logger.warning(
            "⚠ Источник '%s' пока не поддерживает лифтовер координат "
            "(см. main.py::_supports_liftover()) — если ваши координаты "
            "не в сборке %s, результат импутации будет некорректным. "
            "Поддерживается сейчас только для source in ('ftdna', 'myheritage').",
            args.source, panel_cfg["genome_build"].upper(),
        )

    print("[1/7] Парсинг исходных данных")
    parser_fn = SOURCES[args.source]["parser"]
    # csv_for_parsing — результат Этапа 0 (для 'ancestry') либо сам
    # args.csv (для всех остальных источников).
    if _supports_liftover(args.source):
        result = parser_fn(csv_for_parsing, reference, liftover=liftover)
    else:
        result = parser_fn(csv_for_parsing, reference)
    print(f"  Годных вариантов: {len(result.variants)}, сигнатура: {result.chip_signature}")
    if getattr(result, "lift_failed", 0):
        print(f"  ⚠ Не перенесено лифтовером на целевую сборку (lift_failed): {result.lift_failed}")

    # Задача B (доп. пункт): сохраняем позиции чипа тем же парсером, который
    # уже отработал — download_donors.py потом читает этот JSON вместо
    # повторного (и потенциально неверного для MyHeritage/VCF) парсинга CSV.
    # Задача D: выбор строгой/широкой сигнатуры — единая точка для CLI и GUI.
    signature, save_pos_fn = _resolve_chip_signature_mode(
        result, args.source, reuse_donors_across_people=args.reuse_donors_across_people,
    )
    positions_cache_dir = _donor_source_dir(args.source, args.donors_dir, panel=args.panel)
    positions_json = save_pos_fn(positions_cache_dir, result)
    print(f"  Позиции чипа сохранены: {positions_json} (режим: "
          f"{'широкий' if signature == getattr(result, 'chip_signature_broad', None) else 'строгий'}, "
          f"панель: {panel_cfg['display_name']})")
    save_run_info(
        output_dir,
        chip_signature=result.chip_signature,
        chip_signature_broad=getattr(result, "chip_signature_broad", None),
    )

    # ВАЖНО (Задача A): здесь больше НЕТ вызова _save_chip_signature —
    # сигнатура пишется только в download_donors.py после свежего скачивания.

    with (output_dir / "parse_result.pkl").open("wb") as f:
        pickle.dump(result, f)

    print("[2/7] Построение VCF")
    sample_vcf = output_dir / "sample.vcf.gz"
    # Промт "HRC / TopMed", п.3: chrom_prefix из конфигурации панели ("chr"
    # для topmed/GRCh38, "" для hrc/GRCh37) — координаты в result.variants
    # к этому моменту уже лифтованы в целевую сборку (см. liftover=...
    # выше), поэтому CHROM-запись VCF должна использовать те же имена
    # контигов, что и GRCh38-референс/доноры (kgp_sub_*.vcf.gz качаются
    # download_donors.py с тем же префиксом через GENOME_BUILD_CHROM_PREFIX)
    # — иначе bcftools merge/concat не совпадёт по CHROM вообще.
    # Промт "Покрытие X-хромосомы": пол определяется по гетерозиготности
    # nonPAR X ДО записи VCF — от него зависит плоидность мужского X,
    # которую Michigan Imputation Server проверяет отдельным Ploidy Check
    # (см. build_vcf(haploid_x=...)).
    male, x_het_pct, x_calls = infer_male_from_variants(
        result.variants, genome_build=panel_cfg["genome_build"],
    )
    if x_calls:
        print(f"  Пол по X: {'мужской' if male else 'женский'} "
              f"(гетерозиготность nonPAR X {x_het_pct:.2f}% по {x_calls} позициям)")
    else:
        print("  Пол по X не определён: в файле нет калиброванных позиций nonPAR X")
    save_run_info(output_dir, sex_by_x="male" if male else "female",
                  x_het_pct=round(x_het_pct, 3), x_nonpar_calls=x_calls)
    build_vcf(
        result, sample_vcf, sample_name="genotek", bgzip_path=HTSLIB.bgzip_path,
        chrom_prefix=panel_cfg["chrom_prefix"],
        haploid_x=male, genome_build=panel_cfg["genome_build"],
    )
    _index_vcf(sample_vcf)

    print("[3/7] Проверка кэша доноров")
    try:
        donor_vcfs = check_donor_cache(signature, args.source, args.donors_dir, panel=args.panel)
    except RuntimeError as e:
        sys.exit(str(e))

    print(f"[4/7] Объединение доноров ({len(donor_vcfs)} хромосом)")
    kgp_all = output_dir / "kgp_all.vcf.gz"
    _concat_donors(donor_vcfs, kgp_all)

    print("[5/7] Merge sample + donors")
    merged = output_dir / "batch_merged.vcf.gz"
    _merge_with_donors_bcftools(sample_vcf, kgp_all, merged)

    if args.normalize:
        if reference is None:
            logger.warning("--normalize запрошен, но у источника '%s' нет референса — пропускаю", args.source)
        else:
            print("[5b/7] Нормализация (bcftools norm -m-both, опционально)")
            normalized = output_dir / "batch_merged.norm.vcf.gz"
            # Задача C (исправление): раньше здесь был fallback
            # `Path(args.reference)`, который падал с TypeError при
            # автозагрузке референса (--reference не указан явно, значит
            # args.reference is None). ReferenceGenome теперь всегда хранит
            # свой fasta_path (adapters/ftdna_v3.py), так что fallback
            # больше не нужен.
            merged = _normalize_vcf(merged, reference.fasta_path, normalized, HTSLIB.bcftools_path)

    if args.post_merge_intersect:
        print("[5c/7] Post-merge intersect (обязательная фильтрация по донорским позициям)")
        checked = output_dir / "batch_merged.checked.vcf.gz"
        try:
            merged, before_n, after_n = _post_merge_intersect(
                merged, donor_vcfs, checked, HTSLIB.bcftools_path, kgp_all_vcf=kgp_all,
            )
            print(f"  {before_n} → {after_n} позиций ({before_n - after_n} удалено)")
        except (PureCoreError, subprocess.CalledProcessError) as e:
            # ⚠ Больше НЕ безобидная диагностика (см. докстринг
            # _post_merge_intersect): без фильтрации по донорским позициям
            # в отправку на сервер импутации попадут десятки/сотни тысяч
            # позиций, которых физически нет в донорской подвыборке 1000
            # Genomes — сервер (Michigan Imputation Server / BioData
            # Catalyst) гарантированно провалит QC. Подтверждено живым
            # прогоном на panel="topmed" (см. докстринг функции). Поэтому
            # сбой этого шага теперь фатален — останавливаем сборку ДО
            # разбивки по хромосомам и загрузки, а не продолжаем молча.
            sys.exit(
                f"ОШИБКА: post-merge intersect не удался: {e}\n\n"
                f"Это НЕ безобидная диагностика: без фильтрации по донорским "
                f"позициям в отправку на сервер импутации попадут десятки/сотни "
                f"тысяч позиций, которых физически нет в донорской подвыборке "
                f"1000 Genomes — сервер (Michigan Imputation Server / BioData "
                f"Catalyst) гарантированно провалит QC ('Invalid Alleles', "
                f"'SNPs call rate < 90%', 'No chunks passed the QC step').\n"
                f"Запуск остановлен ДО разбивки по хромосомам и загрузки, чтобы "
                f"не тратить время на заведомо провальное задание. Проверьте "
                f"bcftools/наличие kgp_all.vcf.gz и его .tbi и запустите заново."
            )

    print("[6/7] Разбивка по хромосомам для загрузки на MIS")
    upload_dir = output_dir / "upload"
    outputs = split_autosomes(
        merged, upload_dir, bgzip_path=HTSLIB.bgzip_path,
        chrom_prefix=panel_cfg["chrom_prefix"],
        genome_build=panel_cfg["genome_build"],
    )
    print(f"  Создано {len(outputs)} файлов в {upload_dir}")
    print("=" * 70)
    print(f"Загрузите все файлы из {upload_dir} на {panel_cfg['mis_upload_url']}")
    print(f"Reference Panel: {panel_cfg['mis_panel_value']}, Population: EUR")
    print(f"Папка этого запуска: {output_dir} (имя запуска: {run_name!r})")
    print("После получения письма запустите (с ТЕМ ЖЕ --run-name, иначе Этап 7 "
          "не найдёт файлы этого запуска):")
    print(f"  python main.py --source {args.source} --panel {args.panel} --csv \"{args.csv}\" "
          f"--template \"{args.template}\" --bin-dir \"{args.bin_dir or ''}\" "
          f"--run-name \"{run_name}\" "
          f"--finish-with-curl \"<curl-команда>\" --password \"<пароль>\"")
    print("=" * 70)


if __name__ == "__main__":
    main()