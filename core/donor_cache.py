"""
core/donor_cache.py

Сколько ДОНОРСКИХ ОБРАЗЦОВ лежит в кэше доноров — и совпадает ли это с
тем, что просит пользователь.

Зачем отдельный модуль (промт "галочка «все доступные EUR-доноры»
игнорируется"): актуальность кэша доноров проверялась ТОЛЬКО по
chip_signature.txt, то есть по чипу. Размер донорской подвыборки
(eur_sample_count) в сигнатуру не входит, поэтому кэш, собранный когда-то
на 20 образцах, считался актуальным и при включённой галочке "использовать
всех доступных EUR-доноров (~503)": этап скачивания доноров пропускался
целиком, и на MIS уезжали те же 20 образцов. А если скачивание всё-таки
запускалось, process_chromosome() видел готовые kgp_sub_*.vcf.gz и
печатал "chr1 уже готов" — файлы на 20 образцах оставались на месте.

Здесь собраны обе половины проверки, чтобы main.py (проверка кэша) и
download_donors.py (инвалидация кэша перед перекачкой) считали число
образцов одинаково и не расходились:

  * vcf_sample_count() — сколько образцов РЕАЛЬНО в донорском файле
    (колонки после FORMAT в строке #CHROM). Это источник правды: он не
    зависит от служебных eur{N}.txt, которые могли не дожить до текущего
    запуска, и работает для кэшей любой давности.
  * available_eur_count() — сколько EUR-образцов всего есть в панели
    1000 Genomes (файл панели лежит в той же папке доноров). Нужно,
    чтобы понять, что значит "все доступные" в числах.
  * eur_count_verdict() — сравнение "что на диске" с "что запрошено".
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import NamedTuple, Optional

#: Файл панели образцов 1000 Genomes, который download_donors.py кладёт
#: рядом с донорами (там же — SAMPLES_FILENAME). Дублируется здесь
#: намеренно: main.py сознательно не импортирует download_donors.py.
SAMPLES_PANEL_FILENAME = "integrated_call_samples_v3.20130502.ALL.panel"

#: Сигнальное значение для параметров "сколько образцов ожидаем": проверку
#: не выполнять вовсе (обратная совместимость для вызывающего кода, которому
#: размер подвыборки не известен и не важен).
EUR_COUNT_UNCHECKED = -1


def vcf_sample_count(vcf_path: Path) -> Optional[int]:
    """
    Число образцов в VCF по строке #CHROM, или None, если прочитать не
    удалось (файла нет, он битый, заголовка нет).

    errors="replace" — bcftools дописывает в заголовок командную строку с
    полным путём к файлу, и на Windows при не-ASCII символах в пути она
    может оказаться не в UTF-8 (та же причина, что в main.py::
    _donor_sample_order()). Имена образцов всегда чистый ASCII.
    """
    path = Path(vcf_path)
    if not path.exists() or path.stat().st_size == 0:
        return None
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return len(line.rstrip("\r\n").split("\t")[9:])
                if not line.startswith("#"):
                    break
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    return None


def available_eur_count(donors_dir: Path) -> Optional[int]:
    """
    Сколько EUR-образцов есть в файле панели 1000 Genomes, лежащем в папке
    доноров (503 в phase3). None — файла панели нет или он не читается,
    значит "все доступные" не с чем сравнивать.
    """
    panel = Path(donors_dir) / SAMPLES_PANEL_FILENAME
    if not panel.exists():
        return None
    try:
        text = panel.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    count = 0
    for line in text.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3 and parts[2] == "EUR":
            count += 1
    return count or None


class EurCountVerdict(NamedTuple):
    """
    ok=False — кэш собран на другом числе образцов, его нужно перекачать.
    cached/expected — числа для сообщения пользователю (expected может быть
    None, если "все доступные" не удалось выразить в числах).
    """
    ok: bool
    cached: Optional[int]
    expected: Optional[int]


def eur_count_verdict(
    donors_dir: Path,
    reference_vcf: Path,
    requested: Optional[int] = EUR_COUNT_UNCHECKED,
) -> EurCountVerdict:
    """
    Проверяет, что донорский кэш собран на запрошенном числе EUR-образцов.

    requested:
      EUR_COUNT_UNCHECKED (по умолчанию) — не проверять, всегда ok=True.
      None — "все доступные" (галочка в продвинутых настройках): кэш
          считается устаревшим, если в нём МЕНЬШЕ образцов, чем есть EUR
          в панели. Строгое "!=" тут неуместно: панель может пополниться
          не так, как ожидает программа, а больше доступного взять всё
          равно неоткуда.
      int — явное число: кэш годен только при точном совпадении.

    Если число образцов на диске определить не удалось (нет файлов,
    не читается заголовок) — ok=True: пусть решают обычные проверки
    наличия/пустоты файлов, а не эта.
    """
    if requested == EUR_COUNT_UNCHECKED:
        return EurCountVerdict(True, None, None)

    cached = vcf_sample_count(reference_vcf)
    if cached is None:
        return EurCountVerdict(True, None, None)

    if requested is None:
        expected = available_eur_count(donors_dir)
        if expected is None:
            return EurCountVerdict(True, cached, None)
        return EurCountVerdict(cached >= expected, cached, expected)

    expected = int(requested)
    return EurCountVerdict(cached == expected, cached, expected)
