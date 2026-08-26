"""
tests/test_reference_genome_contigs.py

Юнит-тест на adapters.ftdna_v3.ReferenceGenome.base_at() — проверяет
резолвинг контигов с/без префикса "chr" (Промт_TopMed_HRC_v2.md, п.2/9):

  - GRCh37-подобный референс (контиги "1", "2", ..., "MT" без префикса)
    — как сейчас работает HRC.
  - GRCh38-подобный референс (контиги "chr1", "chr2", ..., "chrM") — как
    ожидается для TopMed. Проверяется, что base_at() резолвит ОБА
    варианта запроса позиции: канонический ("1") и с явным "chr"
    ("chr1"), а также что канонический запрос "MT" находит контиг
    "chrM" (не только "chrMT").

Требует pyfaidx (уже зависимость проекта, adapters/ftdna_v3.py его
импортирует безусловно). Создаёт временные .fasta-файлы через pytest
tmp_path — никаких сетевых обращений и реальных многогигабайтных
референсов не требуется.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.ftdna_v3 import ReferenceGenome, FTDNAFormatError  # noqa: E402


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    """Пишет минимальный многозаписной .fasta-файл (без переносов внутри
    последовательности — для теста этого достаточно, длины контигов
    маленькие)."""
    with path.open("w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")
    return path


@pytest.fixture
def grch37_like_fasta(tmp_path: Path) -> Path:
    """Контиги без префикса 'chr' — как в human_g1k_v37.fasta (HRC)."""
    return _write_fasta(tmp_path / "grch37_like.fasta", {
        "1": "ACGTACGTAC",
        "MT": "GGGGCCCCAA",
    })


@pytest.fixture
def grch38_like_fasta(tmp_path: Path) -> Path:
    """Контиги с префиксом 'chr', митохондрия названа 'chrM' (не 'chrMT')
    — типичная особенность части GRCh38-сборок (см. докстринг
    ReferenceGenome._build_contig_alias_map)."""
    return _write_fasta(tmp_path / "grch38_like.fasta", {
        "chr1": "ACGTACGTAC",
        "chrM": "GGGGCCCCAA",
    })


def test_grch37_style_contig_resolves_without_prefix(grch37_like_fasta):
    ref = ReferenceGenome(grch37_like_fasta)
    # Канонический запрос ("1", без префикса) — как приходит от адаптеров
    # (_normalize_chrom()) — должен резолвиться напрямую к контигу "1".
    assert ref.base_at("1", 1) == "A"
    assert ref.base_at("MT", 1) == "G"


def test_grch38_style_contig_resolves_with_and_without_prefix(grch38_like_fasta):
    ref = ReferenceGenome(grch38_like_fasta)
    # Канонический запрос "1" должен найти контиг "chr1" через alias-словарь.
    assert ref.base_at("1", 1) == "A"
    # Явный запрос с префиксом тоже должен сработать (напрямую по ключу
    # fasta, либо через тот же alias-словарь, если контиг переименован).
    assert ref.base_at("chr1", 1) == "A"
    # "MT" (канонический код митохондрии из адаптеров) должен резолвиться
    # к контигу "chrM" — не только к гипотетическому "chrMT".
    assert ref.base_at("MT", 1) == "G"


def test_unknown_chromosome_still_raises_format_error(grch38_like_fasta):
    ref = ReferenceGenome(grch38_like_fasta)
    with pytest.raises(FTDNAFormatError):
        ref.base_at("26", 1)


def test_alias_map_built_once_at_init(grch38_like_fasta):
    """Убеждаемся, что alias-словарь строится ровно один раз при
    инициализации (сравниваем объект словаря до и после нескольких
    вызовов base_at() — если бы он пересчитывался на каждый вызов, это
    были бы разные dict-объекты)."""
    ref = ReferenceGenome(grch38_like_fasta)
    aliases_before = ref._contig_aliases
    ref.base_at("1", 1)
    ref.base_at("MT", 1)
    assert ref._contig_aliases is aliases_before
