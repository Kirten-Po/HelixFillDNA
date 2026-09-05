"""
tests/test_indel_filtering.py

Регрессионный тест на промт "странные генотипы в результате TOPMed".

В итоговом файле, собранном с панели TOPMed, встречались генотипы вида
"TGTGATGTGA", "CACA", "GCG" — вплоть до строк на полсотни символов.
Причина двойная:

  1. Панель TOPMed r3, в отличие от HRC r1.1, содержит не только SNP, но
     и инделы. У инделя REF или ALT — строка из нескольких букв, а
     генотип собирался как a1 + a2, то есть простой склейкой.

  2. Инделя и SNP нередко делят ОДНУ координату, а словарь генотипов
     ключуется по (хромосома, позиция) — значит запись, прочитанная
     позже, затирала прочитанную раньше. Живой пример, chr1:11195977:
         rs17036508   REF=T      ALT=C  -> "TT"          (нужный SNP)
         rs533913726  REF=TGTGA  ALT=T  -> "TGTGATGTGA"  (индель)
     Индель шёл вторым и затирал правильный генотип.

То есть это была не косметика, а подмена значения: 178 позиций из
959 708 в файле были просто неверны. На HRC симптом не проявлялся —
там панель состоит только из SNV.
"""
from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

pytestmark = pytest.mark.skipif(
    not (shutil.which("bcftools") and shutil.which("tabix")),
    reason="нужны bcftools и tabix",
)

HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "##contig=<ID=1>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tgenotek\n"
)


def _write_dose(path: Path, rows: list[str]) -> None:
    bgzf = pytest.importorskip("Bio.bgzf", reason="нужен biopython")
    with bgzf.BgzfWriter(str(path), "wb") as f:
        f.write(HEADER.encode())
        for row in rows:
            f.write((row + "\n").encode())


def test_indel_does_not_overwrite_the_snp_at_the_same_position(tmp_path):
    """Точная реконструкция chr1:11195977 из реального результата."""
    from template.assembler import load_imputed_genotypes

    d = tmp_path / "rerun_results"
    d.mkdir()
    _write_dose(d / "chr1.dose.vcf.gz", [
        "1\t11195977\trs17036508\tT\tC\t.\tPASS\t.\tGT\t0|0",
        "1\t11195977\trs533913726\tTGTGA\tT\t.\tPASS\t.\tGT\t0|0",
    ])
    genotypes = load_imputed_genotypes(d, sample_name="genotek")
    assert genotypes["1_11195977"] == "TT", (
        "должен остаться генотип SNP, а не склейка аллелей инделя"
    )


def test_insertion_is_skipped_too(tmp_path):
    """Многобуквенным может быть и ALT, а не только REF."""
    from template.assembler import load_imputed_genotypes

    d = tmp_path / "rerun_results"
    d.mkdir()
    _write_dose(d / "chr1.dose.vcf.gz", [
        "1\t500\trs1\tC\tCA\t.\tPASS\t.\tGT\t1|1",
        "1\t600\trs2\tA\tG\t.\tPASS\t.\tGT\t1|1",
    ])
    genotypes = load_imputed_genotypes(d, sample_name="genotek")
    assert "1_500" not in genotypes
    assert genotypes["1_600"] == "GG"


def test_every_genotype_is_two_letters(tmp_path):
    """Инвариант формата 23andMe: генотип — ровно две буквы (для
    гаплоидных вызовов сборщик удваивает аллель)."""
    from template.assembler import load_imputed_genotypes

    d = tmp_path / "rerun_results"
    d.mkdir()
    _write_dose(d / "chr1.dose.vcf.gz", [
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1",
        "1\t200\trs2\tACGT\tA\t.\tPASS\t.\tGT\t0|0",
        "1\t300\trs3\tT\tTTTTT\t.\tPASS\t.\tGT\t1|1",
        "1\t400\trs4\tC\tT\t.\tPASS\t.\tGT\t1",
    ])
    genotypes = load_imputed_genotypes(d, sample_name="genotek")
    assert all(len(g) == 2 for g in genotypes.values()), genotypes
    assert set(genotypes) == {"1_100", "1_400"}


def test_snps_are_not_affected(tmp_path):
    """Фильтр не должен трогать обычные SNP — путь HRC не меняется."""
    from template.assembler import load_imputed_genotypes

    d = tmp_path / "rerun_results"
    d.mkdir()
    _write_dose(d / "chr1.dose.vcf.gz", [
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|0",
        "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1",
        "1\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t1|1",
    ])
    genotypes = load_imputed_genotypes(d, sample_name="genotek")
    assert genotypes == {"1_100": "AA", "1_200": "CT", "1_300": "AA"}


# ---------------------------------------------------------------------------
# Последний рубеж: проверка содержимого колонки генотипа
# ---------------------------------------------------------------------------
# 178 испорченных значений прошли ВСЕ прежние проверки незамеченными:
# число строк, число полей, переносы и идентичность первых трёх колонок
# к содержимому четвёртой отношения не имеют. Эта проверка ловит и
# инделы, и любую другую порчу значения — до того, как файл уедет
# в Генотек.
def _mini_template(tmp_path: Path, rows: list[str]) -> Path:
    p = tmp_path / "template.txt"
    p.write_text(
        "# rsid\tchromosome\tposition\tgenotype\n" + "".join(r + "\n" for r in rows),
        encoding="utf-8", newline="",
    )
    return p


def _mini_output(tmp_path: Path, rows: list[str]) -> Path:
    p = tmp_path / "out.txt"
    p.write_text(
        "# rsid\tchromosome\tposition\tgenotype\n" + "".join(r + "\n" for r in rows),
        encoding="utf-8", newline="",
    )
    return p


@pytest.mark.parametrize("genotype", ["AA", "CT", "A", "--", "DD", "DI", "II", "I"])
def test_valid_genotypes_pass(tmp_path, genotype):
    from template.assembler import validate_output

    tpl = _mini_template(tmp_path, ["rs1\t1\t100\tAA"])
    out = _mini_output(tmp_path, [f"rs1\t1\t100\t{genotype}"])
    assert validate_output(out, tpl, "v3").invalid_genotypes == []


@pytest.mark.parametrize("genotype", ["TGTGATGTGA", "CACA", "GCG", "ACGTACGT", "XX", "A1"])
def test_invalid_genotypes_are_caught(tmp_path, genotype):
    from template.assembler import validate_output

    tpl = _mini_template(tmp_path, ["rs1\t1\t100\tAA"])
    out = _mini_output(tmp_path, [f"rs1\t1\t100\t{genotype}"])
    result = validate_output(out, tpl, "v3")
    assert len(result.invalid_genotypes) == 1
    assert result.invalid_genotypes[0][2] == genotype
    assert not result.is_valid, "файл с таким значением не должен считаться валидным"
    assert any("генотип" in e.lower() for e in result.errors), result.errors


# ---------------------------------------------------------------------------
# Повторная сборка из уже скачанных результатов
# ---------------------------------------------------------------------------
# Path.with_suffix() заменяет только ПОСЛЕДНЕЕ расширение: для
# "chr1.dose.vcf.gz" вызов .with_suffix(".vcf.gz.tbi") давал
# "chr1.dose.vcf.vcf.gz.tbi" — несуществующий путь. Проверка "индекс уже
# есть" не срабатывала никогда, tabix запускался всегда, и на ПОВТОРНОЙ
# сборке падал: "[tabix] the index file exists. Please use '-f'".
# Первая сборка проходила, поэтому баг годами не проявлялся — а пересборка
# из уже скачанных rerun_results это штатный сценарий.
def test_reassembly_with_existing_index(tmp_path):
    from template.assembler import load_imputed_genotypes

    d = tmp_path / "rerun_results"
    d.mkdir()
    _write_dose(d / "chr1.dose.vcf.gz", [
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1",
    ])
    first = load_imputed_genotypes(d, sample_name="genotek")
    assert (d / "chr1.dose.vcf.gz.tbi").exists(), "индекс должен лечь рядом с VCF"

    # Вторая сборка тем же кодом — раньше здесь падало.
    second = load_imputed_genotypes(d, sample_name="genotek")
    assert second == first == {"1_100": "AG"}
