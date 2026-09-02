"""
tests/test_x_chromosome.py

Регрессионные тесты на промт "Покрытие X-хромосомы" (жалоба Генотека:
~30% пропусков на X в собранном файле).

Причина была не в качестве импутации, а в том, что X через импутацию
вообще не проходила: весь пайплайн был жёстко зашит на хромосомы 1-22
(доноры, разбивка merged VCF, загрузка на MIS), поэтому в финальный файл
попадали только прямые измерения чипа, а все позиции трафарета, которых
на чипе нет, оставались "--". Живые цифры до фикса: аутосомы 97-99%
заполнения (импутация добавляла по 30-40 тыс. позиций на хромосому), X —
69.6% при 263 добавленных позициях.

Тесты проверяют четыре звена цепочки:
  1. доноры для X вообще запрашиваются, и под ПРАВИЛЬНЫМ именем файла
     (у chrX в 1000 Genomes phase3 свой суффикс, v1c, а не v5a/v5b);
  2. X не выбрасывается при разбивке merged VCF на файлы для загрузки;
  3. пол определяется по гетерозиготности nonPAR X, и мужской nonPAR
     пишется гаплоидно (Ploidy Check на стороне MIS), а PAR — нет;
  4. гаплоидный вызов из ответа сервера не отбрасывается сборщиком.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.base import ParsedVariant, ParseResult  # noqa: E402
from core.pure_python_core import (  # noqa: E402
    UPLOAD_CHROMS, build_vcf, infer_male_from_variants, is_par_position,
    normalise_x_ploidy,
)

# GRCh37: PAR1 60001-2699520, PAR2 154931044-155260560
PAR1_POS = 1_000_000
NONPAR_POS = 50_000_000
PAR2_POS = 155_000_000


# ---------------------------------------------------------------------------
# 1. Доноры
# ---------------------------------------------------------------------------
def test_donor_chroms_include_x():
    import download_donors as dd

    assert "X" in dd.DONOR_CHROMS, (
        "Без X в DONOR_CHROMS доноры для X не качаются, и X выпадает из "
        "impute-пайплайна целиком"
    )
    assert dd.DONOR_CHROMS[:22] == list(range(1, 23))


def test_x_uses_its_own_filename_suffix():
    """У chrX в релизе 20130502 суффикс v1c, у аутосом — v5a/v5b."""
    import download_donors as dd

    autosome = dd._vcf_suffix_candidates_for_build("grch37", 1)
    x = dd._vcf_suffix_candidates_for_build("grch37", "X")
    assert "v5a" in autosome and "v1c" in x
    assert not set(autosome) & set(x), (
        "Суффиксы аутосом и X не должны пересекаться — иначе перебор "
        "уходит на заведомо несуществующий URL"
    )
    url = dd._vcf_template_for_build("grch37").format(chrom="X", suffix=x[0])
    assert url == (
        "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
    )


def test_suffix_cache_key_separates_x_from_autosomes():
    """Кэш «рабочий суффикс на этом зеркале» общий на весь прогон — если
    не разделить ключи, найденный для аутосом v5a увёл бы X на 404."""
    import download_donors as dd

    mirror = "https://example.org/release/"
    assert dd._suffix_cache_key(mirror, 7) != dd._suffix_cache_key(mirror, "X")


def test_chrom_sort_key_handles_mixed_list():
    """sorted() на списке из int и "X" падает с TypeError — ключ нужен."""
    import download_donors as dd

    assert sorted(dd.DONOR_CHROMS, key=dd.chrom_sort_key)[-1] == "X"


# ---------------------------------------------------------------------------
# 2. Разбивка на файлы для загрузки
# ---------------------------------------------------------------------------
def test_upload_chroms_include_x_and_exclude_y_mt():
    assert UPLOAD_CHROMS[-1] == "X"
    assert "Y" not in UPLOAD_CHROMS and "MT" not in UPLOAD_CHROMS, (
        "Y и MT не поддерживаются ни HRC r1.1, ни 1000G Phase 3 как "
        "панели импутации — отправлять их нечего"
    )


# ---------------------------------------------------------------------------
# 3. Пол и плоидность
# ---------------------------------------------------------------------------
def test_par_boundaries_grch37():
    assert is_par_position(60_001)          # начало PAR1
    assert is_par_position(2_699_520)       # конец PAR1
    assert not is_par_position(2_699_521)   # первая позиция nonPAR
    assert not is_par_position(NONPAR_POS)
    assert is_par_position(155_000_000)     # PAR2


def _x_variants(n: int, het: int, pos_start: int = NONPAR_POS) -> list[ParsedVariant]:
    out = []
    for i in range(n):
        gt = "0/1" if i < het else "0/0"
        out.append(ParsedVariant(f"rs{i}", "X", pos_start + i, "A", "G", gt))
    return out


def test_infer_male_from_variants():
    male, pct, calls = infer_male_from_variants(_x_variants(1000, het=2))
    assert male and calls == 1000 and pct == pytest.approx(0.2)

    female, pct, _ = infer_male_from_variants(_x_variants(1000, het=300))
    assert not female and pct == pytest.approx(30.0)


def test_infer_male_needs_enough_calls():
    """Мало данных — не повод объявить образец мужским: диплоидный X
    сервер принимает всегда, гаплоидный по ошибке — испортит результат."""
    male, _, calls = infer_male_from_variants(_x_variants(10, het=0))
    assert not male and calls == 10


def test_infer_male_ignores_par():
    """В PAR у мужчины две копии, гетерозиготы там законны и не должны
    сбивать определение пола."""
    variants = _x_variants(1000, het=0) + _x_variants(200, het=200, pos_start=PAR1_POS)
    male, pct, calls = infer_male_from_variants(variants)
    assert male and calls == 1000 and pct == 0.0


def _build_and_read(tmp_path: Path, variants: list[ParsedVariant], **kwargs) -> list[str]:
    result = ParseResult()
    result.variants = variants
    out = tmp_path / "sample.vcf"
    build_vcf(result, out, sample_name="genotek", compress=False, **kwargs)
    return [
        l.rstrip("\n") for l in out.read_text(encoding="utf-8").splitlines()
        if not l.startswith("#")
    ]


def test_male_nonpar_x_is_haploid_par_is_diploid(tmp_path):
    variants = [
        ParsedVariant("rsPAR1", "X", PAR1_POS, "A", "G", "0/1"),
        ParsedVariant("rsNONPAR_ref", "X", NONPAR_POS, "A", ".", "0/0"),
        ParsedVariant("rsNONPAR_alt", "X", NONPAR_POS + 1, "A", "G", "1/1"),
        ParsedVariant("rsPAR2", "X", PAR2_POS, "C", "T", "0/1"),
        ParsedVariant("rsAuto", "1", 100, "A", "G", "0/1"),
    ]
    lines = _build_and_read(tmp_path, variants, haploid_x=True)
    gt = {l.split("\t")[2]: l.split("\t")[-1] for l in lines}
    assert gt["rsNONPAR_ref"] == "0", "мужской nonPAR X должен быть гаплоидным"
    assert gt["rsNONPAR_alt"] == "1"
    assert gt["rsPAR1"] == "0/1", "PAR остаётся диплоидным"
    assert gt["rsPAR2"] == "0/1"
    assert gt["rsAuto"] == "0/1", "аутосомы не затрагиваются"


def test_haploid_x_drops_heterozygous_nonpar_calls(tmp_path):
    """Гетерозигота в мужском nonPAR — шум чипа. Оставить её диплоидной
    нельзя (смешанная плоидность в nonPAR = провал Ploidy Check и отказ
    всего задания), а выбирать за чип один аллель — выдумывать данные."""
    variants = [
        ParsedVariant("rsHet", "X", NONPAR_POS, "A", "G", "0/1"),
        ParsedVariant("rsHom", "X", NONPAR_POS + 1, "A", ".", "0/0"),
    ]
    lines = _build_and_read(tmp_path, variants, haploid_x=True)
    rsids = [l.split("\t")[2] for l in lines]
    assert rsids == ["rsHom"]


def test_female_x_stays_diploid(tmp_path):
    variants = [ParsedVariant("rsX", "X", NONPAR_POS, "A", "G", "0/1")]
    lines = _build_and_read(tmp_path, variants, haploid_x=False)
    assert lines[0].split("\t")[-1] == "0/1"


# ---------------------------------------------------------------------------
# 4. Сборка результата
# ---------------------------------------------------------------------------
def _write_dose_vcf(path: Path, rows: list[str]) -> None:
    """Настоящий BGZF (не обычный gzip): сборщик индексирует файл через
    tabix, а тот обычный gzip не принимает."""
    bgzf = pytest.importorskip(
        "Bio.bgzf", reason="biopython не установлен — нечем записать BGZF",
    )
    with bgzf.BgzfWriter(str(path), "wb") as f:
        f.write(b"##fileformat=VCFv4.2\n")
        f.write(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tgenotek\n")
        for row in rows:
            f.write((row + "\n").encode("utf-8"))


def test_haploid_dose_genotype_is_accepted(tmp_path):
    """Ответ MIS для мужского X гаплоидный ("0"/"1"). Прежний сборщик
    отбрасывал такие строки жёстким `len(alleles) != 2 → continue`, то
    есть весь импутированный мужской X пропадал бы молча."""
    import shutil

    if not (shutil.which("bcftools") and shutil.which("tabix")):
        pytest.skip("bcftools/tabix не установлены — тест требует bcftools query")

    from template.assembler import load_imputed_genotypes

    imputed_dir = tmp_path / "rerun_results"
    imputed_dir.mkdir()
    # Строки строго по возрастанию позиции — tabix иначе откажется
    # индексировать файл (PAR1 физически идёт раньше nonPAR).
    _write_dose_vcf(imputed_dir / "chrX.dose.vcf.gz", [
        f"X\t{PAR1_POS}\trsDip\tC\tT\t.\tPASS\t.\tGT\t0/1",
        f"X\t{NONPAR_POS}\trsHap\tA\tG\t.\tPASS\t.\tGT\t1",
    ])
    genotypes = load_imputed_genotypes(imputed_dir, sample_name="genotek")
    assert genotypes[f"X_{NONPAR_POS}"] == "GG", (
        "гаплоидный вызов должен записываться удвоенным — так же, как "
        "прямые измерения чипа на X в формате 23andMe v3"
    )
    assert genotypes[f"X_{PAR1_POS}"] == "CT"


def test_x_dose_file_variants_are_found(tmp_path):
    """MIS отдаёт X одним chrX.dose.vcf.gz, но часть версий/зеркал — тремя
    кусками с собственными именами."""
    from template.assembler import _dose_file_pairs

    imputed_dir = tmp_path / "rerun_results"
    imputed_dir.mkdir()
    for name in ("chrX.no.auto_male.dose.vcf.gz", "chrX.par1.dose.vcf.gz"):
        (imputed_dir / name).write_bytes(b"")
    pairs = _dose_file_pairs(imputed_dir, "X")
    assert {p[0].name for p in pairs} == {
        "chrX.no.auto_male.dose.vcf.gz", "chrX.par1.dose.vcf.gz",
    }
    assert [p[1].name for p in pairs] == [
        "chrX.no.auto_male.info.gz", "chrX.par1.info.gz",
    ]


def test_autosome_dose_lookup_unchanged(tmp_path):
    from template.assembler import _dose_file_pairs

    imputed_dir = tmp_path / "rerun_results"
    imputed_dir.mkdir()
    (imputed_dir / "chr7.dose.vcf.gz").write_bytes(b"")
    pairs = _dose_file_pairs(imputed_dir, "7")
    assert [p[0].name for p in pairs] == ["chr7.dose.vcf.gz"]
    assert _dose_file_pairs(imputed_dir, "MT") == []


# ---------------------------------------------------------------------------
# 5. Плоидность после merge (Ploidy Check на стороне MIS)
# ---------------------------------------------------------------------------
# Живой прогон провалился на QC именно здесь:
#   Error: ChrX nonPAR region includes ambiguous samples (haploid and
#   diploid positions). Imputation cannot be started!
# Причина: сервер считает пропуск "./." ДИПЛОИДНОЙ записью, а bcftools
# merge подставлял его нашему гаплоидному мужскому образцу на позициях,
# которые есть у доноров и отсутствуют на чипе. Тринадцати таких записей
# из 28 900 хватило, чтобы задание было отвергнуто целиком.
def _x_line(pos: int, *gts: str) -> str:
    return "\t".join(["X", str(pos), f"rs{pos}", "A", "G", ".", "PASS", ".", "GT", *gts])


def test_missing_gt_follows_sample_own_ploidy():
    """"./." у гаплоидного образца → ".", и наоборот."""
    header = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\tfemale"]
    lines = header + [
        _x_line(NONPAR_POS + i, "0", "0/1") for i in range(10)
    ] + [_x_line(NONPAR_POS + 100, "./.", "./.")]

    out, fixed = normalise_x_ploidy(lines)
    assert fixed == 1
    last = out[-1].split("\t")
    assert last[9] == "."      # гаплоидный образец — гаплоидный пропуск
    assert last[10] == "./."   # диплоидный не тронут


def test_par_is_never_touched():
    """В PAR две копии у всех, проверка сервера туда не распространяется."""
    header = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale"]
    lines = header + [_x_line(NONPAR_POS + i, "0") for i in range(10)] + [
        _x_line(PAR1_POS, "./."), _x_line(PAR2_POS, "0/1"),
    ]
    out, _ = normalise_x_ploidy(lines)
    assert out[-2].split("\t")[9] == "./."
    assert out[-1].split("\t")[9] == "0/1"


def test_homozygous_diploid_call_is_compressed_for_haploid_sample():
    header = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale"]
    lines = header + [_x_line(NONPAR_POS + i, "1") for i in range(10)] + [
        _x_line(NONPAR_POS + 100, "1/1"),
    ]
    out, fixed = normalise_x_ploidy(lines)
    assert fixed == 1 and out[-1].split("\t")[9] == "1"


def test_heterozygous_call_for_haploid_sample_becomes_missing():
    """Гетерозигота гаплоидной быть не может, а выбирать за прибор один
    из двух аллелей — выдумывать данные."""
    header = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale"]
    lines = header + [_x_line(NONPAR_POS + i, "0") for i in range(10)] + [
        _x_line(NONPAR_POS + 100, "0/1"),
    ]
    out, fixed = normalise_x_ploidy(lines)
    assert fixed == 1 and out[-1].split("\t")[9] == "."


def test_no_sample_is_left_ambiguous():
    """Инвариант, который и проверяет сервер: у каждого образца в nonPAR
    ровно одна плоидность, считая пропуски."""
    header = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc"]
    lines = header + [
        _x_line(NONPAR_POS, "0", "0/1", "1"),
        _x_line(NONPAR_POS + 1, "./.", ".", "1/1"),
        _x_line(NONPAR_POS + 2, "1", "1/1", "0"),
    ]
    out, _ = normalise_x_ploidy(lines)
    for col in (9, 10, 11):
        kinds = set()
        for line in out[1:]:
            gt = line.split("\t")[col].rstrip("\n")
            kinds.add("dip" if ("/" in gt or "|" in gt) else "hap")
        assert len(kinds) == 1, f"колонка {col} осталась неоднородной: {kinds}"


def test_autosomes_are_not_affected(tmp_path):
    """normalise_x_ploidy вызывается только для X — аутосомы проходят
    через split_autosomes байт в байт."""
    from core.pure_python_core import split_autosomes

    pytest.importorskip("Bio.bgzf", reason="biopython не установлен")
    merged = tmp_path / "merged.vcf"
    merged.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t./.\n",
        encoding="utf-8", newline="\n",
    )
    split_autosomes(merged, tmp_path / "upload")
    with gzip.open(tmp_path / "upload" / "chr1.vcf.gz", "rt") as f:
        data = [l for l in f if not l.startswith("#")]
    assert data == ["1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t./.\n"]
