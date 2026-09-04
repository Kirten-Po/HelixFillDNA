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


# ---------------------------------------------------------------------------
# 6. Согласование образцов доноров перед bcftools concat
# ---------------------------------------------------------------------------
# Два живых падения на этом шаге:
#   Different sample names in kgp_sub_X.vcf.gz    — порядок колонок
#   Different number of samples in kgp_sub_X...   — размер подвыборки
# Второе к X отношения не имеет: chip_signature.txt описывает ЧИП, а не
# eur_sample_count, поэтому кэш, собранный на 20 образцах, проходит
# проверку сигнатуры, а доскачанная хромосома приходит уже с текущим
# значением по умолчанию (503). Добавление X просто вскрыло расхождение.
def _fake_donor(tmp_path: Path, name: str, samples: list[str]) -> Path:
    path = tmp_path / name
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                     "FILTER", "INFO", "FORMAT", *samples]) + "\n",
        encoding="utf-8", newline="\n",
    )
    return path


def test_donor_sample_order_is_read(tmp_path):
    import main as pipeline

    donor = _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2"])
    assert pipeline._donor_sample_order(donor) == ["HG1", "HG2"]


def test_align_is_noop_when_donors_already_match(tmp_path, monkeypatch):
    import main as pipeline

    calls = []
    monkeypatch.setattr(pipeline, "_run_bcftools", lambda args: calls.append(args))
    donors = [
        _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2"]),
        _fake_donor(tmp_path, "kgp_sub_X.vcf", ["HG1", "HG2"]),
    ]
    pipeline._align_donor_samples(donors)
    assert calls == [], "совпадающие доноры не должны переписываться"


def test_align_fixes_sample_order(tmp_path, monkeypatch):
    """Одинаковый набор, разный порядок колонок — правится только второй
    файл, и по порядку первого."""
    import main as pipeline

    calls = []

    def fake_run(args):
        order_file = Path(args[args.index("-S") + 1])
        calls.append((args[args.index("-S") + 2],
                      order_file.read_text(encoding="utf-8").split()))
        # имитируем результат bcftools, чтобы .replace() нашёл файл
        Path(args[args.index("-o") + 1]).write_text("", encoding="utf-8")

    monkeypatch.setattr(pipeline, "_run_bcftools", fake_run)
    monkeypatch.setattr(pipeline, "_index_vcf", lambda p: None)
    donors = [
        _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2", "HG3"]),
        _fake_donor(tmp_path, "kgp_sub_X.vcf", ["HG3", "HG1", "HG2"]),
    ]
    pipeline._align_donor_samples(donors)
    assert len(calls) == 1
    target_file, order = calls[0]
    assert "kgp_sub_X" in str(target_file)
    assert order == ["HG1", "HG2", "HG3"]


def test_align_takes_intersection_when_sample_counts_differ(tmp_path, monkeypatch):
    """Кэш на 20 образцах + доскачанная хромосома на 503 → общие 20."""
    import main as pipeline

    written = {}

    def fake_run(args):
        order_file = Path(args[args.index("-S") + 1])
        written[str(args[args.index("-S") + 2])] = order_file.read_text(
            encoding="utf-8").split()

    monkeypatch.setattr(pipeline, "_run_bcftools", fake_run)
    monkeypatch.setattr(pipeline, "_index_vcf", lambda p: None)
    monkeypatch.setattr(Path, "replace", lambda self, target: None)

    donors = [
        _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2"]),
        _fake_donor(tmp_path, "kgp_sub_X.vcf", ["HG9", "HG2", "HG1", "HG7"]),
    ]
    monkeypatch.setattr(pipeline, "MIN_DONOR_SAMPLES", 2)
    pipeline._align_donor_samples(donors)
    assert len(written) == 1
    assert list(written.values())[0] == ["HG1", "HG2"], (
        "порядок берётся у первого донора, набор — пересечение"
    )


def test_align_refuses_when_intersection_too_small(tmp_path, monkeypatch):
    """Пересечение почти пустое — это испорченный кэш, а не то, что можно
    молча «починить» урезанием."""
    import main as pipeline

    monkeypatch.setattr(pipeline, "_run_bcftools", lambda args: None)
    donors = [
        _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2", "HG3", "HG4", "HG5"]),
        _fake_donor(tmp_path, "kgp_sub_X.vcf", ["HG8", "HG9", "HG5", "HG6", "HG7"]),
    ]
    with pytest.raises(RuntimeError, match="общих образцов"):
        pipeline._align_donor_samples(donors)


def test_align_does_not_second_guess_a_small_but_consistent_cache(tmp_path, monkeypatch):
    """Порог применяется только когда наборы расходятся. Если все доноры
    согласованы, их число — осознанный выбор пользователя
    (eur_sample_count), и ронять запуск из-за него нельзя."""
    import main as pipeline

    monkeypatch.setattr(pipeline, "_run_bcftools",
                        lambda args: pytest.fail("не должно вызываться"))
    donors = [
        _fake_donor(tmp_path, "kgp_sub_1.vcf", ["HG1", "HG2"]),
        _fake_donor(tmp_path, "kgp_sub_X.vcf", ["HG1", "HG2"]),
    ]
    pipeline._align_donor_samples(donors)   # не бросает


# ---------------------------------------------------------------------------
# 7. Версия формата VCF в файлах для загрузки
# ---------------------------------------------------------------------------
# Живой прогон на панели TOPMed провалился на QC:
#   Task 'Calculating QC Statistics Writing VCF version VCF4_3 is not
#   implemented
# QC-шаг сервера умеет писать VCF только до 4.2, а GRCh38-релиз
# 1000 Genomes (доноры для TOPMed) объявляет себя VCFv4.3 — и bcftools
# merge поднимает версию всего задания до неё. К X отношения не имеет:
# на HRC доноры 4.1, поэтому там это никогда не всплывало.
def test_upload_files_declare_vcf42(tmp_path):
    from core.pure_python_core import UPLOAD_VCF_VERSION, split_autosomes

    pytest.importorskip("Bio.bgzf", reason="biopython не установлен")
    assert UPLOAD_VCF_VERSION == "VCFv4.2"

    merged = tmp_path / "merged.vcf"
    merged.write_text(
        "##fileformat=VCFv4.3\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8", newline="\n",
    )
    split_autosomes(merged, tmp_path / "upload")

    with gzip.open(tmp_path / "upload" / "chr1.vcf.gz", "rt") as f:
        lines = f.read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    # остальной заголовок и данные не тронуты
    assert lines[1] == '##FILTER=<ID=PASS,Description="All filters passed">'
    assert lines[-1] == "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"


# ---------------------------------------------------------------------------
# 8. GRCh38 + chrX: другой набор данных
# ---------------------------------------------------------------------------
# Найдено живым прогоном на TOPMed: донор для X дал 610 позиций из ~30 000.
# Файл ALL.chrX...v2a_27022019.GRCh38.phased.vcf.gz из релиза
# 20190312_biallelic_SNV_and_INDEL содержит ТОЛЬКО псевдоаутосомные
# регионы — 97 875 вариантов в PAR1, 9 088 в PAR2 и ровно ноль в nonPAR.
# Это особенность релиза (тот же файл лежит и на зеркале UCSC), поэтому
# nonPAR для GRCh38 берётся из набора 30x high-coverage, где у X ещё и
# собственный суффикс имени.
def test_grch38_x_uses_high_coverage_dataset():
    import download_donors as dd

    mirror = dd._mirrors_for_build("grch38", "X")[0]
    assert "20220422_3202_phased_SNV_INDEL_SV" in mirror
    name = dd._vcf_template_for_build("grch38", "X").format(
        chrom="X", suffix=dd._vcf_suffix_candidates_for_build("grch38", "X")[0],
    )
    assert name == (
        "1kGP_high_coverage_Illumina.chrX"
        ".filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    )


def test_grch38_autosomes_keep_their_own_dataset():
    """Смена источника касается ТОЛЬКО X: аутосомы GRCh38 остаются на
    релизе 20190312."""
    import download_donors as dd

    assert "20190312" in dd._mirrors_for_build("grch38", 1)[0]
    assert "shapeit2_integrated_snvindels" in dd._vcf_template_for_build("grch38", 1)
    # и без указания хромосомы поведение прежнее
    assert dd._mirrors_for_build("grch38") == dd._mirrors_for_build("grch38", 7)


def test_grch37_x_is_not_affected():
    """У GRCh37 chrX полноценный, менять источник незачем."""
    import download_donors as dd

    assert "20130502" in dd._mirrors_for_build("grch37", "X")[0]
    assert dd._vcf_suffix_candidates_for_build("grch37", "X") == ["v1c", "v1b"]


def test_grch38_x_subset_cache_has_its_own_name():
    """Кэш подвыборки, построенный из старого PAR-файла, не должен
    совпасть по имени с новым — иначе он молча переиспользуется и X
    снова окажется почти пустой."""
    import download_donors as dd

    cache = Path("/cache")
    assert dd._eur_subset_cache_path(cache, "X", 503, "grch38").name == (
        "EUR503.chrX.highcov.vcf.gz"
    )
    assert dd._eur_subset_cache_path(cache, "X", 503, "grch37").name == (
        "EUR503.chrX.vcf.gz"
    )
    assert dd._eur_subset_cache_path(cache, 1, 503, "grch38").name == (
        "EUR503.chr1.vcf.gz"
    )


# ---------------------------------------------------------------------------
# 9. Порядок тегов в заголовке VCF
# ---------------------------------------------------------------------------
# Задание на TOPMed отвергалось ещё на валидации входа:
#   Unable to parse header ... Tag Type in wrong order (was #2, expected #3)
#   in line <ID=END2,Type=Integer,Number=1,Description="...">
# Спецификация VCF порядок тегов не фиксирует, htsjdk на стороне сервера —
# требует. Кривая строка приезжает из набора 30x high-coverage (поле
# структурных вариантов END2), а bcftools concat разносит объединённый
# заголовок по всем 23 файлам — одна строка в донорах chrX роняла весь набор.
def test_header_tag_order_is_normalised():
    from core.pure_python_core import normalise_structured_header_line as fix

    assert fix(
        '##INFO=<ID=END2,Type=Integer,Number=1,Description="Position of breakpoint on CHR2">'
    ) == (
        '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">'
    )


def test_commas_inside_description_are_not_split():
    from core.pure_python_core import normalise_structured_header_line as fix

    line = '##INFO=<ID=AF,Type=Float,Number=A,Description="frequency, range (0,1)">'
    assert fix(line) == (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="frequency, range (0,1)">'
    )


def test_correct_lines_are_left_untouched():
    from core.pure_python_core import normalise_structured_header_line as fix

    for line in (
        '##INFO=<ID=AF,Number=A,Type=Float,Description="freq">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##contig=<ID=chr1,length=248956422>\n',
        '##fileDate=26022019\n',
        '##fileformat=VCFv4.2\n',
    ):
        assert fix(line) == line


def test_extra_tags_are_preserved_after_the_four(tmp_path):
    from core.pure_python_core import normalise_structured_header_line as fix

    line = '##INFO=<ID=X,Type=Integer,Number=1,Source="dbsnp",Version="2">'
    assert fix(line) == (
        '##INFO=<ID=X,Number=1,Type=Integer,Source="dbsnp",Version="2">'
    )


def test_split_autosomes_fixes_the_header(tmp_path):
    """Сквозная проверка: кривая строка не должна доехать до upload/."""
    from core.pure_python_core import split_autosomes

    pytest.importorskip("Bio.bgzf", reason="biopython не установлен")
    merged = tmp_path / "merged.vcf"
    merged.write_text(
        "##fileformat=VCFv4.3\n"
        '##INFO=<ID=END2,Type=Integer,Number=1,Description="breakpoint">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8", newline="\n",
    )
    split_autosomes(merged, tmp_path / "upload")
    with gzip.open(tmp_path / "upload" / "chr1.vcf.gz", "rt") as f:
        head = f.read().splitlines()
    assert head[0] == "##fileformat=VCFv4.2"
    assert head[1] == (
        '##INFO=<ID=END2,Number=1,Type=Integer,Description="breakpoint">'
    )
