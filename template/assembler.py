"""
template/assembler.py
Сборка финального файла по трафарету + 7 проверок из Части 10 гайда.
С автоматической фильтрацией по Rsq и поддержкой хромосом X, Y, MT.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import gzip
import logging
import shutil
import subprocess
from typing import Iterator
from core.pure_python_core import read_rsq_map
from .skeleton import SkeletonRow, extract_skeleton

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Структуры данных
# ---------------------------------------------------------------------------
@dataclass
class ValidationResult:
    total_lines: int
    expected_lines: int
    lines_match: bool
    fields_per_line: set[int]
    fields_valid: bool
    crlf_count: int
    crlf_expected: bool
    crlf_match: bool
    chromosomes: list[str]
    call_rate: float
    structure_identical: bool
    errors: list[str]
    # Промт "Дублирующая позиция после норм. -m-both": позиции, которые
    # встречаются в output файле более одного раза. Доказуемо НЕ могут
    # быть привнесены самой сборкой (assemble_final() пишет ровно по
    # одной строке на каждую запись skeleton, а skeleton строится 1:1 из
    # template_path) — то есть это всегда унаследованная особенность
    # самого трафарета (например, два разных rsid на одной физической
    # координате чипа — известная и легитимная особенность части реальных
    # 23andMe-экспортов), а не ошибка конвертации. Поэтому НЕ входит в
    # errors/is_valid — только информационное поле для прозрачности.
    template_duplicate_positions: list[tuple[str, int]]

    @property
    def is_valid(self) -> bool:
        return (
            self.lines_match
            and self.fields_valid
            and self.crlf_match
            and self.structure_identical
            and not self.errors
        )


class AssemblyError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Загрузка импутированных генотипов (с фильтрацией по Rsq и поддержкой X/Y/MT)
# ---------------------------------------------------------------------------
def load_imputed_genotypes(
    imputed_dir: Path,
    sample_name: str = "genotek",
    panel_pos: list[tuple[str, int]] | None = None,
    rsq_threshold: float = 0.30,
    bcftools_path: str | None = None,
    tabix_path: str | None = None,
) -> dict[str, str]:
    """
    Загружает импутированные генотипы из chr*.dose.vcf.gz.
    Автоматически читает chr*.info.gz и отбрасывает варианты с Rsq < rsq_threshold.
    Обрабатывает хромосомы 1-22, X, Y, MT.
    """
    imputed_dir = Path(imputed_dir)
    if not imputed_dir.is_dir():
        raise FileNotFoundError(f"Папка с импутированными данными не найдена: {imputed_dir}")

    bcftools = bcftools_path or shutil.which("bcftools") or "bcftools"
    tabix = tabix_path or shutil.which("tabix") or "tabix"

    genotypes: dict[str, str] = {}
    panel_set = set(panel_pos) if panel_pos else None

    # Расширенный список хромосом: 1-22 + X + Y + MT
    chrom_names = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

    for chrom in chrom_names:
        for vcf_path, info_path in _dose_file_pairs(imputed_dir, chrom):
            _load_one_dose_file(
                vcf_path, info_path, chrom, genotypes, panel_set,
                rsq_threshold, bcftools, tabix, sample_name, imputed_dir,
            )

    logger.info("Загружено %d импутированных генотипов (Rsq >= %.2f)", len(genotypes), rsq_threshold)
    return genotypes


# Michigan Imputation Server отдаёт X одним файлом chrX.dose.vcf.gz
# (сервер сам склеивает PAR1/nonPAR/PAR2 обратно), но исторические версии
# и часть зеркал отдают её тремя-четырьмя кусками с собственными именами
# (chrX.no.auto_male / chrX.no.auto_female / chrX.par1 / chrX.par2).
# Принимаем оба варианта: иначе результат импутации X, ради которого всё
# и затевалось, молча не попал бы в финальный файл.
_X_DOSE_GLOBS = ("chrX.*.dose.vcf.gz",)


def _dose_file_pairs(imputed_dir: Path, chrom: str) -> list[tuple[Path, Path]]:
    """Пары (dose.vcf.gz, info.gz) для одной хромосомы. Пустой список —
    хромосомы нет в результатах (нормально для Y/MT и для X, если задание
    отправлялось без неё)."""
    pairs: list[tuple[Path, Path]] = []
    main_vcf = imputed_dir / f"chr{chrom}.dose.vcf.gz"
    if main_vcf.exists():
        pairs.append((main_vcf, imputed_dir / f"chr{chrom}.info.gz"))
    if chrom == "X":
        for pattern in _X_DOSE_GLOBS:
            for extra in sorted(imputed_dir.glob(pattern)):
                if extra == main_vcf:
                    continue
                stem = extra.name[: -len(".dose.vcf.gz")]
                pairs.append((extra, imputed_dir / f"{stem}.info.gz"))
    return pairs


def _load_one_dose_file(
    vcf_path: Path,
    info_path: Path,
    chrom: str,
    genotypes: dict[str, str],
    panel_set,
    rsq_threshold: float,
    bcftools: str,
    tabix: str,
    sample_name: str,
    imputed_dir: Path,
) -> None:
    """Тело прежнего цикла по хромосомам, вынесенное в отдельную функцию:
    для X их теперь может быть несколько файлов на одну хромосому (см.
    _dose_file_pairs())."""

    # 1. Загружаем Rsq для этой хромосомы
    rsq_map = read_rsq_map(info_path)

    # 2. Индексируем VCF, если нужно
    tbi_path = vcf_path.with_suffix(".vcf.gz.tbi")
    if not tbi_path.exists():
        tabix_result = subprocess.run(
            [tabix, "-p", "vcf", str(vcf_path)], capture_output=True, text=True,
        )
        if tabix_result.returncode != 0:
            # Раньше здесь стоял check=True без вывода stderr — до
            # пользователя долетал только бесполезный
            # "Command '[...]' returned non-zero exit status 1"
            # без единого слова о РЕАЛЬНОЙ причине (файл не BGZF,
            # повреждён/усечён, не отсортирован и т.п.), хотя tabix
            # эту причину печатает в stderr. Теперь она попадает в
            # текст исключения и видна в логе GUI/CLI.
            raise AssemblyError(
                f"tabix не смог проиндексировать {vcf_path.name} "
                f"(код {tabix_result.returncode}):\n"
                f"{tabix_result.stderr.strip() or '(tabix не вывел никакого сообщения об ошибке)'}\n"
                f"Обычно это значит, что файл повреждён/усечён (неудачная "
                f"докачка/распаковка) или не является настоящим BGZF-"
                f"сжатым VCF. Попробуйте скачать результаты MIS заново "
                f"для этой хромосомы."
            )

    # 3. Извлекаем генотипы
    cmd = [
        bcftools, "query",
        "-s", sample_name,
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n",
        str(vcf_path),
    ]
    if panel_set:
        panel_file = imputed_dir / f"_panel_{chrom}.txt"
        # Промт "встроить лифтовер HRC/TopMed в gui/app.py", точечный
        # фикс (НЕ связанный с самим лифтовером координат — тот уже
        # решает свою часть задачи, перенос ПОЗИЦИЙ; здесь отдельная,
        # независимая проблема — ИМЕНОВАНИЕ контига).
        #
        # panel_pos приходит сюда с КАНОНИЧЕСКИМИ именами хромосом (без
        # префикса "chr") — так их всегда отдают и
        # template/skeleton.py::extract_skeleton() (GRCh37/HRC), и
        # main.py::liftover_positions_forward()/
        # core/liftover.py::ChainLiftover.lift() (после форвард-
        # лифтовера под TopMed — ChainLiftover тоже всегда возвращает
        # канонический вид). Но сам VCF-результат Michigan Imputation
        # Server для сборок с REFERENCE_PANELS[panel]["chrom_prefix"]
        # == "chr" (сейчас — TopMed/GRCh38, см. main.py) почти наверняка
        # использует CHROM="chr1".."chr22" — то есть regions-файл с
        # именами БЕЗ префикса не совпал бы ни с одной записью VCF, и
        # bcftools view -R молча вернул бы 0 строк для каждой
        # хромосомы: все импутированные генотипы для TopMed исчезли бы
        # целиком, без единой ошибки.
        #
        # Фикс агностичен к лифтоверу и не требует знать о нём здесь:
        # пишем в regions-файл ОБЕ формы имени хромосомы ("1" и "chr1").
        # bcftools -R не требует, чтобы каждое имя контига в файле
        # реально существовало в читаемом VCF — несовпавшие строки
        # просто не находят совпадений и безвредны. Для HRC/GRCh37
        # (где CHROM в самом VCF без "chr") совпадёт только форма без
        # префикса — поведение не меняется, идемпотентно.
        #
        # ⚠ Не покрывает возможную путаницу "chrMT" vs "chrM" для
        # митохондриальной хромосомы на некоторых GRCh38-релизах
        # (иногда митохондрию называют "chrM", а не "chr" + "MT" =
        # "chrMT") — это отдельный, неподтверждённый вопрос вне рамок
        # этого точечного фикса.
        # ⚠ Фикс (та же ловушка, что была найдена в
        # main.py::_post_merge_intersect(), "Failed to read the
        # regions"): без newline="" питоновский текстовый режим на
        # Windows транслирует "\n" в "\r\n" при записи — regions-файл
        # для bcftools -R получал CRLF-переносы на каждой строке. Это
        # не всегда роняет bcftools с ошибкой (в отличие от
        # common_pos.txt), но может тихо портить часть совпадений при
        # чтении позиций — прямой риск для ЭТОГО файла особенно
        # велик, так как он используется в bcftools query -R на Этапе
        # 7 (загрузка импутированных генотипов) — если часть позиций
        # молча не совпадёт, часть результата импутации попадёт в
        # финальный файл как "--" вместо реального генотипа, без
        # единой видимой ошибки.
        with panel_file.open("w", newline="\n") as f:
            for c, p in panel_set:
                if c == chrom:
                    f.write(f"{c}\t{p}\n")
                    f.write(f"chr{c}\t{p}\n")
        cmd.extend(["-R", str(panel_file)])

    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue

        c, p, ref, alt, gt = parts[0], parts[1], parts[2], parts[3], parts[4]
        c_norm = c.replace("chr", "")
        p_int = int(p)

        # === ГЛАВНОЕ: ПРОВЕРКА Rsq ===
        rsq = rsq_map.get((c_norm, p_int), 1.0)
        if rsq < rsq_threshold:
            continue  # Отбрасываем низкокачественные варианты

        if "," in alt:
            continue
        if gt in ("./.", "."):
            continue

        gt_norm = gt.replace("|", "/")
        alleles = gt_norm.split("/")
        # ГАПЛОИДНЫЙ вызов ("0" / "1") — нормальная и ожидаемая форма
        # для мужского nonPAR X: именно так мы отправляем эти позиции
        # на сервер (см. core/pure_python_core.py::build_vcf(haploid_x=))
        # и так же их возвращает Michigan Imputation Server. Раньше
        # здесь стояло жёсткое `len(alleles) != 2 → continue`, то есть
        # весь мужской X был бы отброшен целиком и молча.
        # В формат 23andMe v3 такой вызов записывается удвоенным
        # (как и прямые измерения чипа на X, см. load_measured_genotypes)
        # — так же, как это делал прежний файл, принятый Генотеком.
        if len(alleles) == 1:
            alleles = [alleles[0], alleles[0]]
        elif len(alleles) != 2:
            continue

        try:
            a1 = ref if alleles[0] == "0" else alt
            a2 = ref if alleles[1] == "0" else alt
            genotype = a1 + a2
            key = f"{c_norm}_{p_int}"
            genotypes[key] = genotype
        except (IndexError, ValueError):
            continue

    if panel_set:
        panel_file.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Загрузка реальных измерений
# ---------------------------------------------------------------------------
def load_measured_genotypes(variants: list) -> dict[str, str]:
    genotypes: dict[str, str] = {}
    for v in variants:
        if v.gt == "0/0":
            genotype = v.ref + v.ref
        elif v.gt in ("0/1", "1/0"):
            genotype = v.ref + v.alt
        elif v.gt == "1/1":
            genotype = v.alt + v.alt
        else:
            continue
        key = f"{v.chrom}_{v.pos}"
        genotypes[key] = genotype
    logger.info("Загружено %d реальных измерений", len(genotypes))
    return genotypes


# ---------------------------------------------------------------------------
# Слияние словарей
# ---------------------------------------------------------------------------
def merge_dictionaries(imputed: dict[str, str], measured: dict[str, str]) -> dict[str, str]:
    merged = imputed.copy()
    merged.update(measured)  # measured перезаписывает imputed
    logger.info("Слияние словарей: импутировано=%d, реальных=%d, итого=%d",
                len(imputed), len(measured), len(merged))
    return merged


# ---------------------------------------------------------------------------
# Заголовок трафарета
# ---------------------------------------------------------------------------
_FALLBACK_HEADER = (
    "# This data file generated by Converter",
    "#",
    "# rsid\tchromosome\tposition\tgenotype",
)


def extract_template_header(template_path: Path) -> list[str]:
    """
    Читает '#'-строки трафарета (как есть, без \\r\\n/\\n на конце) —
    эквивалент шага 9.4 гайда:
        grep '^#' template.txt | tr -d '\\r' > header.txt
    Файл Генотека принимает как "настоящий" 23andMe-экспорт именно по этому
    заголовку (20 строк с описанием формата), поэтому копировать его 1-в-1
    из трафарета, а не генерировать заново, — критично.
    """
    template_path = Path(template_path)
    header_lines: list[str] = []
    with template_path.open("r", encoding="utf-8-sig", newline="") as f:
        for line in f:
            if not line.startswith("#"):
                # Заголовок всегда идёт первым блоком строк в файле —
                # как только встретили не-'#' строку, заголовок кончился.
                break
            header_lines.append(line.rstrip("\r\n"))
    if not header_lines:
        raise AssemblyError(f"В трафарете {template_path} не найден '#'-заголовок")
    return header_lines


# ---------------------------------------------------------------------------
# Сборка финального файла
# ---------------------------------------------------------------------------
def assemble_final(skeleton: list[SkeletonRow], genotypes: dict[str, str],
                   output_path: Path, format_version: str = "v3",
                   template_path: Path | None = None) -> None:
    """
    template_path: путь к оригинальному трафарету (template.txt). Если
    указан, заголовок ('#'-строки) копируется из него дословно — так же,
    как это делает гайд (Часть 9.4). Если не указан (обратная совместимость),
    используется старый 3-строчный заглушечный заголовок — но именно он,
    предположительно, был причиной отказа Генотека, так как не совпадает
    с родным заголовком 23andMe-экспорта.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    line_ending = "\r\n" if format_version == "v5" else "\n"

    header_lines = (
        extract_template_header(template_path) if template_path is not None
        else list(_FALLBACK_HEADER)
    )

    with output_path.open("w", encoding="utf-8", newline="") as f:
        for header_line in header_lines:
            f.write(f"{header_line}{line_ending}")
        for row in skeleton:
            key = f"{row.chrom}_{row.pos}"
            genotype = genotypes.get(key, "--")
            f.write(f"{row.rsid}\t{row.chrom}\t{row.pos}\t{genotype}{line_ending}")

    logger.info("Собран финальный файл: %d строк в %s (формат %s, заголовок из %s)",
                len(skeleton), output_path, format_version,
                template_path if template_path is not None else "встроенного шаблона")


# ---------------------------------------------------------------------------
# 7 проверок из Части 10 гайда
# ---------------------------------------------------------------------------
def validate_output(output_path: Path, template_path: Path, format_version: str = "v3") -> ValidationResult:
    output_path = Path(output_path)
    template_path = Path(template_path)
    errors: list[str] = []

    with output_path.open("r", encoding="utf-8-sig", newline="") as f:
        output_lines = f.readlines()
    with template_path.open("r", encoding="utf-8-sig", newline="") as f:
        template_lines = f.readlines()

    output_data = [l for l in output_lines if not l.startswith("#")]
    template_data = [l for l in template_lines if not l.startswith("#")]

    total_lines = len(output_data)
    expected_lines = len(template_data)
    lines_match = total_lines == expected_lines
    if not lines_match:
        errors.append(f"Число строк не совпадает: {total_lines} vs {expected_lines}")

    fields_per_line: set[int] = set()
    for line in output_data:
        parts = line.rstrip("\r\n").split("\t")
        fields_per_line.add(len(parts))
    fields_valid = fields_per_line == {4}
    if not fields_valid:
        errors.append(f"Некорректное число полей: {fields_per_line}")

    crlf_count = sum(1 for l in output_data if l.endswith("\r\n"))
    crlf_expected = format_version == "v5"
    crlf_match = (crlf_count == total_lines) if crlf_expected else (crlf_count == 0)
    if not crlf_match:
        errors.append(f"CRLF/LF не соответствует формату {format_version}")

    chromosomes: list[str] = []
    for line in output_data:
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) >= 2:
            chromosomes.append(parts[1])

    total = len(output_data)
    missing = sum(1 for l in output_data if l.rstrip("\r\n").split("\t")[3] == "--")
    call_rate = 100.0 * (total - missing) / total if total > 0 else 0.0

    structure_identical = True
    for i, (out_line, tmpl_line) in enumerate(zip(output_data, template_data)):
        out_parts = out_line.rstrip("\r\n").split("\t")[:3]
        tmpl_parts = tmpl_line.rstrip("\r\n").split("\t")[:3]
        if out_parts != tmpl_parts:
            structure_identical = False
            errors.append(f"Структура не идентична на строке {i+1}")
            break

    positions_seen: set[tuple[str, int]] = set()
    template_duplicate_positions: list[tuple[str, int]] = []
    for line in output_data:
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) >= 3:
            key = (parts[1], int(parts[2]))
            if key in positions_seen:
                # НЕ добавляем в errors и НЕ прерываем цикл (break) —
                # это позиция, унаследованная из template_path (см.
                # докстринг ValidationResult.template_duplicate_positions
                # выше), а не признак сбоя сборки. Собираем ВСЕ такие
                # позиции для прозрачности, а не только первую.
                template_duplicate_positions.append(key)
                continue
            positions_seen.add(key)

    if template_duplicate_positions:
        logger.info(
            "В выходном файле %d раз(а) встречаются позиции, дублирующиеся "
            "уже в самом трафарете %s (не ошибка сборки — унаследованная "
            "особенность template_path): %s",
            len(template_duplicate_positions), template_path,
            template_duplicate_positions[:10],
        )

    return ValidationResult(
        total_lines=total_lines, expected_lines=expected_lines, lines_match=lines_match,
        fields_per_line=fields_per_line, fields_valid=fields_valid, crlf_count=crlf_count,
        crlf_expected=crlf_expected, crlf_match=crlf_match, chromosomes=chromosomes,
        call_rate=call_rate, structure_identical=structure_identical, errors=errors,
        template_duplicate_positions=template_duplicate_positions,
    )