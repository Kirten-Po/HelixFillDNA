"""
adapters/vcf_source.py
Адаптер для готового VCF (свой файл, экспорт другого сервиса, WGS и т.д.)
-> ParseResult. В отличие от ftdna_v3/myheritage_v5, здесь не нужен
референсный геном и не нужно разрешать ориентацию по чипу — REF/ALT/GT
уже даны в самом VCF, мы их только нормализуем и валидируем.

Интерфейс функции сделан совместимым с остальными парсерами
(parser_fn(path, reference, both_non_ref_threshold_pct=...)), чтобы
main.py мог вызывать её точно так же, как parse_ftdna_v3 /
parse_myheritage_v5 — параметр reference у VCF-источника просто
игнорируется (можно передавать None).
"""
from __future__ import annotations
import gzip
import hashlib
import json
import logging
from pathlib import Path
from typing import IO, Optional

from .base import ParsedVariant, ParseResult

logger = logging.getLogger(__name__)

ACGT = {"A", "C", "G", "T"}
CHROM_NORMALIZE = {
    "XY": "X", "M": "MT",
    "chrX": "X", "chrY": "Y", "chrM": "MT", "chrMT": "MT",
}

# Дефолт держим тем же именем, что и в build_vcf/остальном пайплайне,
# чтобы --format vcf можно было гонять через prepare-mis/assemble без
# лишних телодвижений.
DEFAULT_SAMPLE_NAME = "genotek"


class VCFFormatError(ValueError):
    pass


# ---------------------------------------------------------------------------
def _normalize_chrom(chrom: str) -> str:
    if chrom.startswith("chr") and chrom not in CHROM_NORMALIZE:
        chrom = chrom[3:]
    return CHROM_NORMALIZE.get(chrom, chrom)


def _open_text(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _find_sample_column(header_fields: list[str], sample_name: Optional[str]) -> tuple[int, str]:
    """
    header_fields — колонки строки '#CHROM\tPOS\t...\tFORMAT\tSAMPLE1\t...'.
    Возвращает (индекс_колонки_образца, имя_образца).
    Если sample_name не задан — берём первый (и обычно единственный)
    образец после FORMAT.
    """
    try:
        format_idx = header_fields.index("FORMAT")
    except ValueError as exc:
        raise VCFFormatError(
            "В заголовке VCF не найдена колонка FORMAT — это не похоже на "
            "валидный VCF с генотипами (нужен как минимум один образец)."
        ) from exc

    sample_cols = header_fields[format_idx + 1:]
    if not sample_cols:
        raise VCFFormatError("В VCF нет ни одной колонки с образцом (после FORMAT).")

    if sample_name is None:
        if len(sample_cols) > 1:
            logger.warning(
                "В VCF несколько образцов (%s) — sample_name не указан, "
                "беру первый: %s", ", ".join(sample_cols), sample_cols[0]
            )
        return format_idx + 1, sample_cols[0]

    if sample_name not in sample_cols:
        raise VCFFormatError(
            f"Образец '{sample_name}' не найден в VCF. "
            f"Доступные образцы: {', '.join(sample_cols)}"
        )
    return header_fields.index(sample_name), sample_name


# ---------------------------------------------------------------------------
def _chip_signature(positions: list[tuple[str, int]]) -> str:
    h = hashlib.sha256()
    for chrom, pos in positions:
        h.update(f"{chrom}:{pos}\n".encode("utf-8"))
    return h.hexdigest()[:16]


def save_position_cache(cache_dir: Path, result: ParseResult) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"{result.chip_signature}.positions.json"
    if not out_path.exists():
        payload = [(v.chrom, v.pos) for v in result.variants]
        out_path.write_text(json.dumps(payload), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
def parse_vcf_source(
    vcf_path: Path,
    reference=None,  # игнорируется; параметр только ради общего интерфейса
    both_non_ref_threshold_pct: float = 0.1,  # игнорируется, оставлен для совместимости
    sample_name: Optional[str] = None,
) -> ParseResult:
    """
    Разбирает готовый VCF (обычный или .gz) и возвращает ParseResult.

    Правила:
      - мультиаллельные позиции (ALT с запятой) и инделы (len(REF)!=1 или
        len(ALT)!=1) отбрасываются в invalid_codes;
      - './.' , '.' и любой GT без ровно двух аллелей — missing;
      - REF/ALT, не входящие в ACGT — ref_non_acgt;
      - GT нормализуется в 0/0 (alt='.'), 0/1, 1/1 — как в остальных
        адаптерах, чтобы дальше по пайплайну (build_vcf, assembler) всё
        работало без изменений.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF не найден: {vcf_path}")

    result = ParseResult()
    positions_for_signature: list[tuple[str, int]] = []

    with _open_text(vcf_path) as f:
        header_fields: Optional[list[str]] = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\r\n").lstrip("#").split("\t")
                break
        if header_fields is None:
            raise VCFFormatError("В файле не найдена строка заголовка '#CHROM...' — это не VCF.")

        sample_idx, resolved_sample_name = _find_sample_column(header_fields, sample_name)
        logger.info("VCF-источник: беру генотипы образца '%s'", resolved_sample_name)

        for raw_line in f:
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) <= sample_idx:
                result.malformed_rows += 1
                continue

            chrom, pos_str, rsid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            format_col = fields[8] if len(fields) > 8 else ""
            sample_col = fields[sample_idx]

            if "," in alt:
                result.invalid_codes += 1
                continue
            if len(ref) != 1 or len(alt) != 1:
                # индели вне контракта ParsedVariant (только SNP)
                result.invalid_codes += 1
                continue

            ref = ref.upper()
            alt = alt.upper()

            format_keys = format_col.split(":") if format_col else ["GT"]
            sample_values = sample_col.split(":")
            try:
                gt_idx = format_keys.index("GT")
                gt_raw = sample_values[gt_idx]
            except (ValueError, IndexError):
                result.malformed_rows += 1
                continue

            gt_norm = gt_raw.replace("|", "/")
            if gt_norm in (".", "./.", ".|."):
                result.missing += 1
                continue

            alleles = gt_norm.split("/")
            if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
                # не-биаллельный / multiallelic GT (2, 3...) или мусор
                result.total_measured += 1
                result.invalid_codes += 1
                continue

            result.total_measured += 1

            if ref not in ACGT or alt not in ACGT:
                result.ref_non_acgt += 1
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                result.malformed_rows += 1
                result.total_measured -= 1
                continue

            chrom_norm = _normalize_chrom(chrom)
            a1, a2 = alleles
            if a1 == "0" and a2 == "0":
                resolved_ref, resolved_alt, gt = ref, ".", "0/0"
            elif a1 == a2 == "1":
                resolved_ref, resolved_alt, gt = ref, alt, "1/1"
            else:
                resolved_ref, resolved_alt, gt = ref, alt, "0/1"

            final_rsid = rsid if rsid and rsid != "." else f"{chrom_norm}_{pos}"

            result.variants.append(ParsedVariant(
                rsid=final_rsid, chrom=chrom_norm, pos=pos,
                ref=resolved_ref, alt=resolved_alt, gt=gt,
            ))
            positions_for_signature.append((chrom_norm, pos))

    if result.total_measured == 0:
        raise VCFFormatError("В VCF не найдено ни одной пригодной для использования позиции")

    positions_for_signature.sort()
    result.chip_signature = _chip_signature(positions_for_signature)
    logger.info(
        "VCF-источник: измерено=%d, пропущено=%d, invalid=%d, "
        "ref_non_acgt=%d, битых_строк=%d, годных=%d, signature=%s",
        result.total_measured, result.missing, result.invalid_codes,
        result.ref_non_acgt, result.malformed_rows, len(result.variants),
        result.chip_signature,
    )
    return result
