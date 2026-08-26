"""
template/skeleton.py
Извлекает позиции из трафарета (template.txt) — структуру файла без генотипов.
Эквивалент шага 9.1 гайда:
grep -v '^#' template.txt | awk ... > panel_pos.txt
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SkeletonRow:
    rsid: str
    chrom: str
    pos: int
    raw_line: str  # исходная строка из трафарета (с \r\n или \n)


class SkeletonError(ValueError):
    pass


def extract_skeleton(
    template_path: Path,
    autosomes_only: bool = True,
) -> list[SkeletonRow]:
    """
    Извлекает позиции из трафарета (template.txt).
    template_path: путь к файлу template.txt (формат 23andMe v3 или v5).
    autosomes_only: если True, берёт только хромосомы 1-22 (как в гайде).
    Возвращает список SkeletonRow в порядке трафарета.
    """
    template_path = Path(template_path)
    if not template_path.exists():
        raise FileNotFoundError(f"Трафарет не найден: {template_path}")

    skeleton: list[SkeletonRow] = []
    allowed_chroms = {str(i) for i in range(1, 23)} if autosomes_only else None

    with template_path.open("r", encoding="utf-8-sig", newline="") as f:
        for line in f:
            # Пропускаем заголовок
            if line.startswith("#"):
                continue
            # Убираем \r\n или \n для парсинга
            stripped = line.rstrip("\r\n")
            parts = stripped.split("\t")
            if len(parts) < 3:
                logger.warning("Строка трафарета имеет меньше 3 полей: %r", stripped)
                continue
            rsid, chrom, pos_str = parts[0], parts[1], parts[2]
            # Фильтруем по хромосомам
            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                logger.warning("Некорректная позиция в трафарете: %r", pos_str)
                continue
            skeleton.append(SkeletonRow(
                rsid=rsid,
                chrom=chrom,
                pos=pos,
                raw_line=line,
            ))

    if not skeleton:
        raise SkeletonError(
            f"В трафарете {template_path} не найдено ни одной позиции"
        )

    logger.info(
        "Извлечено %d позиций из трафарета %s",
        len(skeleton), template_path,
    )
    return skeleton


def save_skeleton(
    skeleton: list[SkeletonRow],
    output_path: Path,
    line_ending: str = "\n",
) -> None:
    """
    Сохраняет skeleton в файл (только первые 3 колонки, без генотипа).
    Используется для промежуточной проверки.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as f:
        for row in skeleton:
            f.write(f"{row.rsid}\t{row.chrom}\t{row.pos}{line_ending}")
    logger.info("Сохранён skeleton: %d строк в %s", len(skeleton), output_path)