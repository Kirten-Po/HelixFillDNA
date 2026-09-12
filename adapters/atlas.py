"""
adapters/atlas.py
Адаптер файлов компании «Атлас» -> ParseResult.

Файл Атласа оформлен точно так же, как сырые данные 23andMe: одна
'#'-строка с названиями колонок и четыре колонки через табуляцию
(rsid, chromosome, position, genotype), хромосомы как X/Y/MT, гаплоидные
вызовы одной буквой, инделы как I/D. Ровно этот формат уже читает
adapters/ancestry_v2.py (его LAYOUT_CONVERTED — «файл в оформлении
23andMe»): те же 4 колонки, то же удвоение гаплоидной буквы для X/Y/MT,
та же отбраковка инделов в QC. Поэтому здесь не дублируется разбор, а
переиспользуется готовый парсер — иначе появилась бы вторая копия логики
разрешения ориентации аллелей, которую пришлось бы править дважды.

Отличие Атласа — не в оформлении, а в СБОРКЕ ГЕНОМА: координаты в GRCh38.
Приводит их к GRCh37 отдельный Этап 0 (core/atlas_convert.py), до вызова
этого парсера; здесь на входе уже GRCh37, как у любого другого источника.
Параметр liftover остаётся в сигнатуре и означает то же, что у остальных
адаптеров: ПРЯМОЙ перенос GRCh37 -> сборка панели для panel="topmed".
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from .ancestry_v2 import (
    DEFAULT_BOTH_NON_REF_THRESHOLD_PCT,
    ReferenceGenome,
    parse_ancestry_v2,
    save_position_cache as _save_position_cache,
    save_position_cache_broad as _save_position_cache_broad,
)
from .base import ParseResult
from core.liftover import ChainLiftover

logger = logging.getLogger(__name__)


class AtlasFormatError(ValueError):
    pass


def parse_atlas(
    csv_path: Path,
    reference: ReferenceGenome,
    both_non_ref_threshold_pct: float = DEFAULT_BOTH_NON_REF_THRESHOLD_PCT,
    liftover: Optional[ChainLiftover] = None,
) -> ParseResult:
    """
    Парсит файл Атласа (после Этапа 0 — уже в GRCh37, в оформлении
    23andMe). Сигнатура и семантика совпадают с parse_ftdna_v3()/
    parse_myheritage_v5()/parse_ancestry_v2().
    """
    logger.info("Атлас: разбор файла %s (оформление 23andMe)", Path(csv_path).name)
    return parse_ancestry_v2(
        csv_path, reference,
        both_non_ref_threshold_pct=both_non_ref_threshold_pct,
        liftover=liftover,
    )


# Кэш позиций/сигнатура чипа считаются ровно так же, как у остальных
# источников — единый интерфейс main.py::SOURCES[...].
save_position_cache = _save_position_cache
save_position_cache_broad = _save_position_cache_broad
