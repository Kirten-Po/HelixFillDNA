"""
core/liftover.py
Промт "HRC / TopMed", лифтовер координат: чистый Python перенос координат
одной точечной позиции (SNP) между сборками генома (GRCh37 -> GRCh38 и
обратно) по UCSC chain-файлу — без внешних зависимостей (liftOver/CrossMap/
`bcftools +liftover` не требуются).

Почему не bcftools +liftover / Picard LiftoverVcf:
  - `bcftools +liftover` — не входит в основной репозиторий bcftools (это
    сторонний плагин, github.com/freeseek/score), его наличие в конкретном
    бандле htslib-бинарников проекта не гарантировано и не проверено.
  - Picard LiftoverVcf требует Java, которой в проекте сейчас нет вообще
    как зависимости.
  - Наши данные — ТОЛЬКО SNP (одна позиция, один референсный/альтернативный
    аллель; adapters/*.py в принципе не порождают инделы, ParsedVariant.ref/
    alt всегда одна буква ACGT или "."), поэтому не нужна вся сложность
    Picard/CrossMap с многоаллельными сайтами, VCF-заголовками, инделами и
    сдвигом REF на границах блоков — достаточно перенести (chrom, 1-based
    pos) в целевую сборку по тем же alignment-блокам, что использует
    оригинальный liftOver.

Формат UCSC chain-файла (https://genome.ucsc.edu/goldenPath/help/chain.html):

    chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    size dt dq
    size dt dq
    ...
    size

- Одна "chain"-запись описывает один блок выравнивания между исходной (t,
  target — в наших chain-файлах это ВСЕГДА исходная сборка, например hg19
  для hg19ToHg38.over.chain.gz) и целевой (q, query) последовательностями.
- Внутри chain-записи идёт список ungapped-блоков ("size") — участков без
  вставок/делеций между t и q, разделённых величиной сдвига в t (dt) и в q
  (dq) до начала следующего блока. Последняя строка блока содержит только
  "size" (без dt/dq) — она завершает chain-запись, после неё в файле идёт
  пустая строка перед следующей "chain".
- tStrand по соглашению UCSC для liftOver-chain-файлов ВСЕГДА "+" — если
  это не так, что-то нестандартное в файле, и код это явно проверяет.
- qStrand может быть "+" или "-". Если "-", то qStart/qEnd/координаты блоков
  в файле заданы относительно РЕВЕРС-КОМПЛЕМЕНТА query-последовательности
  (т.е. отсчитаны от конца исходной plus-strand хромосомы), а не от её
  начала — при переносе координаты в plus-strand систему координат целевой
  сборки нужно явно развернуть: plus_pos = qSize - qPos_reverse - 1.

Единицы координат: UCSC chain-файл — 0-based, полуоткрытые интервалы
[start, end). ParsedVariant.pos (adapters/*.py) — 1-based. ChainLiftover.lift()
принимает и возвращает 1-based позицию — конвертация в/из 0-based происходит
внутри метода, вызывающему коду (adapters/*.py) не нужно об этом думать.

Множественные перекрывающиеся chain-записи для одной и той же t-хромосомы
(бывает, если в файле несколько chain-блоков с разным score покрывают
пересекающиеся участки) НЕ разрешаются здесь через выбор "лучшего" score —
используется первый найденный по позиции блок среди объединённого
отсортированного списка блоков всех chain-записей этой хромосомы. Для
подавляющего большинства позиций UCSC chain-файлы содержат ровно один
покрывающий блок, так что это приближение (тот же принцип, которым
пользуется свободная библиотека pyliftover) даёт результат, идентичный
полноценному liftOver, за редким исключением сильно перестроенных участков
генома — что для точечных позиций чипа ДНК-теста статистически ничтожно.

=============================================================================
ИЗМЕНЕНИЯ (промт "TopMed/HRC v4", баг "Неподдерживаемая хромосома:
1_KI270766v1_alt"):
=============================================================================

Найденный баг: официальные UCSC chain-файлы (hg19ToHg38.over.chain.gz)
содержат chain-записи, где ЦЕЛЕВОЙ (qName) контиг — не только основная
хромосома GRCh38 ("chr1" и т.п.), но и альтернативные/случайные контиги
той же сборки: "chr1_KI270766v1_alt", "chr1_KI270706v1_random" и другие —
это официальная часть сборки GRCh38, представляющая альтернативные
гаплотипы участков, уже покрытых основной хромосомой. Раньше _parse()
складывал ВСЕ блоки в индекс без разбора, к какому контигу они ведут, а
_canonical_chrom() только снимала префикс "chr" и нормализовала "M"/"MT",
не проверяя, что результат — один из канонических 1-22/X/Y/MT. В
результате ChainLiftover.lift() иногда честно находил по bisect ближайший
подходящий блок — но этот блок вёл на альт-контиг — и возвращал
("1_KI270766v1_alt", pos). Это попадало в ParsedVariant.chrom, и
core/pure_python_core.py::validate_variants() затем справедливо отклоняло
всю сборку VCF с ошибкой "Неподдерживаемая хромосома: 1_KI270766v1_alt".

Исправление: блоки, чей qName после канонизации НЕ входит в множество
_CANONICAL_TARGET_CHROMS (1-22, X, Y, MT), теперь ПОЛНОСТЬЮ исключаются
из индекса ещё на этапе _parse() — они никогда не участвуют в bisect и
никогда не могут быть возвращены lift(). Позиция, которая раньше
ошибочно переносилась бы на такой блок, теперь либо:
  - находит другой, канонический блок, покрывающий ту же позицию t —
    переносится на него как обычно (типичный случай — альт-контиг и
    основная хромосома часто перекрываются по t-диапазону), либо
  - не находит ни одного канонического блока на этом участке —
    lift() возвращает None (in_gap), позиция учитывается в
    result.lift_failed в адаптерах, как и любая другая нелифтуемая
    позиция — что и является корректным поведением: раз каноничная
    хромосома этот участок не покрывает, переносить позицию НЕКУДА.

Ничего в самом lift() менять не потребовалось — раз "плохие" блоки просто
не попадают в self._blocks, вся остальная логика (bisect, границы блока,
reverse-strand математика) продолжает работать без изменений.
"""
from __future__ import annotations

import bisect
import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class LiftoverError(RuntimeError):
    """Chain-файл скачался/присутствует на диске, но не парсится, либо
    имеет структуру, которую этот модуль сознательно не поддерживает
    (tStrand != '+', синтаксически некорректная строка блока и т.п.)."""


@dataclass
class LiftoverStats:
    """Накопительная статистика вызовов ChainLiftover.lift() за время
    жизни объекта — для диагностики (например, в будущем можно вывести в
    лог по завершении парсинга файла, аналогично result.lift_failed в
    adapters/*.py, который считается независимо самим адаптером)."""
    total_calls: int = 0
    lifted: int = 0
    no_chain_for_chrom: int = 0   # для этой хромосомы нет ни одного блока
    in_gap: int = 0               # позиция между alignment-блоками
    out_of_target_range: int = 0  # результат вне [0, qSize) целевой сборки

    @property
    def failed(self) -> int:
        return self.no_chain_for_chrom + self.in_gap + self.out_of_target_range


@dataclass(frozen=True)
class _Block:
    t_start: int   # 0-based, включительно
    t_end: int     # 0-based, исключительно
    q_base: int    # см. ChainLiftover._parse() — уже учитывает qStrand
    is_reverse: bool
    q_size: int
    q_chrom: str   # каноническое имя целевой хромосомы ("1".."22","X","Y","MT")


def _canonical_chrom(name: str) -> str:
    """
    Приводит имя хромосомы UCSC chain-файла (например "chr1", "chrX",
    "chrM", "chr1_KI270706v1_random" — альтернативные контиги тоже
    встречаются в chain-файлах) к тому же каноническому виду, что
    используют adapters/*.py::_normalize_chrom() — без префикса "chr",
    "M"/"MT" -> "MT".

    ⚠ Эта функция НЕ проверяет, что результат входит в канонический набор
    1-22/X/Y/MT — она просто снимает префикс/нормализует митохондрию.
    Для "chr1_KI270766v1_alt" она вернёт "1_KI270766v1_alt" как есть.
    Отсеивание нехромосомных/альтернативных контигов — отдельный шаг,
    см. _CANONICAL_TARGET_CHROMS и его использование в _parse().
    """
    body = name[3:] if name.lower().startswith("chr") else name
    if body in ("M", "MT"):
        return "MT"
    return body


# Промт "TopMed/HRC v4": множество ЕДИНСТВЕННО допустимых канонических
# целевых хромосом — то же самое, что принимает validate_variants() в
# core/pure_python_core.py и что производят _normalize_chrom() в
# adapters/ftdna_v3.py/adapters/myheritage_v5.py. Блоки chain-файла,
# чей канонизированный qName не входит сюда (альт-контиги вида
# "1_KI270766v1_alt", случайные контиги вида "1_KI270706v1_random" и
# любые другие нестандартные имена), исключаются из индекса ещё на
# этапе _parse() — см. докстринг модуля, раздел про баг
# "Неподдерживаемая хромосома".
_CANONICAL_TARGET_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def _open_chain_text(path: Path):
    """
    Открывает chain-файл в текстовом режиме. UCSC раздаёт chain-файлы
    сжатыми gzip (.chain.gz) — определяем формат по магическим байтам
    файла, а не по расширению, чтобы синтетические/тестовые chain-файлы
    (обычный текст, без .gz) тоже читались без специальных ухищрений в
    тестах.
    """
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


class ChainLiftover:
    """
    Переносит SNP-позиции между сборками генома по UCSC chain-файлу.
    Парсинг файла и построение bisect-индекса происходят ОДИН раз в
    конструкторе — объект рассчитан на переиспользование для всех
    ~700-900 тыс. позиций одного чипа за один вызов parse_ftdna_v3()/
    parse_myheritage_v5() (повторный парсинг файла на каждую позицию не
    происходит.
    """

    def __init__(self, chain_path: Path):
        self.chain_path = Path(chain_path)
        if not self.chain_path.exists():
            raise LiftoverError(f"Chain-файл не найден: {self.chain_path}")

        # {канонический tName: [_Block, ...]} — отсортировано по t_start
        # после парсинга, плюс параллельный список t_start для bisect.
        self._blocks: dict[str, list[_Block]] = {}
        self._starts: dict[str, list[int]] = {}
        self.stats = LiftoverStats()
        # Промт "TopMed/HRC v4": сколько блоков было пропущено при парсинге
        # как ведущих на альт/random-контиги целевой сборки — только для
        # информационного лога после _parse(), на bisect/lift() не влияет.
        self._skipped_alt_blocks = 0

        self._parse()

    # -----------------------------------------------------------------
    def _parse(self) -> None:
        raw_blocks: dict[str, list[_Block]] = {}
        skipped_alt_blocks = 0

        with _open_chain_text(self.chain_path) as f:
            in_chain = False
            t_chrom = ""
            t_pos = 0
            q_pos = 0
            q_size = 0
            q_chrom = ""
            is_reverse = False

            for line_no, raw_line in enumerate(f, start=1):
                line = raw_line.strip()

                if not line:
                    # Пустая строка — конец текущей chain-записи.
                    in_chain = False
                    continue

                if line.startswith("chain"):
                    parts = line.split()
                    if len(parts) != 13:
                        raise LiftoverError(
                            f"{self.chain_path.name}:{line_no}: некорректная "
                            f"строка заголовка chain (ожидалось 13 полей, "
                            f"получено {len(parts)}): {line!r}"
                        )
                    (_kw, _score, tName, _tSize, tStrand, tStart, _tEnd,
                     qName, qSize, qStrand, qStart, _qEnd, _chain_id) = parts

                    if tStrand != "+":
                        raise LiftoverError(
                            f"{self.chain_path.name}:{line_no}: tStrand={tStrand!r} "
                            f"не поддерживается (ожидается '+' — так всегда "
                            f"оформлены официальные UCSC liftOver chain-файлы; "
                            f"нестандартный tStrand означает, что это не "
                            f"обычный *.over.chain.gz с сайта UCSC)."
                        )
                    if qStrand not in ("+", "-"):
                        raise LiftoverError(
                            f"{self.chain_path.name}:{line_no}: некорректный "
                            f"qStrand={qStrand!r} (ожидается '+' или '-')."
                        )

                    try:
                        t_pos = int(tStart)
                        q_pos = int(qStart)
                        q_size = int(qSize)
                    except ValueError as exc:
                        raise LiftoverError(
                            f"{self.chain_path.name}:{line_no}: нечисловые "
                            f"координаты в заголовке chain: {line!r}"
                        ) from exc

                    t_chrom = _canonical_chrom(tName)
                    q_chrom = _canonical_chrom(qName)
                    is_reverse = (qStrand == "-")
                    in_chain = True
                    continue

                if not in_chain:
                    # Строка вне chain-записи и не заголовок chain — по
                    # спецификации такого быть не должно (единственное, что
                    # тут ожидаемо, — комментарии, но формат chain их не
                    # определяет; чтобы не падать на возможных пустых
                    # хвостах/BOM, просто пропускаем такие строки, ничего
                    # не накапливая).
                    continue

                fields = line.split()
                if len(fields) not in (1, 3):
                    raise LiftoverError(
                        f"{self.chain_path.name}:{line_no}: некорректная "
                        f"строка данных блока (ожидалось 1 или 3 поля, "
                        f"получено {len(fields)}): {line!r}"
                    )
                try:
                    size = int(fields[0])
                except ValueError as exc:
                    raise LiftoverError(
                        f"{self.chain_path.name}:{line_no}: нечисловой размер "
                        f"блока: {line!r}"
                    ) from exc

                if size > 0:
                    if is_reverse:
                        q_base = q_size - q_pos - 1
                    else:
                        q_base = q_pos

                    # Промт "TopMed/HRC v4": блоки, ведущие на альт/random-
                    # контиги целевой сборки (q_chrom не входит в
                    # _CANONICAL_TARGET_CHROMS — например
                    # "1_KI270766v1_alt", "1_KI270706v1_random"), НЕ
                    # добавляются в индекс вовсе. См. докстринг модуля.
                    if q_chrom in _CANONICAL_TARGET_CHROMS:
                        block = _Block(
                            t_start=t_pos, t_end=t_pos + size,
                            q_base=q_base, is_reverse=is_reverse,
                            q_size=q_size, q_chrom=q_chrom,
                        )
                        raw_blocks.setdefault(t_chrom, []).append(block)
                    else:
                        skipped_alt_blocks += 1

                if len(fields) == 3:
                    try:
                        dt = int(fields[1])
                        dq = int(fields[2])
                    except ValueError as exc:
                        raise LiftoverError(
                            f"{self.chain_path.name}:{line_no}: нечисловой "
                            f"сдвиг (dt/dq) в строке блока: {line!r}"
                        ) from exc
                    t_pos += size + dt
                    q_pos += size + dq
                else:
                    # Последняя строка блока в chain-записи — сдвигов нет,
                    # запись логически завершена (пустая строка ниже это
                    # только подтвердит).
                    t_pos += size
                    in_chain = False

        if not raw_blocks:
            raise LiftoverError(
                f"В chain-файле {self.chain_path} не найдено ни одной "
                f"валидной chain-записи с блоками выравнивания, ведущими на "
                f"канонические хромосомы (1-22/X/Y/MT) — файл пуст, "
                f"повреждён, либо содержит только блоки на альтернативные/"
                f"случайные контиги."
            )

        for chrom, blocks in raw_blocks.items():
            blocks.sort(key=lambda b: b.t_start)
            self._blocks[chrom] = blocks
            self._starts[chrom] = [b.t_start for b in blocks]

        self._skipped_alt_blocks = skipped_alt_blocks
        logger.info(
            "ChainLiftover: %s — загружено %d хромосом(ы), %d блоков "
            "выравнивания всего%s",
            self.chain_path.name, len(self._blocks),
            sum(len(v) for v in self._blocks.values()),
            (
                f" (пропущено {skipped_alt_blocks} блок(ов), ведущих на "
                f"альтернативные/случайные контиги целевой сборки — "
                f"позиции, которые перенеслись бы только на них, будут "
                f"считаться нелифтуемыми)"
                if skipped_alt_blocks else ""
            ),
        )

    # -----------------------------------------------------------------
    def lift(self, chrom: str, pos: int) -> Optional[tuple[str, int]]:
        """
        Переносит 1-based позицию (chrom, pos) в целевую сборку по
        chain-файлу. Возвращает (target_chrom, target_pos) — тоже
        1-based, канонический вид хромосомы (без "chr"), ВСЕГДА из
        множества 1-22/X/Y/MT (см. _CANONICAL_TARGET_CHROMS — блоки,
        ведущие на альт/random-контиги, отфильтрованы ещё в _parse() и
        физически не могут быть возвращены отсюда). Возвращает None,
        если перенос невозможен (нет chain для этой хромосомы, позиция
        попадает в разрыв между блоками выравнивания — в том числе если
        единственный покрывающий блок вёл на альт/random-контиг и был
        отфильтрован, — либо результат оказался бы вне границ целевой
        хромосомы).
        """
        self.stats.total_calls += 1

        blocks = self._blocks.get(chrom)
        if not blocks:
            self.stats.no_chain_for_chrom += 1
            return None

        t0 = pos - 1  # 1-based -> 0-based
        starts = self._starts[chrom]
        idx = bisect.bisect_right(starts, t0) - 1
        if idx < 0:
            self.stats.in_gap += 1
            return None

        block = blocks[idx]
        if not (block.t_start <= t0 < block.t_end):
            # Позиция между двумя блоками (в разрыве/делеции/инсерции
            # относительно целевой сборки, либо единственный покрывающий
            # блок вёл на альт/random-контиг и был отфильтрован в
            # _parse()) — не переносится.
            self.stats.in_gap += 1
            return None

        offset = t0 - block.t_start
        if block.is_reverse:
            q0 = block.q_base - offset
        else:
            q0 = block.q_base + offset

        if q0 < 0 or q0 >= block.q_size:
            self.stats.out_of_target_range += 1
            return None

        self.stats.lifted += 1
        return block.q_chrom, q0 + 1  # 0-based -> 1-based
