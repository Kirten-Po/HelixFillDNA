"""
core/panel_advisor.py
Рекомендация референсной панели (HRC r1.1 против TopMed r3) по СОСТАВУ
ЧИПА, а не по качеству файла.

Почему это отдельно от core/preflight.py
----------------------------------------
"Плохое качество -> берём TopMed" не работает как правило. Если файл
действительно плохой (низкий call rate чипа, контаминация, перепутанный
пол), большая панель его не спасёт — она аккуратно доимпутирует мусор.
Такой файл надо помечать и предупреждать, чем и занимается preflight.

Панель выбирается по другим признакам:

                    HRC r1.1                TopMed r3
    Сборка          GRCh37, лифтовер не     GRCh38, нужен лифтовер
                    нужен                   туда и обратно
    Образцов        ~32 тыс., почти все     ~97 тыс., разнообразная
                    европейцы
    Сайтов          ~39 млн                 ~300 млн
    Редкий хвост    обрезан по MAC>=5       глубокий

То есть TopMed стоит брать, когда образец неевропейский ИЛИ когда чип
насыщен редкими маркерами. Ни то ни другое не имеет отношения к
"качеству". Цена — лифтовер: он в пайплайне есть (core/liftover.py::
ChainLiftover), но это лишний шаг с лишним риском.

Что считает этот модуль
-----------------------
Один сигнал из двух, и оба измеримы тем, что в приложении уже есть, БЕЗ
единой закачки — донорские VCF (donors/<source>/<panel>/kgp_sub_*.vcf.gz)
уже лежат в кэше и несут в INFO частоты 1000 Genomes по суперпопуляциям:

  * доля маркеров чипа с EUR_AF ниже 0,5 % — "редкий хвост". На разборе
    падения 31.08 у FTDNA это было 13 %, у AncestryDNA — 3 %;
  * доля позиций чипа, для которых в 1000G phase3 записи нет вовсе (у
    FTDNA — 18 837 из 724 937).

Второй сигнал из исходной постановки — грубая оценка ПРОИСХОЖДЕНИЯ
образца проекцией на частоты суперпопуляций — здесь сознательно НЕ
реализован: он требует файлов частот всех суперпопуляций, а не только
того, что уже лежит в кэше. Место под него оставлено явно
(Recommendation.ancestry_hint), и правило ниже сформулировано так, чтобы
происхождение можно было добавить как второе слагаемое, ничего не
переписывая.

Как это встроено
----------------
Не молчаливым переключением. Экран с отчётом и строкой вида "рекомендую
TopMed: доля маркеров с MAF<0,5 % — 14 %", радиокнопка предвыбрана, но
переключается вручную. Молчаливый автовыбор плох тем, что потом нельзя
разобрать, почему два запуска разошлись.

Пороги
------
Их придётся калибровать самому: по двум файлам границу не вывести.
Значения ниже — осознанно грубые ориентиры; настоящая граница
нарисуется по накопленному runs_metrics.csv (core/metrics_log.py), где
рядом лежат метрики исходника, панель, Rsq, финальный call rate и
вердикт приёмки.
"""
from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

logger = logging.getLogger(__name__)

#: Порог "редкого" маркера: MAF в европейской подвыборке 1000 Genomes.
RARE_MAF = 0.005

#: Доля редких маркеров, выше которой имеет смысл платить за TopMed
#: лифтовером. ПРЕДВАРИТЕЛЬНЫЙ ориентир, посаженный между двумя
#: известными точками (AncestryDNA ~3 %, FTDNA ~13 %) ближе к верхней:
#: ошибиться в сторону HRC дешевле — он быстрее и без лифтовера.
RARE_TAIL_PCT_THRESHOLD = 10.0

#: Доля позиций чипа, которых в 1000G phase3 нет вовсе. ВСПОМОГАТЕЛЬНЫЙ
#: сигнал, а не самостоятельный триггер: отсутствие маркера в 1000G ещё
#: не значит, что он есть в TopMed, а HRC частично построен на том же
#: 1000G. Замеры: FTDNA 2,7 %, AncestryDNA 9,9 % — при этом редкий хвост
#: у них 13,4 % и 3,3 % соответственно, то есть эти две метрики
#: расходятся, и решать по второй было бы прямо неверно. Поэтому высокая
#: доля отсутствующих ЛИШЬ УСИЛИВАЕТ уже принятое по редкому хвосту
#: решение, а сама панель по ней не переключается.
MISSING_PCT_THRESHOLD = 8.0

PANEL_HRC = "hrc"
PANEL_TOPMED = "topmed"


class PanelAdvisorError(RuntimeError):
    """Посчитать состав чипа не по чему (нет кэша доноров)."""


# ---------------------------------------------------------------------------
# Состав чипа
# ---------------------------------------------------------------------------
@dataclass
class ChipComposition:
    donors_dir: str = ""
    chip_positions: int = 0        # всего позиций чипа в списке фильтра
    matched: int = 0               # из них найдено в 1000G phase3
    rare: int = 0                  # из найденных: MAF(EUR) < RARE_MAF
    monomorphic_eur: int = 0       # EUR_AF ровно 0 или 1 — для HRC мёртвый груз
    per_chrom_matched: dict[str, int] = field(default_factory=dict)

    @property
    def missing(self) -> int:
        return max(0, self.chip_positions - self.matched)

    @property
    def missing_pct(self) -> float:
        if not self.chip_positions:
            return 0.0
        return 100.0 * self.missing / self.chip_positions

    @property
    def rare_pct(self) -> float:
        """Доля редких — от НАЙДЕННЫХ в 1000G, а не от всех позиций чипа."""
        return 100.0 * self.rare / self.matched if self.matched else 0.0

    @property
    def monomorphic_eur_pct(self) -> float:
        return 100.0 * self.monomorphic_eur / self.matched if self.matched else 0.0


def find_donor_cache(donors_root: Path, source: str) -> Optional[Path]:
    """
    Папка с донорскими VCF для этого источника — любая панель, какая
    скачана. Частоты 1000 Genomes от выбора ПАНЕЛИ не зависят, поэтому
    считать состав чипа можно по тому кэшу, который уже есть, и до того,
    как пользователь выбрал панель для нового запуска (иначе получилась
    бы курица и яйцо: чтобы посоветовать панель, надо скачать доноров для
    панели).
    """
    base = Path(donors_root) / source
    if not base.is_dir():
        return None
    candidates = [d for d in sorted(base.iterdir())
                  if d.is_dir() and (d / "kgp_sub_1.vcf.gz").is_file()]
    return candidates[0] if candidates else None


def _count_chip_positions(donors_dir: Path) -> int:
    """
    Число позиций чипа — из списка, которым фильтровались доноры
    (`*_pos.txt`, две колонки: хромосома и позиция). Это ровно тот
    знаменатель, о котором идёт речь в "18 837 из 724 937".
    """
    for name in ("ftdna_pos.txt", "chip_pos.txt"):
        f = donors_dir / name
        if f.is_file():
            with open(f, "r", encoding="utf-8", errors="replace") as fh:
                return sum(1 for line in fh if line.strip())
    matches = sorted(donors_dir.glob("*_pos.txt"))
    if matches:
        with open(matches[0], "r", encoding="utf-8", errors="replace") as fh:
            return sum(1 for line in fh if line.strip())
    return 0


def _iter_eur_af_bcftools(vcf: Path, bcftools: str) -> Iterator[str]:
    """
    Быстрый путь: bcftools query отдаёт только нужное поле, не разворачивая
    строку с двумя тысячами генотипов. На кэше FTDNA (23 файла, 52 МБ)
    это разница в разы по сравнению с чтением gzip средствами Python.
    """
    proc = subprocess.run(
        [bcftools, "query", "-f", "%INFO/EUR_AF\n", str(vcf)],
        capture_output=True, text=True, errors="replace",
    )
    if proc.returncode != 0:
        raise PanelAdvisorError(
            f"bcftools query по {vcf.name} завершился с кодом "
            f"{proc.returncode}: {proc.stderr.strip()[:200]}"
        )
    for line in proc.stdout.splitlines():
        yield line.strip()


def _iter_eur_af_python(vcf: Path) -> Iterator[str]:
    """Запасной путь без bcftools — разбираем INFO сами."""
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 8:
                continue
            value = "."
            for field_ in parts[7].split(";"):
                if field_.startswith("EUR_AF="):
                    value = field_[7:]
                    break
            yield value


def _maf(value: str) -> Optional[float]:
    """
    EUR_AF -> минорная частота. Мультиаллельные записи приходят списком
    ("0.01,0.002") — берём максимальную альтернативную частоту: именно
    она определяет, насколько сайт вообще полиморфен в европейцах.
    Отсутствующее значение ('.') -> None, такие записи в статистику
    редкости не идут.
    """
    if not value or value == ".":
        return None
    best: Optional[float] = None
    for chunk in value.split(","):
        try:
            af = float(chunk)
        except ValueError:
            continue
        if best is None or af > best:
            best = af
    if best is None:
        return None
    return min(best, 1.0 - best)


def analyse_chip(
    donors_dir: Path,
    *,
    bcftools_path: Optional[str] = None,
    chip_positions: Optional[int] = None,
) -> ChipComposition:
    """
    Считает состав чипа по уже скачанному кэшу доноров. Ничего не качает.

    chip_positions — знаменатель (сколько всего позиций у чипа). По
    умолчанию берётся из списка фильтра рядом с донорами; передать явно
    имеет смысл, когда пре-флайт уже посчитал число маркеров исходника.
    """
    donors_dir = Path(donors_dir)
    vcfs = sorted(donors_dir.glob("kgp_sub_*.vcf.gz"))
    if not vcfs:
        raise PanelAdvisorError(
            f"В {donors_dir} нет донорских VCF — состав чипа считать не по "
            f"чему. Он посчитается сам после первого скачивания доноров."
        )

    comp = ChipComposition(donors_dir=str(donors_dir))
    comp.chip_positions = chip_positions or _count_chip_positions(donors_dir)

    bcftools = bcftools_path or shutil.which("bcftools")
    for vcf in vcfs:
        chrom = vcf.stem.replace("kgp_sub_", "").replace(".vcf", "")

        # Читатели пробуются по очереди: сначала быстрый (bcftools query
        # отдаёт одно поле, не разворачивая строку с двумя тысячами
        # генотипов), при любой его неудаче — свой разбор gzip.
        #
        # Отказ bcftools здесь не экзотика: его может не быть на машине,
        # он может не запуститься, а может и отказаться от конкретного
        # файла, если в заголовке нет описания тега EUR_AF. И ни один из
        # этих случаев не повод выбросить файл из статистики: тогда молча
        # уехал бы вниз знаменатель "найдено в 1000G", а за ним и доля
        # отсутствующих позиций — то есть отчёт показал бы неправду вместо
        # того, чтобы просто читать помедленнее.
        readers = []
        if bcftools:
            readers.append(("bcftools", lambda v=vcf: _iter_eur_af_bcftools(v, bcftools)))
        readers.append(("gzip", lambda v=vcf: _iter_eur_af_python(v)))

        rare = mono = found = 0
        ok = False
        for name, reader in readers:
            # Счётчики обнуляются перед КАЖДОЙ попыткой: если быстрый путь
            # отвалился на середине файла, часть записей уже посчитана, и
            # запасной путь посчитал бы их второй раз.
            rare = mono = found = 0
            try:
                for raw in reader():
                    found += 1
                    maf = _maf(raw)
                    if maf is None:
                        continue
                    if maf <= 0.0:
                        mono += 1
                    if maf < RARE_MAF:
                        rare += 1
                ok = True
                break
            except (PanelAdvisorError, OSError, EOFError) as e:
                logger.info("Чтение %s способом '%s' не удалось (%s)",
                            vcf.name, name, e)
                continue

        if not ok:
            logger.warning(
                "Не удалось прочитать %s ни одним способом — файл исключён "
                "из подсчёта состава чипа, числа ниже занижены.", vcf.name,
            )
            rare = mono = found = 0
        comp.rare += rare
        comp.monomorphic_eur += mono
        comp.per_chrom_matched[chrom] = found
        comp.matched += found

    if not comp.chip_positions:
        # Знаменатель не нашёлся — тогда "отсутствующих" честно 0, а не
        # выдуманное число.
        comp.chip_positions = comp.matched
    return comp


# ---------------------------------------------------------------------------
# Рекомендация
# ---------------------------------------------------------------------------
@dataclass
class Recommendation:
    panel: str                   # PANEL_HRC | PANEL_TOPMED
    reasons: list[str] = field(default_factory=list)
    counter_reasons: list[str] = field(default_factory=list)
    ancestry_hint: str = ""      # место под оценку происхождения
    composition: Optional[ChipComposition] = None

    @property
    def panel_title(self) -> str:
        return {PANEL_HRC: "HRC r1.1", PANEL_TOPMED: "TopMed r3"}.get(
            self.panel, self.panel
        )

    @property
    def headline(self) -> str:
        if self.reasons:
            return f"Рекомендую {self.panel_title}: " + "; ".join(self.reasons)
        return f"Рекомендую {self.panel_title}"


def recommend_panel(
    comp: ChipComposition,
    *,
    sample_is_european: Optional[bool] = None,
) -> Recommendation:
    """
    Читаемое правило: неевропейское происхождение ИЛИ тяжёлый редкий
    хвост -> TopMed; иначе HRC как более быстрый и без лифтовера.

    sample_is_european=None означает "происхождение не оценивалось" — это
    штатный режим: оценка ancestry в этой версии не реализована (см.
    докстринг модуля), и правило работает на одном сигнале из двух.
    """
    reasons: list[str] = []
    counter: list[str] = []
    ancestry_hint = ""

    heavy_tail = comp.rare_pct >= RARE_TAIL_PCT_THRESHOLD
    many_missing = comp.missing_pct >= MISSING_PCT_THRESHOLD

    if sample_is_european is False:
        ancestry_hint = "образец оценён как неевропейский"
        reasons.append(
            "происхождение образца не европейское — подвыборка HRC "
            "(почти сплошь европейцы) для него заведомо хуже"
        )
    elif sample_is_european is True:
        ancestry_hint = "образец оценён как европейский"
    else:
        ancestry_hint = ("происхождение не оценивалось — правило работает "
                         "по составу чипа")

    if heavy_tail:
        reasons.append(
            f"доля маркеров с MAF(EUR) ниже {RARE_MAF * 100:.1f}% — "
            f"{comp.rare_pct:.1f}% (порог {RARE_TAIL_PCT_THRESHOLD:.0f}%), "
            f"а редкий хвост HRC обрезан по MAC>=5"
        )
    # Доля отсутствующих в 1000G — только подпорка к уже принятому
    # решению (см. комментарий у MISSING_PCT_THRESHOLD): сама по себе
    # панель по ней не переключается.
    missing_note = ""
    if many_missing:
        missing_note = (
            f"{comp.missing_pct:.1f}% позиций чипа отсутствуют в 1000G "
            f"phase3 ({comp.missing:,} из {comp.chip_positions:,}) — "
            f"чип насыщен нестандартными маркерами"
        ).replace(",", " ")
        if reasons:
            reasons.append(missing_note)

    if reasons and sample_is_european is not False:
        counter.append(
            "цена TopMed — лифтовер GRCh37 -> GRCh38 и обратно: лишний шаг "
            "с лишним риском, и результат приходит в другой сборке"
        )
    panel = PANEL_TOPMED if reasons else PANEL_HRC
    if panel == PANEL_HRC:
        reasons.append(
            f"редкий хвост умеренный ({comp.rare_pct:.1f}% маркеров с "
            f"MAF(EUR) ниже {RARE_MAF * 100:.1f}%), позиций вне 1000G "
            f"{comp.missing_pct:.1f}% — выигрыш TopMed не окупает лифтовер"
        )
        if missing_note:
            counter.append(missing_note)
    return Recommendation(
        panel=panel, reasons=reasons, counter_reasons=counter,
        ancestry_hint=ancestry_hint, composition=comp,
    )


def format_recommendation(rec: Recommendation) -> str:
    comp = rec.composition
    lines = [
        "=" * 70,
        "СОСТАВ ЧИПА И ВЫБОР РЕФЕРЕНСНОЙ ПАНЕЛИ",
        "=" * 70,
    ]
    if comp is not None:
        lines += [
            f"Позиций чипа:            {comp.chip_positions:,}".replace(",", " "),
            f"Найдено в 1000G phase3:  {comp.matched:,}".replace(",", " "),
            f"Нет в 1000G вовсе:       {comp.missing:,} ({comp.missing_pct:.2f}%)".replace(",", " "),
            f"MAF(EUR) ниже {RARE_MAF * 100:.1f}%:      {comp.rare:,} "
            f"({comp.rare_pct:.2f}% от найденных)".replace(",", " "),
            f"Мономорфны в EUR:        {comp.monomorphic_eur:,} "
            f"({comp.monomorphic_eur_pct:.2f}%)".replace(",", " "),
            "",
        ]
    lines.append(rec.headline)
    for c in rec.counter_reasons:
        lines.append(f"  ⚠ {c}")
    lines.append(f"  ℹ {rec.ancestry_hint}")
    lines.append(
        "  ℹ Пороги предварительные и требуют калибровки: копите запуски в "
        "runs_metrics.csv (метрики исходника, панель, Rsq, финальный call "
        "rate, вердикт приёмки) — граница нарисуется сама."
    )
    lines.append("=" * 70)
    return "\n".join(lines)


def metrics_row(comp: Optional[ChipComposition],
                rec: Optional[Recommendation]) -> dict:
    """Плоский набор чисел для core/metrics_log.py."""
    if comp is None:
        return {}
    return {
        "chip_positions": comp.chip_positions,
        "chip_matched_1000g": comp.matched,
        "chip_missing_1000g_pct": round(comp.missing_pct, 4),
        "chip_rare_eur_pct": round(comp.rare_pct, 4),
        "chip_monomorphic_eur_pct": round(comp.monomorphic_eur_pct, 4),
        "panel_recommended": rec.panel if rec else "",
    }
