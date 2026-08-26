"""
adapters/base.py
Единый контракт для всех адаптеров источников (FTDNA, MyHeritage и т.д.).
Все адаптеры возвращают ParseResult — дальше с ним работает common_core.py,
не зная, из какого формата пришли данные.
"""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass(frozen=True)
class ParsedVariant:
    """Один нормализованный вариант после разрешения ориентации."""
    rsid: str
    chrom: str          # уже нормализован: "1".."22", "X", "Y", "MT"
    pos: int            # 1-based
    ref: str            # одна буква ACGT
    alt: str            # одна буква ACGT или "." для 0/0
    gt: str             # "0/0" | "0/1" | "1/1"


@dataclass
class ParseResult:
    """Результат работы адаптера: список вариантов + счётчики QC."""
    variants: list[ParsedVariant] = field(default_factory=list)
    total_measured: int = 0
    missing: int = 0
    het_self_complementary: int = 0
    both_non_ref: int = 0
    invalid_codes: int = 0
    malformed_rows: int = 0
    ref_non_acgt: int = 0
    duplicate_positions: int = 0  # позиция (chrom,pos) уже встречалась в файле
    # Промт "HRC / TopMed", лифтовер координат: сколько позиций НЕ удалось
    # перенести из GRCh37 в целевую сборку панели (нет chain-блока для
    # хромосомы, позиция попала в разрыв между блоками выравнивания, либо
    # результат оказался бы вне границ целевой хромосомы — см.
    # core/liftover.py::ChainLiftover.lift()). Такие позиции полностью
    # исключаются из result.variants и из широкой сигнатуры чипа (Задача D) —
    # они физически не существуют в целевой сборке. Остаётся 0 для panel="hrc"
    # (liftover=None, лифтовер не применяется, поведение не меняется) и для
    # вызовов без параметра liftover вовсе (обратная совместимость).
    lift_failed: int = 0
    chip_signature: str = ""

    @property
    def het_aa(self) -> int:
        return (
            self.het_self_complementary
            + self.both_non_ref
            + self.invalid_codes
        )

    @property
    def both_non_ref_pct(self) -> float:
        if self.total_measured == 0:
            return 0.0
        return 100.0 * self.both_non_ref / self.total_measured

    @property
    def het_self_complementary_pct(self) -> float:
        if self.total_measured == 0:
            return 0.0
        return 100.0 * self.het_self_complementary / self.total_measured
