"""
tools/check_output.py — проверка качества собранного файла 23andMe.

bcftools здесь бесполезен: итоговый файл — текстовый TSV, а не VCF.
Скрипт считает три независимые вещи (см. README раздела «Что нового»):

  1. Call rate — доля непустых генотипов, всего и по хромосомам.
     Именно по этой цифре Генотек отклоняет файлы, и именно на X она
     проседала до 69,6 % до появления импутации X.

  2. Сверку с исходным файлом чипа на ОБЩИХ позициях. Прямые измерения
     имеют приоритет над импутацией, поэтому здесь ожидается ~100 %
     совпадений; любое расхождение — признак ошибки конвертации
     (ориентация аллелей, лифтовер, сдвиг позиций), а не «неточной
     импутации».

  3. Если передан третий файл — РЕАЛЬНЫЙ экспорт того же человека с
     другого чипа (например, Генотека) — считается точность именно
     ИМПУТИРОВАННЫХ позиций: тех, которых на исходном чипе не было.
     Это единственная честная оценка качества импутации, потому что
     сверять импутацию с тем же файлом, из которого она построена,
     бессмысленно.

Использование:
    python tools/check_output.py results/out.txt 37_S_Polomoshnov...csv
    python tools/check_output.py results/out.txt источник.csv genotek.txt
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

NO_CALL = {"--", "-", "00", "", "__"}
CHROM_ORDER = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "MT": 25}


def read_23andme(path: Path) -> dict[tuple[str, str], str]:
    """{(хромосома, позиция): генотип} из файла в оформлении 23andMe."""
    out = {}
    with path.open("r", encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 4:
                continue
            out[(p[1], p[2])] = p[3].upper()
    return out


def read_ftdna(path: Path) -> dict[tuple[str, str], str]:
    """То же из CSV FamilyTreeDNA. XY (псевдоаутосомный регион) сводится
    к X — так же, как это делает адаптер (adapters/ftdna_v3.py)."""
    import csv

    out = {}
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) < 4:
                continue
            chrom = {"XY": "X", "23": "X", "24": "Y", "25": "X", "26": "MT"}.get(
                row[1], row[1])
            out[(chrom, row[2])] = row[3].upper()
    return out


def read_any(path: Path) -> dict[tuple[str, str], str]:
    return read_ftdna(path) if path.suffix.lower() == ".csv" else read_23andme(path)


def called(gt: str) -> bool:
    return gt not in NO_CALL and "-" not in gt


def same(a: str, b: str) -> bool:
    """Генотипы без учёта порядка аллелей: AG == GA. Гаплоидная запись
    ("A") сравнивается с гомозиготной ("AA") как равная — на X у мужчин
    оба варианта означают одно и то же."""
    a = a * 2 if len(a) == 1 else a
    b = b * 2 if len(b) == 1 else b
    return sorted(a) == sorted(b)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    out_path = Path(sys.argv[1])
    result = read_23andme(out_path)

    # --- 1. Call rate ----------------------------------------------------
    per = defaultdict(lambda: [0, 0])   # хромосома -> [всего, заполнено]
    for (chrom, _pos), gt in result.items():
        per[chrom][0] += 1
        per[chrom][1] += called(gt)

    total = sum(v[0] for v in per.values())
    filled = sum(v[1] for v in per.values())
    print(f"\n=== {out_path.name} ===")
    print(f"Позиций: {total}, заполнено: {filled} ({100 * filled / total:.2f} %)\n")
    print(f"{'хр.':>4} {'позиций':>9} {'заполнено':>10} {'call rate':>10}")
    for chrom in sorted(per, key=lambda c: CHROM_ORDER.get(c, 99)):
        n, ok = per[chrom]
        print(f"{chrom:>4} {n:>9} {ok:>10} {100 * ok / n:>9.2f} %")

    if len(sys.argv) < 3:
        return 0

    # --- 2. Сверка с исходным чипом --------------------------------------
    src = read_any(Path(sys.argv[2]))
    shared = match = 0
    mismatches = []
    for key, gt_src in src.items():
        gt_out = result.get(key)
        if gt_out is None or not called(gt_src) or not called(gt_out):
            continue
        shared += 1
        if same(gt_src, gt_out):
            match += 1
        elif len(mismatches) < 10:
            mismatches.append((key, gt_src, gt_out))

    print(f"\n=== Сверка с исходным чипом ({Path(sys.argv[2]).name}) ===")
    print(f"Общих измеренных позиций: {shared}")
    if shared:
        print(f"Совпадает: {match} ({100 * match / shared:.3f} %)")
        print("Ожидается ~100 %: прямые измерения перекрывают импутацию.")
        print("Расхождения здесь — ошибка конвертации, а не качество импутации.")
    for key, a, b in mismatches:
        print(f"  {key[0]}:{key[1]}  чип={a}  файл={b}")

    if len(sys.argv) < 4:
        return 0

    # --- 3. Точность импутации по независимому файлу ---------------------
    truth = read_any(Path(sys.argv[3]))
    imp_n = imp_ok = 0
    per_chrom = defaultdict(lambda: [0, 0])
    for key, gt_truth in truth.items():
        if key in src:            # позиция была на чипе — не импутирована
            continue
        gt_out = result.get(key)
        if gt_out is None or not called(gt_truth) or not called(gt_out):
            continue
        imp_n += 1
        ok = same(gt_truth, gt_out)
        imp_ok += ok
        per_chrom[key[0]][0] += 1
        per_chrom[key[0]][1] += ok

    print(f"\n=== Точность импутации против {Path(sys.argv[3]).name} ===")
    print("Только позиции, которых НЕ было на исходном чипе.")
    if not imp_n:
        print("Общих импутированных позиций не нашлось.")
        return 0
    print(f"Сравнено: {imp_n}, совпало: {imp_ok} ({100 * imp_ok / imp_n:.2f} %)\n")
    print(f"{'хр.':>4} {'сравнено':>9} {'совпало':>9} {'точность':>10}")
    for chrom in sorted(per_chrom, key=lambda c: CHROM_ORDER.get(c, 99)):
        n, ok = per_chrom[chrom]
        print(f"{chrom:>4} {n:>9} {ok:>9} {100 * ok / n:>9.2f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
