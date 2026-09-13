"""
Отсечка импутированных вызовов по max(GP) вместо Rsq.

Почему метрика меняется, а не только порог: Rsq — качество ПОЗИЦИИ по
всей выборке (метрика из GWAS), GP — уверенность вызова У ЭТОГО человека
В ЭТОЙ позиции. На ультраредких вариантах (MAF < 0,1 %) разброса
дозировок нет, Rsq уходит в ноль по статистической причине — и позиция
отбрасывается, хотя модель ставит там гомозиготу по референсу с
max(GP) >= 0,99. Размер выигрыша зависит от числа доноров (при 20 Rsq
почти не работает, при 503 работает) и от покрытия трафарета чипом:
замеры дают от +1 до +15 п.п., на полном наборе доноров — около нуля.

bcftools/tabix в этом окружении могут отсутствовать — subprocess.run()
подменяется моком, который отвечает на два вызова:
  * `bcftools view -h` — шапка VCF (по ней _has_format_tag() решает,
    можно ли вообще запрашивать GP: на неизвестном теге настоящий
    bcftools не возвращает точки, а падает с ошибкой);
  * `bcftools query -f ...` — строки с генотипами.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from template.assembler import (  # noqa: E402
    QUALITY_GP, QUALITY_RSQ, _max_gp, load_imputed_genotypes,
)

HEADER_WITH_GP = (
    '##fileformat=VCFv4.2\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posterior">\n'
)
HEADER_WITHOUT_GP = (
    '##fileformat=VCFv4.2\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)

#: Три позиции chr22. Первая — ультраредкая: Rsq почти ноль, но вызов
#: уверенный (ровно тот случай, ради которого метрика и меняется).
#: Вторая — обратная: Rsq приличный, а вызов размазан между генотипами.
#: Третья проходит по обеим метрикам.
ROWS_WITH_GP = (
    "22\t100\tA\tG\t0|0\t0.999,0.001,0\n"
    "22\t200\tC\tT\t0|1\t0.30,0.45,0.25\n"
    "22\t300\tG\tA\t1|1\t0.001,0.009,0.99\n"
)
ROWS_NO_GP = (
    "22\t100\tA\tG\t0|0\n"
    "22\t200\tC\tT\t0|1\n"
    "22\t300\tG\tA\t1|1\n"
)
INFO_LINES = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "22\t100\t.\tA\tG\t.\t.\tIMPUTED;MAF=0.0004;R2=0.04\n"
    "22\t200\t.\tC\tT\t.\t.\tIMPUTED;MAF=0.21;R2=0.61\n"
    "22\t300\t.\tG\tA\t.\t.\tIMPUTED;MAF=0.33;R2=0.97\n"
)


def _make_dir(tmp_path: Path, with_gp: bool = True) -> Path:
    d = tmp_path / "rerun_results"
    d.mkdir()
    vcf = d / "chr22.dose.vcf.gz"
    vcf.write_bytes(b"")
    (vcf.with_suffix(vcf.suffix + ".tbi")).write_bytes(b"")
    import gzip
    with gzip.open(d / "chr22.info.gz", "wt", encoding="utf-8") as f:
        f.write(INFO_LINES)
    return d


def _fake_run(with_gp: bool):
    def run(cmd, *a, **kw):
        if "view" in cmd and "-h" in cmd:
            out = HEADER_WITH_GP if with_gp else HEADER_WITHOUT_GP
        else:
            fmt = cmd[cmd.index("-f") + 1] if "-f" in cmd else ""
            out = ROWS_WITH_GP if "[%GP]" in fmt else ROWS_NO_GP
        return subprocess.CompletedProcess(cmd, 0, stdout=out, stderr="")
    return run


def test_gp_keeps_ultra_rare_and_drops_uncertain(tmp_path):
    """GP возвращает уверенный вызов с нулевым Rsq и убирает размазанный."""
    d = _make_dir(tmp_path)
    with patch("template.assembler.subprocess.run", side_effect=_fake_run(True)):
        got = load_imputed_genotypes(d, quality=QUALITY_GP, gp_threshold=0.90)
    assert got == {"22_100": "AA", "22_300": "AA"}


def test_rsq_does_the_opposite(tmp_path):
    """Прежняя метрика на тех же данных делает ровно обратный выбор."""
    d = _make_dir(tmp_path)
    with patch("template.assembler.subprocess.run", side_effect=_fake_run(True)):
        got = load_imputed_genotypes(d, quality=QUALITY_RSQ, rsq_threshold=0.30)
    assert got == {"22_200": "CT", "22_300": "AA"}


def test_falls_back_to_rsq_when_no_gp_in_header(tmp_path):
    """Без FORMAT/GP в шапке GP не запрашивается — иначе bcftools упал бы."""
    d = _make_dir(tmp_path)
    with patch("template.assembler.subprocess.run", side_effect=_fake_run(False)):
        got = load_imputed_genotypes(d, quality=QUALITY_GP, gp_threshold=0.90,
                                     rsq_threshold=0.30)
    assert got == {"22_200": "CT", "22_300": "AA"}


def test_quality_out_carries_the_metric_actually_used(tmp_path):
    """
    В rsq_out кладётся значение ТОЙ метрики, которой отсекали: по нему
    core/rsq_tuner.py строит кривую, и смешивать шкалы нельзя.
    """
    d = _make_dir(tmp_path)
    out: dict[str, float] = {}
    with patch("template.assembler.subprocess.run", side_effect=_fake_run(True)):
        load_imputed_genotypes(d, quality=QUALITY_GP, gp_threshold=0.90,
                               rsq_out=out)
    assert out == {"22_100": 0.999, "22_300": 0.99}


def test_max_gp_parsing():
    assert _max_gp("0.999,0.001,0") == 0.999
    assert _max_gp("0.2,0.8") == 0.8          # гаплоидный вызов: две компоненты
    assert _max_gp(".") is None               # значения нет -> откат на Rsq
    assert _max_gp("") is None
    assert _max_gp("не число") is None


def test_unreadable_header_is_not_reported_as_missing_gp(tmp_path, caplog):
    """
    Если сам вызов bcftools не состоялся, это НЕ значит "в дозах нет GP".
    Сообщение должно винить вызов, а не файл: иначе диагностика уводит в
    сторону — ровно так и случилось, когда в функцию передали "bcftools"
    вместо настроенного пути, и предупреждение обвинило выгрузку TopMed,
    в которой поле GP на самом деле есть.
    """
    import logging

    d = _make_dir(tmp_path)

    def run(cmd, *a, **kw):
        if "view" in cmd and "-h" in cmd:
            raise FileNotFoundError(2, "Не удается найти указанный файл")
        return subprocess.CompletedProcess(cmd, 0, stdout=ROWS_NO_GP, stderr="")

    with caplog.at_level(logging.WARNING):
        with patch("template.assembler.subprocess.run", side_effect=run):
            got = load_imputed_genotypes(d, quality=QUALITY_GP, gp_threshold=0.90,
                                         rsq_threshold=0.30)

    # Откат на Rsq состоялся — работа не потеряна.
    assert got == {"22_200": "CT", "22_300": "AA"}

    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "Не удалось прочитать шапку" in text
    assert "в дозах нет поля FORMAT/GP" not in text


def test_has_format_tag_returns_none_when_call_fails(tmp_path):
    """Три исхода: True (есть), False (нет), None (проверить не удалось)."""
    from template.assembler import _has_format_tag

    vcf = tmp_path / "chr22.dose.vcf.gz"
    vcf.write_bytes(b"")

    def ok(cmd, *a, **kw):
        return subprocess.CompletedProcess(cmd, 0, stdout=HEADER_WITH_GP, stderr="")

    def missing(cmd, *a, **kw):
        return subprocess.CompletedProcess(cmd, 0, stdout=HEADER_WITHOUT_GP, stderr="")

    def broken(cmd, *a, **kw):
        raise OSError("bcftools не найден")

    with patch("template.assembler.subprocess.run", side_effect=ok):
        assert _has_format_tag(vcf, "GP", "bcftools") is True
    with patch("template.assembler.subprocess.run", side_effect=missing):
        assert _has_format_tag(vcf, "GP", "bcftools") is False
    with patch("template.assembler.subprocess.run", side_effect=broken):
        assert _has_format_tag(vcf, "GP", "bcftools") is None
