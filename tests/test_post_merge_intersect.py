"""
tests/test_post_merge_intersect.py

Регрессионный тест на промт "post_merge_intersect не просто диагностика,
а критическое условие для TopMed":

  1. main.py::_post_merge_intersect() принимает kgp_all_vcf (уже готовый
     BGZF+tabix-индексированный результат _concat_donors()) и использует
     его НАПРЯМУЮ как аргумент -R вместо хрупкого текстового
     common_pos.txt (фикс "Failed to read the regions").
  2. Если kgp_all_vcf не передан / не существует, либо путь через него
     не сработал — используется fallback на старый текстовый путь, с
     тем же результатом фильтрации.
  3. Если ОБА пути проваливаются, функция бросает исключение
     (PureCoreError), а НЕ возвращает "тихий" неотфильтрованный merged —
     раньше сбой этого шага считался диагностикой и вызывающий код
     (main()/gui/app.py) продолжал сборку молча; теперь он фатален.

Тесты 1 и 2 требуют реального bcftools/tabix/bgzip в PATH (строят
настоящие BGZF+tabix-индексированные VCF) — пропускаются, если бинарники
не найдены. Тест 3 не требует bcftools вовсе: subprocess.run замокан
целиком на всех путях.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from main import _post_merge_intersect, PureCoreError  # noqa: E402

BCFTOOLS = shutil.which("bcftools")
TABIX = shutil.which("tabix")
BGZIP = shutil.which("bgzip")
_HAS_HTSLIB = bool(BCFTOOLS and TABIX and BGZIP)

pytestmark = pytest.mark.skipif(
    False, reason="module import always allowed; individual tests skip as needed",
)

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


def _write_bgzipped_vcf(tmp_path: Path, name: str, records: list[str]) -> Path:
    """Пишет обычный VCF, затем сжимает его в настоящий BGZF через bgzip
    и строит tabix-индекс — то, что реально ожидает bcftools view -R."""
    plain = tmp_path / f"{name}.plain.vcf"
    with plain.open("w", encoding="utf-8", newline="\n") as f:
        f.write(_VCF_HEADER)
        for line in records:
            f.write(line + "\n")

    gz_path = tmp_path / f"{name}.vcf.gz"
    result = subprocess.run(
        [BGZIP, "-c", str(plain)], capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    gz_path.write_bytes(result.stdout)

    tabix_res = subprocess.run(
        [TABIX, "-p", "vcf", str(gz_path)], capture_output=True, text=True,
    )
    assert tabix_res.returncode == 0, tabix_res.stderr
    return gz_path


def _record(pos: int, gt: str = "0/1") -> str:
    return f"1\t{pos}\trs{pos}\tA\tG\t.\tPASS\t.\tGT\t{gt}"


@pytest.mark.skipif(not _HAS_HTSLIB, reason="требуется bcftools/tabix/bgzip в PATH")
def test_post_merge_intersect_via_kgp_all_vcf_filters_correctly(tmp_path):
    """merged_vcf содержит позиции 100/200/300; donor-подвыборка (и,
    соответственно, kgp_all_vcf) покрывает только 100/300 — после
    intersect должна остаться только пересекающаяся часть, removed > 0."""
    merged_vcf = _write_bgzipped_vcf(
        tmp_path, "merged", [_record(100), _record(200), _record(300)],
    )
    donor_vcf = _write_bgzipped_vcf(
        tmp_path, "kgp_all", [_record(100), _record(300)],
    )

    output_vcf = tmp_path / "checked.vcf.gz"
    result_path, before, after = _post_merge_intersect(
        merged_vcf, [donor_vcf], output_vcf,
        bcftools_path=BCFTOOLS, kgp_all_vcf=donor_vcf,
    )

    assert before == 3
    assert after == 2
    assert result_path == output_vcf

    count_res = subprocess.run(
        [BCFTOOLS, "view", "-H", str(output_vcf)], capture_output=True, text=True,
    )
    remaining_positions = {
        int(line.split("\t")[1]) for line in count_res.stdout.splitlines() if line.strip()
    }
    assert remaining_positions == {100, 300}


@pytest.mark.skipif(not _HAS_HTSLIB, reason="требуется bcftools/tabix/bgzip в PATH")
def test_post_merge_intersect_falls_back_to_text_regions_when_kgp_all_missing(tmp_path):
    """kgp_all_vcf передан, но не существует на диске — функция должна
    тихо откатиться на старый текстовый путь (common_pos.txt, собранный
    через bcftools query по donor_vcfs) и дать тот же корректный результат."""
    merged_vcf = _write_bgzipped_vcf(
        tmp_path, "merged", [_record(100), _record(200), _record(300)],
    )
    donor_vcf = _write_bgzipped_vcf(
        tmp_path, "donor", [_record(100), _record(300)],
    )
    missing_kgp_all = tmp_path / "does_not_exist.vcf.gz"

    output_vcf = tmp_path / "checked.vcf.gz"
    result_path, before, after = _post_merge_intersect(
        merged_vcf, [donor_vcf], output_vcf,
        bcftools_path=BCFTOOLS, kgp_all_vcf=missing_kgp_all,
    )

    assert before == 3
    assert after == 2

    common_pos = output_vcf.parent / "common_pos.txt"
    assert common_pos.exists(), "fallback должен был построить common_pos.txt"
    lines = set(common_pos.read_text(encoding="utf-8").strip().splitlines())
    assert lines == {"1\t100", "1\t300"}


def test_post_merge_intersect_raises_when_both_paths_fail(tmp_path):
    """Если и путь через kgp_all_vcf, и текстовый fallback проваливаются
    (bcftools возвращает ненулевой код на обоих), функция должна бросить
    PureCoreError, а не молча вернуть неотфильтрованный merged. Не требует
    реального bcftools — subprocess.run замокан целиком."""
    merged_vcf = tmp_path / "merged.vcf.gz"
    merged_vcf.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")  # валидный gzip magic, содержимое не важно
    donor_vcf = tmp_path / "donor.vcf.gz"
    donor_vcf.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")
    kgp_all = tmp_path / "kgp_all.vcf.gz"
    kgp_all.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")
    output_vcf = tmp_path / "checked.vcf.gz"

    call_count = {"n": 0}

    def fake_run(cmd, capture_output=True, text=True, check=False):
        call_count["n"] += 1
        if cmd[1] == "view" and "-H" in cmd:
            # _count_records(merged_vcf) — нужно, чтобы дойти до самого
            # intersect (before считается ДО intersect и не должен сам
            # по себе провоцировать раннее исключение).
            return subprocess.CompletedProcess(cmd, returncode=0, stdout="1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n", stderr="")
        if cmd[1] == "query":
            # Построение common_pos.txt для текстового fallback.
            return subprocess.CompletedProcess(cmd, returncode=0, stdout="1\t100\n", stderr="")
        if cmd[1] == "view" and "-R" in cmd:
            # И путь через kgp_all_vcf, И путь через common_pos.txt должны
            # провалиться — оба вызова -R возвращают ненулевой код.
            return subprocess.CompletedProcess(cmd, returncode=1, stdout="", stderr="bcftools: Failed to read the regions")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    with patch("main.subprocess.run", side_effect=fake_run):
        with pytest.raises(PureCoreError):
            _post_merge_intersect(
                merged_vcf, [donor_vcf], output_vcf,
                bcftools_path="bcftools", kgp_all_vcf=kgp_all,
            )


if __name__ == "__main__":
    raise SystemExit(
        "Этот файл рассчитан на запуск через pytest: "
        "python -m pytest tests/test_post_merge_intersect.py -v"
    )
