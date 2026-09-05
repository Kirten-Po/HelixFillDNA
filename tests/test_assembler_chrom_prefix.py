"""
Проверка точечного фикса template/assembler.py::load_imputed_genotypes():
regions-файл (_panel_{chrom}.txt), передаваемый в `bcftools query -R`,
теперь содержит ОБЕ формы имени хромосомы ("1" и "chr1") — вне зависимости
от того, в каком виде реальный VCF результата MIS называет свой CHROM.

bcftools/tabix в этом окружении недоступны — subprocess.run() внутри
load_imputed_genotypes() подменяется моком, который:
  1. на вызов индексации (tabix -p vcf ...) — не используется вовсе,
     поскольку .tbi-заглушка создаётся заранее (ветка "индекс уже есть"
     пропускается);
  2. на вызов bcftools query -R <regions_file> ... — перехватывает
     аргументы, СЧИТЫВАЕТ содержимое regions_file (он ещё не удалён на
     этот момент — unlink() происходит уже после subprocess.run()) и
     возвращает пустой CompletedProcess (сам факт совпадения строк VCF
     здесь не проверяется — assembler.py делегирует это bcftools).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

# Под pytest путь к корню проекта уже добавлен в sys.path через
# tests/conftest.py (см. соседний файл) — этот блок нужен ТОЛЬКО для
# прямого запуска `python tests/test_assembler_chrom_prefix.py` без
# pytest, где conftest.py не подхватывается автоматически. Условие
# "не в sys.path" делает вставку безвредной и при пуске под pytest —
# просто ничего не добавит повторно.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from template.assembler import load_imputed_genotypes  # noqa: E402


def _make_fake_imputed_dir(tmp_path: Path, chrom: str) -> Path:
    imputed_dir = tmp_path / "rerun_results"
    imputed_dir.mkdir()
    vcf_path = imputed_dir / f"chr{chrom}.dose.vcf.gz"
    # ⚠ Именно так, а не .with_suffix(".vcf.gz.tbi"): with_suffix заменяет
    # только ПОСЛЕДНЕЕ расширение, и второй вариант дал бы
    # "chr1.dose.vcf.vcf.gz.tbi". Ровно эта ошибка годами жила в
    # assembler.py и роняла повторную сборку — тест повторял её и потому
    # не ловил.
    tbi_path = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    # Содержимое не важно — subprocess.run замокан целиком, реальный
    # bcftools/tabix эти файлы не читает. Существование vcf_path нужно,
    # чтобы load_imputed_genotypes() не пропустил хромосому как
    # отсутствующую (см. `if not vcf_path.exists(): continue`).
    vcf_path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")  # валидный gzip magic
    tbi_path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00")
    return imputed_dir


def test_regions_file_contains_both_chrom_forms(tmp_path):
    chrom = "1"
    imputed_dir = _make_fake_imputed_dir(tmp_path, chrom)

    captured_regions_content = {}

    def fake_run(cmd, capture_output=True, text=True, check=False):
        if "-R" in cmd:
            regions_path = Path(cmd[cmd.index("-R") + 1])
            # Файл на этот момент ещё должен существовать — unlink()
            # в assembler.py происходит уже ПОСЛЕ этого вызова.
            assert regions_path.exists(), "regions-файл должен существовать во время вызова bcftools"
            captured_regions_content["text"] = regions_path.read_text()
            stdout = ""  # пустой результат — сами данные VCF здесь не важны
        else:
            stdout = ""
        return subprocess.CompletedProcess(cmd, returncode=0, stdout=stdout, stderr="")

    # panel_pos — КАНОНИЧЕСКИЕ имена хромосом (без "chr"), ровно так, как
    # их отдают extract_skeleton()/ChainLiftover.lift() в реальном коде.
    panel_pos = [(chrom, 12345), (chrom, 67890)]

    with patch("template.assembler.subprocess.run", side_effect=fake_run):
        genotypes = load_imputed_genotypes(
            imputed_dir, sample_name="genotek", panel_pos=panel_pos,
            rsq_threshold=0.30, bcftools_path="bcftools", tabix_path="tabix",
        )

    assert genotypes == {}  # stdout был пуст — генотипов и не ожидается

    content = captured_regions_content["text"]
    lines = set(content.strip().splitlines())

    # Обе формы для каждой позиции — это и есть сам фикс.
    assert "1\t12345" in lines, f"канонической формы нет в regions-файле: {lines}"
    assert "chr1\t12345" in lines, f"chr-формы нет в regions-файле (регрессия фикса): {lines}"
    assert "1\t67890" in lines
    assert "chr1\t67890" in lines
    assert len(lines) == 4, f"ожидалось ровно 4 строки (2 позиции x 2 формы), получено: {lines}"

    print("OK: regions-файл содержит обе формы имени хромосомы")


def test_hrc_style_positions_are_idempotent_no_extra_matches_lost(tmp_path):
    """
    Для HRC (без лифтовера) поведение не должно меняться по существу:
    в regions-файле по-прежнему присутствует каноническая форма (та,
    которая реально совпадёт с CHROM в HRC-VCF без префикса "chr") —
    добавление "chr1"-строки лишь дополняет файл безвредной лишней
    строкой, не заменяет и не теряет исходную.
    """
    chrom = "7"
    imputed_dir = _make_fake_imputed_dir(tmp_path, chrom)
    captured = {}

    def fake_run(cmd, capture_output=True, text=True, check=False):
        if "-R" in cmd:
            regions_path = Path(cmd[cmd.index("-R") + 1])
            captured["text"] = regions_path.read_text()
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    panel_pos = [(chrom, 555)]
    with patch("template.assembler.subprocess.run", side_effect=fake_run):
        load_imputed_genotypes(
            imputed_dir, sample_name="genotek", panel_pos=panel_pos,
            bcftools_path="bcftools", tabix_path="tabix",
        )

    lines = set(captured["text"].strip().splitlines())
    assert "7\t555" in lines, "каноническая (HRC-совместимая) форма обязана остаться"
    assert "chr7\t555" in lines

    print("OK: каноническая форма (нужная для HRC) сохранена")


if __name__ == "__main__":
    import tempfile

    for fn in (test_regions_file_contains_both_chrom_forms,
               test_hrc_style_positions_are_idempotent_no_extra_matches_lost):
        with tempfile.TemporaryDirectory() as td:
            fn(Path(td))
    print("\nВСЕ ТЕСТЫ ПРОШЛИ")
