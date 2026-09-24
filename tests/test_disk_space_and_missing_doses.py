"""
Нехватка места при распаковке и пропавшая хромосома на Этапе 7.

Живой случай: 7z упал на chr_4.zip с «Cannot set length for output file»
(на диске кончилось место), а текст ошибки пришёл кракозябрами. Если
после такого файл chr4.dose.vcf.gz отсутствует, сборка раньше молча
оставляла всю хромосому пустой.
"""
from __future__ import annotations

import zipfile
from collections import namedtuple
from pathlib import Path

import pytest

from core import archive_utils
from core.archive_utils import ArchiveExtractionError
from template.assembler import (
    AssemblyError, load_imputed_genotypes, missing_autosome_doses,
)

_Usage = namedtuple("_Usage", "total used free")


def _make_zip(path: Path, name: str, size: int) -> Path:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(name, b"\0" * size)
    return path


def test_required_space_counts_uncompressed_minus_existing(tmp_path):
    z = _make_zip(tmp_path / "chr_4.zip", "chr4.dose.vcf.gz", 10_000)
    dest = tmp_path / "out"
    dest.mkdir()
    assert archive_utils.required_extract_space([z], dest) == 10_000
    (dest / "chr4.dose.vcf.gz").write_bytes(b"\0" * 4_000)
    assert archive_utils.required_extract_space([z], dest) == 6_000


def test_check_free_space_refuses_before_extracting(tmp_path, monkeypatch):
    z = _make_zip(tmp_path / "chr_4.zip", "chr4.dose.vcf.gz", 10_000)
    dest = tmp_path / "out"
    monkeypatch.setattr(archive_utils.shutil, "disk_usage",
                        lambda _p: _Usage(0, 0, 1024))
    with pytest.raises(ArchiveExtractionError, match="Недостаточно места"):
        archive_utils.extract_all([z], dest, "pwd")
    assert not (dest / "chr4.dose.vcf.gz").exists()


def test_check_free_space_passes_when_enough(tmp_path, monkeypatch):
    z = _make_zip(tmp_path / "chr_4.zip", "chr4.dose.vcf.gz", 10_000)
    monkeypatch.setattr(archive_utils.shutil, "disk_usage",
                        lambda _p: _Usage(0, 0, 10 * 1024 ** 3))
    archive_utils.check_free_space([z], tmp_path / "out")


@pytest.mark.parametrize("text", [
    "ERROR: Cannot set length for output file : Недостаточно места на диске.",
    "ERROR: There is not enough space on the disk.",
    "No space left on device",
])
def test_disk_full_is_recognised(text):
    assert archive_utils._is_disk_full(text)


def test_wrong_password_is_not_disk_full():
    assert not archive_utils._is_disk_full("ERROR: Wrong password : chr4.dose.vcf.gz")


def _touch_doses(d: Path, chroms):
    for c in chroms:
        (d / f"chr{c}.dose.vcf.gz").write_bytes(b"x")


def test_missing_autosome_doses_reports_gap(tmp_path):
    _touch_doses(tmp_path, [str(i) for i in range(1, 23) if i != 4])
    assert missing_autosome_doses(tmp_path) == ["4"]


def test_empty_dose_file_counts_as_missing(tmp_path):
    _touch_doses(tmp_path, [str(i) for i in range(1, 23)])
    (tmp_path / "chr4.dose.vcf.gz").write_bytes(b"")
    assert missing_autosome_doses(tmp_path) == ["4"]


def test_x_is_not_required(tmp_path):
    _touch_doses(tmp_path, [str(i) for i in range(1, 23)])
    assert missing_autosome_doses(tmp_path) == []


def test_load_refuses_to_assemble_without_a_chromosome(tmp_path):
    _touch_doses(tmp_path, [str(i) for i in range(1, 23) if i != 4])
    with pytest.raises(AssemblyError, match="chr4"):
        load_imputed_genotypes(tmp_path, require_autosomes=True)
