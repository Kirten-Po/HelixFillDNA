"""
tests/test_download_watch.py

Панель «что сейчас качается» дважды подряд показала пользователю строку,
означавшую не то, что она читалась. Оба случая закреплены тестами.

  1. Мигание. Файл показывался активным ТОЛЬКО в тот опрос, когда его
     размер вырос. На 11 КБ/с при опросе раз в 1,5 с это ~16 КБ, которые
     Windows часто держит в буфере, — строка то появлялась, то
     исчезала, и это читалось как «закачка прыгает между хромосомами».
  2. «0 Б/с — данных нет 20 с» у ДОСКАЧАННОГО файла. Починка первого
     бага заводила счётчик скорости при первой же встрече с файлом, со
     значением 0. А готовый архив появляется под своим именем именно как
     новый путь: во время закачки он зовётся <имя>.part и
     переименовывается только после проверки целостности.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from core.download_watch import (
    ACTIVE_GRACE_SECONDS, STALLED_SECONDS, DownloadWatcher,
)

ZIP = Path("chr_15.zip")
PART = Path("chr_15.zip.part")


def test_first_sighting_is_not_a_download():
    """
    Главная регрессия: файл, впервые попавший в поле зрения, активным не
    считается. Ни доскачанный архив после переименования, ни файл,
    оставшийся от прошлого запуска.
    """
    w = DownloadWatcher()
    assert w.poll([(ZIP, 67_000_000)], now=0.0) == []
    # и на следующем опросе тоже — размер не растёт, значит не качается
    assert w.poll([(ZIP, 67_000_000)], now=1.5) == []
    assert w.poll([(ZIP, 67_000_000)], now=30.0) == []


def test_rename_from_part_does_not_look_like_a_stalled_download():
    """
    Сквозной сценарий из жалобы: .part растёт, потом исчезает, вместо
    него появляется .zip того же размера. Готовый архив НЕ должен
    показаться как «0 Б/с — данных нет».
    """
    w = DownloadWatcher()
    w.poll([(PART, 10_000_000)], now=0.0)
    active = w.poll([(PART, 20_000_000)], now=1.0)
    assert len(active) == 1 and active[0].speed == pytest.approx(10_000_000)

    # закачка завершилась: .part переименован в .zip
    assert w.poll([(ZIP, 20_000_000)], now=2.0) == []
    assert w.poll([(ZIP, 20_000_000)], now=5.0) == []


def test_growth_makes_it_active_and_speed_is_measured():
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    active = w.poll([(PART, 3_000_000)], now=2.0)
    assert len(active) == 1
    assert active[0].speed == pytest.approx(1_500_000)
    assert active[0].size == 3_000_000
    assert active[0].idle == pytest.approx(0.0)


def test_stays_listed_between_polls_without_growth():
    """Против мигания: размер обновляется рывками, строка держится."""
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    w.poll([(PART, 100_000)], now=1.5)
    for t in (3.0, 4.5, 6.0, 10.0, 20.0):
        active = w.poll([(PART, 100_000)], now=t)
        assert len(active) == 1, f"строка пропала на {t} с"


def test_drops_out_after_grace_period():
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    w.poll([(PART, 100_000)], now=1.0)
    assert w.poll([(PART, 100_000)], now=1.0 + ACTIVE_GRACE_SECONDS - 0.1)
    assert w.poll([(PART, 100_000)], now=1.0 + ACTIVE_GRACE_SECONDS + 0.1) == []


def test_stalled_flag():
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    w.poll([(PART, 100_000)], now=1.0)
    fresh = w.poll([(PART, 100_000)], now=1.0 + STALLED_SECONDS - 0.1)[0]
    assert not fresh.stalled
    stale = w.poll([(PART, 100_000)], now=1.0 + STALLED_SECONDS + 0.1)[0]
    assert stale.stalled
    assert stale.idle == pytest.approx(STALLED_SECONDS + 0.1)


def test_restart_from_zero_is_tracked_correctly():
    """
    Повтор без докачки: .part удаляется и начинается с нуля. Со старой
    записью (46 МБ) новый файл не считался бы растущим, пока снова не
    дорос бы до 46 МБ — минуты «мёртвого» прогресса.
    """
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    w.poll([(PART, 46_000_000)], now=10.0)
    # файл исчез (удалён перед повторной попыткой)
    assert w.poll([], now=11.0) == []
    # начали заново с нуля
    assert w.poll([(PART, 0)], now=12.0) == []
    active = w.poll([(PART, 1_000_000)], now=13.0)
    assert len(active) == 1 and active[0].size == 1_000_000


def test_display_name_hides_part_suffix():
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    item = w.poll([(PART, 1)], now=1.0)[0]
    assert item.display_name == "chr_15.zip"
    assert item.path.name == "chr_15.zip.part"


def test_shrinking_file_is_not_negative_speed():
    """Файл может укоротиться (перезапись с нуля) — скорость не считаем."""
    w = DownloadWatcher()
    w.poll([(PART, 5_000_000)], now=0.0)
    assert w.poll([(PART, 1_000)], now=1.0) == []


def test_sorted_by_size_descending():
    w = DownloadWatcher()
    a, b = Path("chr_1.zip.part"), Path("chr_2.zip.part")
    w.poll([(a, 0), (b, 0)], now=0.0)
    active = w.poll([(a, 1_000), (b, 9_000)], now=1.0)
    assert [x.path for x in active] == [b, a]


def test_reset_forgets_everything():
    w = DownloadWatcher()
    w.poll([(PART, 0)], now=0.0)
    w.poll([(PART, 100)], now=1.0)
    w.reset()
    assert w.poll([(PART, 100)], now=2.0) == []
