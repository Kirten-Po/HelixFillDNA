"""
core/download_watch.py
Учёт «что сейчас качается» по размерам файлов на диске.

Почему это отдельный модуль, а не десяток строк внутри gui/app.py:
логика жила прямо в методе Tk-класса, тестировать её было нечем — и она
дважды подряд выдала пользователю строку, которая означала не то, что
читалась:

  1. Файл показывался активным ТОЛЬКО в тот опрос, когда его размер
     вырос. На скорости 11 КБ/с при опросе раз в 1,5 с это ~16 КБ,
     которые Windows часто ещё держит в буфере, — строка мигала, и это
     читалось как «закачка прыгает между хромосомами», хотя архивы
     качаются строго по одному.
  2. Починка через «держать строку ещё N секунд после последнего роста»
     завела счётчик скорости при ПЕРВОЙ же встрече с файлом, со
     значением 0. А готовый архив появляется под своим именем именно как
     новый путь: пока он качается, он зовётся <имя>.part и
     переименовывается только после проверки целостности. В итоге
     ДОСКАЧАННЫЙ файл 20 секунд висел как «0 Б/с — данных нет».

Отсюда правило, которое модуль соблюдает буквально: **файл считается
активным только после того, как рост его размера увидели своими
глазами.** Ни первая встреча, ни просто присутствие файла на диске
активностью не считаются.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

#: Сколько секунд после последнего роста файл ещё показывается активным.
#: Размер файла на Windows обновляется рывками (буферизация), и без этого
#: запаса строка мигает на медленной закачке.
ACTIVE_GRACE_SECONDS = 25.0

#: После этого простоя к строке дописывается «данных нет N с» — честнее,
#: чем молча показывать последнюю измеренную скорость.
STALLED_SECONDS = 5.0

#: Рабочее расширение файла во время закачки (mis_adapter.download_one
#: пишет в него и переименовывает только после проверки целостности).
PART_SUFFIX = ".part"


@dataclass(frozen=True)
class ActiveFile:
    path: Path
    size: int
    speed: float      # байт/с, по последнему замеренному приросту
    idle: float       # секунд с последнего роста

    @property
    def display_name(self) -> str:
        """Имя без служебного .part — пользователю он ничего не говорит."""
        name = self.path.name
        return name[: -len(PART_SUFFIX)] if name.endswith(PART_SUFFIX) else name

    @property
    def stalled(self) -> bool:
        return self.idle > STALLED_SECONDS


class DownloadWatcher:
    """
    Хранит историю размеров и отвечает на вопрос «что качается прямо
    сейчас». Ничего не знает ни про Tk, ни про файловую систему: на вход
    ему дают уже собранный список (путь, размер).
    """

    def __init__(self, grace_seconds: float = ACTIVE_GRACE_SECONDS):
        self.grace_seconds = grace_seconds
        self._sizes: dict[Path, tuple[float, int]] = {}   # путь -> (когда вырос, размер)
        self._speeds: dict[Path, tuple[float, float]] = {}  # путь -> (когда вырос, скорость)

    def reset(self) -> None:
        self._sizes.clear()
        self._speeds.clear()

    def poll(self, entries: Iterable[tuple[Path, int]], now: float) -> list[ActiveFile]:
        """
        Возвращает активные закачки, от большей к меньшей.

        entries — все файлы, которые сейчас видно на диске. Пути, которых
        в этом списке нет, забываются: готовый архив переименовывается из
        <имя>.part в <имя>.zip, и запись про .part иначе висела бы вечно;
        а при повторной попытке без докачки .part удаляется и начинается
        с нуля — со старой записью (скажем, 46 МБ) новый файл не считался
        бы растущим, пока снова не дорос бы до 46 МБ, то есть минуты
        «мёртвого» прогресса.
        """
        seen: set[Path] = set()
        active: list[ActiveFile] = []

        for path, size in entries:
            seen.add(path)
            previous = self._sizes.get(path)
            if previous is None:
                # Первая встреча — это ещё не закачка (см. докстринг
                # модуля, пункт 2). Просто запоминаем отправную точку.
                self._sizes[path] = (now, size)
                continue

            prev_time, prev_size = previous
            if size > prev_size:
                speed = (size - prev_size) / max(0.001, now - prev_time)
                self._sizes[path] = (now, size)
                self._speeds[path] = (now, speed)

            measured = self._speeds.get(path)
            if measured is None:
                continue
            last_growth, speed = measured
            idle = now - last_growth
            if idle > self.grace_seconds:
                continue
            active.append(ActiveFile(path=path, size=size, speed=speed, idle=idle))

        for forgotten in [p for p in self._sizes if p not in seen]:
            self._sizes.pop(forgotten, None)
            self._speeds.pop(forgotten, None)

        active.sort(key=lambda a: -a.size)
        return active
