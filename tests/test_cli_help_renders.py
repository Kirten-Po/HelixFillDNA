"""
`main.py --help` должен строиться и печататься на windows-консоли.

Два независимых способа уронить справку, оба реально случались:

  1. **Одиночный `%`.** argparse прогоняет каждую help-строку через
     оператор `%` (ArgumentDefaultsHelpFormatter._expand_help), поэтому
     "98,5 % их чипа" превращается в попытку форматирования по
     спецификатору "% и" -> ValueError: unsupported format character.
     Так справка была сломана с 1.3.0, когда в help --format попали доли
     покрытия чипа.

  2. **Символ вне cp1251.** argparse печатает справку одним куском, и на
     русской windows-консоли (кодовая страница 1251) первый же символ
     вроде "⚠" (U+26A0) роняет вывод через UnicodeEncodeError. В логах и
     print() такие символы безвредны — там поток уже в UTF-8, — но в
     справке недопустимы.

Тест проверяет обе поломки сразу и на всех аргументах, а не только на
тех, что правились последними.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))


def _build_parser() -> argparse.ArgumentParser:
    """
    Достаёт готовый парсер из main._parse_args(), не вызывая сам разбор:
    parse_args() подменяется перехватчиком, который запоминает self и
    немедленно выходит.
    """
    main = importlib.import_module("main")
    captured: dict[str, argparse.ArgumentParser] = {}
    original = argparse.ArgumentParser.parse_args

    def fake(self, *args, **kwargs):
        captured["parser"] = self
        raise SystemExit(0)

    argparse.ArgumentParser.parse_args = fake
    saved_argv = sys.argv
    try:
        sys.argv = ["main.py"]
        try:
            main._parse_args()
        except SystemExit:
            pass
    finally:
        argparse.ArgumentParser.parse_args = original
        sys.argv = saved_argv
    assert "parser" in captured, "не удалось перехватить парсер main._parse_args()"
    return captured["parser"]


def test_help_builds_without_percent_error():
    """format_help() не должен падать на одиночных знаках процента."""
    parser = _build_parser()
    text = parser.format_help()          # здесь и происходил ValueError
    assert "--quality" in text
    assert "--gp-threshold" in text


def test_help_is_printable_on_cp1251_console():
    """Вся справка должна кодироваться в cp1251 — иначе `--help` упадёт."""
    text = _build_parser().format_help()
    try:
        text.encode("cp1251")
    except UnicodeEncodeError as e:
        bad = text[e.start:e.end]
        # Показываем сам символ и его окружение: без этого искать
        # виновный аргумент в сотне строк справки неудобно.
        context = text[max(0, e.start - 70):e.end + 70].replace("\n", " ")
        pytest.fail(
            f"в справке есть символ {bad!r} (U+{ord(bad):04X}), которого нет "
            f"в cp1251 — на русской windows-консоли `main.py --help` упадёт "
            f"с UnicodeEncodeError.\nОкружение: ...{context}..."
        )


@pytest.mark.parametrize("forbidden", ["⚠", "→", "✓", "✗", "≥", "≤"])
def test_help_avoids_common_non_cp1251_symbols(forbidden):
    """
    Отдельно — по самым ходовым в этом проекте символам: они щедро
    рассыпаны по логам и комментариям, и перенести такую строку в help
    копипастой очень легко.
    """
    text = _build_parser().format_help()
    assert forbidden not in text, (
        f"символ {forbidden!r} попал в справку argparse; в логах он "
        f"безвреден, в справке — роняет --help на cp1251-консоли"
    )
