"""
tests/test_truncated_body.py

Регрессия на порчу референса TopMed, найденную по ошибке pyfaidx:
"Line length of fasta file is not consistent! Inconsistent line found in
>chr19 at line 38152550".

В середине файла (смещение 2 708 831 339) оказался XML прокси:
"<Error><Code>ConnectionClosedException</Code><Message>Premature end of
Content-Length delimited message body (expected: 572360640; received:
17508937)...". Прокси оборвал Range-ответ и дописал в ТЕЛО свой отчёт об
ошибке; для read() это обычные данные, они попали в файл, дальше докачка
считала смещение от раздутого размера и продолжала уже за мусором.
Итоговый размер при этом совпал с Content-Length всего файла, поэтому
единственная проверка (final_size != total) порчу пропустила.

Проверяем два рубежа:
  1) короткое тело отдельного соединения обнаруживается сразу, а всё
     записанное этим соединением откатывается (мусор в файл не попадает);
  2) если битый файл уже лежит на диске (наш случай — вместе с
     .sha256-сайдкаром, который совпадает), он опознаётся при проверке
     целостности, а не позже внутри pyfaidx.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import main


class _FakeResponse:
    """Минимальный стенд вместо http-ответа urlopen()."""

    def __init__(self, body: bytes, *, status: int, content_length: int):
        self._body = body
        self._pos = 0
        self.status = status
        self.headers = {"Content-Length": str(content_length)}

    def read(self, size: int = -1) -> bytes:
        if size is None or size < 0:
            size = len(self._body) - self._pos
        chunk = self._body[self._pos:self._pos + size]
        self._pos += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_short_body_rolls_back_and_does_not_keep_garbage(tmp_path, monkeypatch):
    dest = tmp_path / "ref.fa"
    good = b"A" * 1000
    dest.write_bytes(good)

    # Сервер обещает ещё 1000 байт, а отдаёт 20 байт данных и свой XML.
    garbage = b"ACGT" * 5 + b"<Error><Code>ConnectionClosedException</Code></Error>"

    def fake_urlopen(req, timeout=None, context=None):
        assert req.headers.get("Range") == "bytes=1000-"
        return _FakeResponse(garbage, status=206, content_length=1000)

    monkeypatch.setattr(main.urllib.request, "urlopen", fake_urlopen)

    with pytest.raises(main._IncompleteDownloadError) as exc:
        main._download_attempt("https://example.invalid/ref.fa", dest)

    assert "оборвал поток" in str(exc.value)
    # Главное: файл откатился до того, что было скачано корректно.
    assert dest.read_bytes() == good


def test_full_body_is_accepted(tmp_path, monkeypatch):
    dest = tmp_path / "ref.fa"
    dest.write_bytes(b"A" * 1000)
    rest = b"C" * 1000

    def fake_urlopen(req, timeout=None, context=None):
        return _FakeResponse(rest, status=206, content_length=len(rest))

    monkeypatch.setattr(main.urllib.request, "urlopen", fake_urlopen)
    main._download_attempt("https://example.invalid/ref.fa", dest)
    assert dest.read_bytes() == b"A" * 1000 + rest


def test_sha256_scan_finds_injected_error_document(tmp_path):
    path = tmp_path / "ref.fa"
    payload = (
        b">chr1\n"
        + b"ACGT" * 100
        + b"<Error><Code>ConnectionClosedException</Code>"
          b"<Message>Premature end of Content-Length delimited message body"
          b"</Message></Error>"
        + b"ACGT" * 100
    )
    path.write_bytes(payload)

    found: list = []
    main._sha256_of_file(path, injected_out=found)
    assert found, "мусор прокси должен обнаруживаться на проходе SHA-256"
    offset, snippet = found[0]
    assert payload[offset:offset + 7] == b"<Error>"
    assert b"ConnectionClosedException" in snippet


def test_sha256_scan_clean_file_reports_nothing(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_bytes(b">chr1\n" + b"ACGT" * 10_000)
    found: list = []
    digest = main._sha256_of_file(path, injected_out=found)
    assert found == []
    assert len(digest) == 64


def test_marker_spanning_chunk_boundary_is_found(tmp_path, monkeypatch):
    """Мусор может лежать ровно на границе читаемых кусков."""
    path = tmp_path / "ref.fa"
    filler = b"A" * (4 * 1024 * 1024 - 3)
    path.write_bytes(filler + b"<Error>xx" + b"A" * 100)
    found: list = []
    main._sha256_of_file(path, injected_out=found)
    assert found and found[0][0] == len(filler)


def test_injected_error_message_names_file_and_offset(tmp_path):
    path = tmp_path / "GRCh38.fa"
    msg = main._injected_error_message(path, [(2_708_831_339, b"<Error>oops</Error>")])
    assert "GRCh38.fa" in msg
    assert "2708831339" in msg
    assert "<Error>oops</Error>" in msg
