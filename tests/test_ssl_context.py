"""
tests/test_ssl_context.py

Регрессия на падение "ssl.SSLError: not enough data: cadata does not
contain a certificate".

Пять попыток скачивания референсного генома тогда не открыли ни одного
сокета: падение происходило в ssl.create_default_context(), ДО
подключения — на Windows он читает системное хранилище сертификатов, и
если хранилище битое, ошибка детерминирована. Поэтому проверяем два
свойства: (1) при битом хранилище контекст всё равно строится, через
явный cafile; (2) ошибка не заворачивается в "проверьте подключение к
интернету" и не уходит в ретраи, где сгорели бы 50 секунд на заведомо
одинаковый результат.
"""
from __future__ import annotations

import ssl

import pytest

import main
from core import network_utils


@pytest.fixture(autouse=True)
def _clear_ssl_cache():
    """make_ssl_context кэширован на процесс — сбрасываем между тестами."""
    network_utils.make_ssl_context.cache_clear()
    yield
    network_utils.make_ssl_context.cache_clear()


def _break_windows_store(monkeypatch):
    """Эмулирует битое хранилище: create_default_context() падает всегда."""
    def boom(*a, **kw):
        if kw.get("cafile"):
            # С явным cafile ветка load_default_certs() не выполняется —
            # ровно то поведение ssl, на которое опирается обход.
            return ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
        raise ssl.SSLError("not enough data: cadata does not contain a certificate")
    monkeypatch.setattr(ssl, "create_default_context", boom)


def test_context_falls_back_to_cafile(monkeypatch, tmp_path):
    ca = tmp_path / "cacert.pem"
    ca.write_bytes(b"x" * (network_utils.CA_BUNDLE_MIN_SIZE + 1))
    monkeypatch.setattr(network_utils, "_LAST_CA_BUNDLE", ca)
    _break_windows_store(monkeypatch)

    ctx = network_utils.make_ssl_context()
    assert isinstance(ctx, ssl.SSLContext)


def test_broken_store_without_bundle_raises_named_error(monkeypatch, tmp_path):
    """
    Нет ни cacert.pem, ни certifi — тогда нужна отдельная, узнаваемая
    ошибка, а не голый RuntimeError: по ней вызывающий код понимает, что
    ретраи бессмысленны.
    """
    monkeypatch.setattr(network_utils, "_LAST_CA_BUNDLE", None)
    monkeypatch.setattr(network_utils, "find_ca_bundle", lambda bin_dir=None: None)
    _break_windows_store(monkeypatch)

    with pytest.raises(network_utils.BrokenCertStoreError):
        network_utils.make_ssl_context()


def test_find_ca_bundle_rejects_truncated_file(tmp_path, monkeypatch):
    """
    Обрезанный .pem как cafile даст ту же SSLError, только в другом
    месте — такой файл хуже, чем его отсутствие.
    """
    monkeypatch.setattr(network_utils, "_LAST_CA_BUNDLE", None)
    monkeypatch.delenv("CURL_CA_BUNDLE", raising=False)
    monkeypatch.delenv("SSL_CERT_FILE", raising=False)
    tiny = tmp_path / "cacert.pem"
    tiny.write_text("обрезано", encoding="utf-8")
    assert network_utils.find_ca_bundle(tmp_path) != tiny


def test_download_does_not_retry_on_ssl_error(monkeypatch, tmp_path):
    """
    Главное свойство фикса: одна попытка вместо пяти, и текст ошибки про
    сертификаты, а не про роутер.
    """
    calls = {"n": 0}

    def boom(url, dest):
        calls["n"] += 1
        raise ssl.SSLError("not enough data: cadata does not contain a certificate")

    monkeypatch.setattr(main, "_download_attempt", boom)
    monkeypatch.setattr(main.time, "sleep", lambda s: None)

    with pytest.raises(RuntimeError) as exc:
        main._download_with_resume(
            ["https://a/x.gz", "https://b/x.gz"], tmp_path / "x.gz",
            max_retries=5,
        )

    assert calls["n"] == 1, "SSL-ошибка детерминирована — повторять её бессмысленно"
    text = str(exc.value)
    assert "сертификат" in text.lower()
    assert "интернет" not in text.lower(), (
        "нельзя отправлять пользователя чинить роутер, когда проблема "
        "в хранилище сертификатов"
    )
    assert "certutil" in text, "в тексте должна быть готовая команда починки"


def test_incomplete_download_still_retries(monkeypatch, tmp_path):
    """
    Обычный обрыв связи, в отличие от SSL-ошибки, ретраить по-прежнему
    нужно — проверяем, что фикс не выключил ретраи вообще.
    """
    calls = {"n": 0}

    def flaky(url, dest):
        calls["n"] += 1
        raise main._IncompleteDownloadError("обрыв")

    monkeypatch.setattr(main, "_download_attempt", flaky)
    monkeypatch.setattr(main.time, "sleep", lambda s: None)

    with pytest.raises(RuntimeError):
        main._download_with_resume("https://a/x.gz", tmp_path / "x.gz", max_retries=3)
    assert calls["n"] == 3
