"""
tests/test_mis_download.py

Регрессии на скачивание архивов результатов Michigan Imputation Server.

Разбор реального отказа у пользователя (три сцепленных бага):

  1. «chr_16.zip качается по 11 КБ/с и то появляется, то пропадает» —
     скачивание здесь строго последовательное, поэтому второй «качающийся»
     файл мог писать только ЧУЖОЙ процесс: осиротевший curl.exe от
     предыдущего запуска, который никто не убивал при закрытии окна.
  2. «[WinError 32] файл занят другим процессом» — тот же осиротевший
     curl держал chr_16.zip, а `dest.unlink()` не был защищён и ронял
     весь Шаг 3 сообщением, по которому невозможно понять, что делать.
  3. «Wrong password» СРАЗУ ПО ВСЕМ 23 архивам при верном пароле — файлы
     от ПРЕДЫДУЩЕГО задания MIS лежали в папке, честно проходили проверку
     «валидный ZIP» и пропускались как уже скачанные; пароль из нового
     письма к ним, естественно, не подходил.

Тесты ниже закрывают все три, плюс докачку и ограничение числа попыток.
"""
from __future__ import annotations

import zipfile

import pytest

import main as pipeline
import mis_adapter
from mis_adapter import MISAdapter, MISAdapterError


# ---------------------------------------------------------------------------
# Фикстуры
# ---------------------------------------------------------------------------
def _make_zip(path, payload=b"data"):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("chr1.dose.vcf.gz", payload)
    return path


def _adapter(tmp_path):
    return MISAdapter(upload_dir=tmp_path / "upload",
                      results_dir=tmp_path / "results")


def _script(urls):
    return "\n".join(
        f"curl -sL {u} -o {u.rsplit('/', 1)[-1]}" for u in urls
    )


@pytest.fixture
def stub_script(monkeypatch):
    """Подменяет получение скрипта скачивания у MIS — без сети."""
    def _install(urls):
        class _Res:
            stdout = _script(urls)
        monkeypatch.setattr(mis_adapter.subprocess, "run",
                            lambda *a, **kw: _Res())
    return _install


# ---------------------------------------------------------------------------
# Расшифровка кодов curl
# ---------------------------------------------------------------------------
def test_curl_exit_code_is_explained():
    """
    Пользователь видел голое "returned non-zero exit status 56" — код
    без единого слова о том, что произошло.
    """
    text = mis_adapter.curl_error_text(56)
    assert "56" in text
    assert "обрыв" in text.lower()

    assert "ссылк" in mis_adapter.curl_error_text(22).lower(), (
        "22 — это истёкшая ссылка, а не проблема со связью"
    )
    assert "сертификат" in mis_adapter.curl_error_text(60).lower()
    # Незнакомый код не должен ломать форматирование.
    assert "999" in mis_adapter.curl_error_text(999)


# ---------------------------------------------------------------------------
# Занятый файл
# ---------------------------------------------------------------------------
def test_locked_file_gives_actionable_error(tmp_path, monkeypatch):
    """
    [WinError 32] сам по себе не говорит НИ ЧТО за процесс, НИ что делать.
    Ошибка должна называть виновника (curl от прошлого запуска) и давать
    инструкцию.
    """
    target = tmp_path / "chr_16.zip"
    target.write_bytes(b"x")

    def boom(self, missing_ok=False):
        raise PermissionError(
            32, "Процесс не может получить доступ к файлу, "
                "так как этот файл занят другим процессом")
    monkeypatch.setattr("pathlib.Path.unlink", boom)

    with pytest.raises(mis_adapter.FileLockedError) as exc:
        mis_adapter._safe_unlink(target)

    text = str(exc.value)
    assert "chr_16.zip" in text
    assert "curl" in text.lower()
    assert "Диспетчер задач" in text


# ---------------------------------------------------------------------------
# Отпечаток задания
# ---------------------------------------------------------------------------
def test_job_id_is_stable_and_order_independent():
    a = ["https://mis.example.org/job1/chr_1.zip", "https://mis.example.org/job1/chr_2.zip"]
    assert mis_adapter.job_id_from_urls(a) == mis_adapter.job_id_from_urls(a[::-1])
    b = ["https://mis.example.org/job2/chr_1.zip", "https://mis.example.org/job2/chr_2.zip"]
    assert mis_adapter.job_id_from_urls(a) != mis_adapter.job_id_from_urls(b)


# ---------------------------------------------------------------------------
# Вердикт по уже лежащему файлу
# ---------------------------------------------------------------------------
def test_manifest_hit_skips_without_network(tmp_path, monkeypatch):
    adapter = _adapter(tmp_path)
    dest = _make_zip(adapter.results_dir / "chr_1.zip")
    url = "https://mis.example.org/job1/chr_1.zip"
    known = {"chr_1.zip": {"url": url, "size": dest.stat().st_size}}

    monkeypatch.setattr(mis_adapter, "remote_size", lambda u: pytest.fail(
        "при совпадении манифеста сеть трогать не нужно"))
    assert adapter._existing_file_verdict(dest, url, known) == "keep"


def test_file_from_other_job_is_rejected(tmp_path):
    """Главная регрессия: архив от другого задания — не «уже скачанный»."""
    adapter = _adapter(tmp_path)
    dest = _make_zip(adapter.results_dir / "chr_1.zip")
    known = {"chr_1.zip": {"url": "https://mis.example.org/СТАРОЕ_задание/chr_1.zip",
                           "size": dest.stat().st_size}}
    verdict = adapter._existing_file_verdict(
        dest, "https://mis.example.org/новое_задание/chr_1.zip", known,
    )
    assert verdict != "keep"
    assert "задание" in verdict


def test_unknown_file_adopted_when_size_matches_server(tmp_path, monkeypatch):
    """
    Файлы от предыдущей версии программы манифеста не имеют. Перекачивать
    гигабайты вслепую нельзя, поэтому сверяем длину с сервером.
    """
    adapter = _adapter(tmp_path)
    dest = _make_zip(adapter.results_dir / "chr_1.zip")
    size = dest.stat().st_size

    monkeypatch.setattr(mis_adapter, "remote_size", lambda u: size)
    assert adapter._existing_file_verdict(dest, "https://mis.example.org/j/chr_1.zip", {}) == "adopt"

    monkeypatch.setattr(mis_adapter, "remote_size", lambda u: size + 1000)
    verdict = adapter._existing_file_verdict(dest, "https://mis.example.org/j/chr_1.zip", {})
    assert verdict != "adopt" and "размер" in verdict


def test_truncated_archive_is_rejected(tmp_path):
    adapter = _adapter(tmp_path)
    dest = adapter.results_dir / "chr_1.zip"
    dest.write_bytes(b"PK\x03\x04 truncated")
    verdict = adapter._existing_file_verdict(dest, "https://mis.example.org/j/chr_1.zip", {})
    assert "повреждён" in verdict


# ---------------------------------------------------------------------------
# Цикл скачивания целиком
# ---------------------------------------------------------------------------
def test_stale_job_archives_are_redownloaded(tmp_path, stub_script, monkeypatch):
    """
    Сквозная проверка бага «Wrong password по всем 23»: в папке лежат
    целые архивы от прошлого задания, приходит скрипт нового задания —
    все файлы должны быть перекачаны, а не пропущены.
    """
    adapter = _adapter(tmp_path)
    old_urls = [f"https://mis.example.org/OLD/chr_{i}.zip" for i in (1, 2)]
    new_urls = [f"https://mis.example.org/NEW/chr_{i}.zip" for i in (1, 2)]

    # Имитируем состояние после скачивания старого задания.
    for url in old_urls:
        name = url.rsplit("/", 1)[-1]
        _make_zip(adapter.results_dir / name, "старое задание".encode("utf-8"))
    mis_adapter.write_manifest(adapter.results_dir, {
        "job": mis_adapter.job_id_from_urls(old_urls),
        "files": {url.rsplit("/", 1)[-1]: {"url": url, "size": 1}
                  for url in old_urls},
    })

    downloaded = []

    def fake_download(url, dest, cancel_check=None):
        downloaded.append(url)
        _make_zip(dest, "новое задание".encode("utf-8"))

    monkeypatch.setattr(mis_adapter, "download_one", fake_download)
    monkeypatch.setattr(mis_adapter, "remote_size",
                        lambda u: pytest.fail("должно решаться по манифесту"))
    stub_script(new_urls)

    adapter.download_results("curl -sL https://mis.example.org/get/NEW")

    assert sorted(downloaded) == sorted(new_urls), (
        "архивы от другого задания обязаны быть перекачаны"
    )
    manifest = mis_adapter.read_manifest(adapter.results_dir)
    assert manifest["job"] == mis_adapter.job_id_from_urls(new_urls)


def test_same_job_second_run_downloads_nothing(tmp_path, stub_script, monkeypatch):
    """Повторный запуск того же задания не должен качать заново."""
    adapter = _adapter(tmp_path)
    urls = [f"https://mis.example.org/JOB/chr_{i}.zip" for i in (1, 2)]

    def fake_download(url, dest, cancel_check=None):
        _make_zip(dest, b"payload")

    monkeypatch.setattr(mis_adapter, "download_one", fake_download)
    stub_script(urls)
    adapter.download_results("curl -sL https://mis.example.org/get/JOB")

    calls = []
    monkeypatch.setattr(mis_adapter, "download_one",
                        lambda *a, **kw: calls.append(a))
    adapter.download_results("curl -sL https://mis.example.org/get/JOB")
    assert calls == []


def test_retries_are_bounded(tmp_path, stub_script, monkeypatch):
    """
    Раньше цикл повторов был `while True`: пользователь мог до
    бесконечности жать «Да» в диалоге, каждый раз начиная файл с нуля.
    """
    adapter = _adapter(tmp_path)
    urls = ["https://mis.example.org/JOB/chr_1.zip"]
    attempts = {"n": 0}

    def always_fails(url, dest, cancel_check=None):
        attempts["n"] += 1
        raise MISAdapterError("обрыв связи")

    monkeypatch.setattr(mis_adapter, "download_one", always_fails)
    stub_script(urls)

    with pytest.raises(MISAdapterError) as exc:
        adapter.download_results(
            "curl -sL https://mis.example.org/get/JOB",
            on_file_error=lambda name, err: True,   # пользователь всегда «Да»
        )
    assert attempts["n"] == mis_adapter.MAX_FILE_ATTEMPTS
    assert "с места обрыва" in str(exc.value).lower()


def test_cancel_stops_the_loop(tmp_path, stub_script, monkeypatch):
    adapter = _adapter(tmp_path)
    stub_script([f"https://mis.example.org/JOB/chr_{i}.zip" for i in (1, 2, 3)])
    monkeypatch.setattr(mis_adapter, "download_one",
                        lambda *a, **kw: pytest.fail("не должно дойти до закачки"))
    with pytest.raises(mis_adapter.DownloadCancelled):
        adapter.download_results("curl -sL https://mis.example.org/get/JOB",
                                 cancel_check=lambda: True)


# ---------------------------------------------------------------------------
# Один файл: .part и переименование
# ---------------------------------------------------------------------------
def test_partial_download_never_becomes_the_final_name(tmp_path, monkeypatch):
    """
    Раньше curl писал сразу в финальное имя, и оборванная закачка
    оставалась на диске под именем настоящего архива.
    """
    dest = tmp_path / "chr_1.zip"

    def fake_curl(args, cancel_check=None):
        out = tmp_path / (dest.name + ".part")
        out.write_bytes(b"PK\x03\x04 half a file")
        return 0, ""

    monkeypatch.setattr(mis_adapter, "_run_curl", fake_curl)
    with pytest.raises(MISAdapterError) as exc:
        mis_adapter.download_one("https://mis.example.org/JOB/chr_1.zip", dest)

    assert not dest.exists(), "битый файл не должен получить финальное имя"
    assert not (tmp_path / "chr_1.zip.part").exists()
    assert "ZIP" in str(exc.value)


def test_successful_download_renames_part(tmp_path, monkeypatch):
    dest = tmp_path / "chr_1.zip"

    def fake_curl(args, cancel_check=None):
        _make_zip(tmp_path / (dest.name + ".part"))
        return 0, ""

    monkeypatch.setattr(mis_adapter, "_run_curl", fake_curl)
    mis_adapter.download_one("https://mis.example.org/JOB/chr_1.zip", dest)
    assert dest.exists() and mis_adapter._is_valid_zip(dest)
    assert not (tmp_path / "chr_1.zip.part").exists()


def test_curl_flags_protect_against_stalled_transfer():
    """
    Соединение на 11 КБ/с должно обрываться и переподключаться, а не
    тянуть 250 МБ шесть часов; HTTP-ошибка должна быть ошибкой, а не
    молча скачанной страницей вместо архива.
    """
    args = mis_adapter._curl_base_args()
    joined = " ".join(args)
    assert "--speed-limit" in joined and "--speed-time" in joined
    assert "-sSfL" in args, "-f обязателен, -S возвращает текст ошибки"
    assert "--connect-timeout" in joined


def test_resume_is_used_when_part_exists(tmp_path, monkeypatch):
    dest = tmp_path / "chr_1.zip"
    part = tmp_path / "chr_1.zip.part"
    part.write_bytes(b"x" * 1000)
    seen = {}

    def fake_curl(args, cancel_check=None):
        seen["args"] = args
        _make_zip(part)
        return 0, ""

    monkeypatch.setattr(mis_adapter, "_run_curl", fake_curl)
    mis_adapter.download_one("https://mis.example.org/JOB/chr_1.zip", dest)
    assert "-C" in seen["args"], "оборванная закачка должна докачиваться"


def test_server_without_resume_restarts_once(tmp_path, monkeypatch):
    """Код 33 = докачка не поддерживается: качаем заново, но без цикла."""
    dest = tmp_path / "chr_1.zip"
    part = tmp_path / "chr_1.zip.part"
    part.write_bytes(b"x" * 1000)
    calls = {"n": 0}

    def fake_curl(args, cancel_check=None):
        calls["n"] += 1
        if "-C" in args:
            return 33, "curl: (33) Range requests not supported"
        _make_zip(part)
        return 0, ""

    monkeypatch.setattr(mis_adapter, "_run_curl", fake_curl)
    mis_adapter.download_one("https://mis.example.org/JOB/chr_1.zip", dest)
    assert calls["n"] == 2
    assert dest.exists()


# ---------------------------------------------------------------------------
# Подсказка про пароль
# ---------------------------------------------------------------------------
def test_password_hint_added_only_for_password_errors(tmp_path):
    wrong = pipeline._password_failure_hint(
        "ERROR: Wrong password : chr1.dose.vcf.gz", tmp_path,
    )
    assert "разных заданий" in wrong or "разных задани" in wrong
    assert str(tmp_path) in wrong

    other = "tabix не смог проиндексировать файл"
    assert pipeline._password_failure_hint(other, tmp_path) == other
