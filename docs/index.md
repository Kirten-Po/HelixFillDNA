---
layout: default
title: HelixFillDNA
---

<style>
.hero {
  text-align: center;
  padding: 2.5rem 1rem;
}
.hero h1 {
  margin-bottom: 0.3rem;
}
.hero p.tagline {
  color: #555;
  font-size: 1.1rem;
  max-width: 640px;
  margin: 0 auto 1.5rem;
}
.btn-download {
  display: inline-block;
  padding: 0.85rem 1.8rem;
  background: #2a7ae2;
  color: #fff !important;
  border-radius: 6px;
  font-weight: 600;
  text-decoration: none;
}
.btn-download:hover {
  background: #1e5fb8;
}
.screenshot-placeholder {
  border: 2px dashed #ccc;
  border-radius: 8px;
  padding: 3rem 1rem;
  text-align: center;
  color: #888;
  margin: 1.5rem 0;
}
.steps {
  counter-reset: step;
  padding: 0;
  list-style: none;
}
.steps li {
  counter-increment: step;
  position: relative;
  padding-left: 2.8rem;
  margin-bottom: 1.2rem;
}
.steps li::before {
  content: counter(step);
  position: absolute;
  left: 0;
  top: 0;
  width: 2rem;
  height: 2rem;
  border-radius: 50%;
  background: #2a7ae2;
  color: #fff;
  display: flex;
  align-items: center;
  justify-content: center;
  font-weight: 600;
}
</style>

<div class="hero">
  <h1 style="display:flex; align-items:center; justify-content:center; gap:0.6rem;">
    <img src="logo.png" alt="Логотип HelixFillDNA" style="height:1.2em; width:auto; vertical-align:middle;">
    HelixFillDNA
  </h1>
  <p class="tagline">
    Конвертирует сырые данные ДНК (FTDNA, MyHeritage или готовый VCF)
    в формат 23andMe — для импутации и последующей загрузки в Генотек.
  </p>
  <a class="btn-download"
     href="https://github.com/{{ site.repository }}/releases/latest">
    ⬇ Скачать последнюю версию
  </a>
  <p><small>Windows 10+, 64-бит. Установщик не требует прав администратора.</small></p>
</div>

<figure style="margin:1.5rem 0; text-align:center;">
  <img src="main_page.jpg" alt="Вкладка Подготовка: выбор источника данных и референсной панели" style="max-width:100%; border-radius:8px;">
  <figcaption style="color:#777; font-size:0.9rem; margin-top:0.4rem;">
    Вкладка «Подготовка» — выбор источника данных (FTDNA/MyHeritage/VCF), референсной панели импутации и нужных файлов.
  </figcaption>
</figure>

<figure style="margin:1.5rem 0; text-align:center;">
  <img src="main_page2.jpg" alt="Вкладка Подготовка: параметры вывода и порог Rsq" style="max-width:100%; border-radius:8px;">
  <figcaption style="color:#777; font-size:0.9rem; margin-top:0.4rem;">
    Там же — настройка формата вывода, порога качества импутации (Rsq) и параметров донорской панели.
  </figcaption>
</figure>

<figure style="margin:1.5rem 0; text-align:center;">
  <img src="run_page.jpg" alt="Вкладка Запуск: история запусков и старт этапов 1-6" style="max-width:100%; border-radius:8px;">
  <figcaption style="color:#777; font-size:0.9rem; margin-top:0.4rem;">
    Вкладка «Запуск» — история предыдущих запусков и кнопка запуска этапов 1-6 (подготовка файлов для загрузки на сервер импутации).
  </figcaption>
</figure>

<figure style="margin:1.5rem 0; text-align:center;">
  <img src="run_page2.jpg" alt="Этап 2: загрузка результатов Michigan Imputation Server" style="max-width:100%; border-radius:8px;">
  <figcaption style="color:#777; font-size:0.9rem; margin-top:0.4rem;">
    Этап 2 — вставьте curl-команду и пароль из письма Michigan Imputation Server, чтобы скачать результаты и собрать финальный файл.
  </figcaption>
</figure>

## Зачем это нужно

FTDNA и Генотек используют разные чипы — наборы позиций ДНК, которые
они читают. Пересечение между ними — около 49%, а Генотек требует не
меньше 85%, поэтому сырой файл FTDNA он отклоняет напрямую.

HelixFillDNA автоматизирует весь процесс — импутацию (статистическое
восстановление недостающих позиций) и сборку итогового файла по
трафарету — вместо ручного набора команд в Linux/WSL.

## Как начать

<ol class="steps">
  <li>
    <strong>Скачайте установщик</strong> из
    <a href="https://github.com/{{ site.repository }}/releases/latest">GitHub Releases</a>
    и запустите его — права администратора не нужны.
  </li>
  <li>
    <strong>Подготовьте два файла</strong>: ваш сырой экспорт
    (FTDNA/MyHeritage/VCF) и трафарет — экспорт 23andMe, который
    Генотек уже принимал у кого-либо ранее.
  </li>
  <li>
    <strong>Запустите конвертацию</strong> в приложении и загрузите
    готовый файл в Генотек. Полный процесс занимает около двух часов
    (в основном — ожидание закачек), повторный запуск для следующего
    теста — около 20 минут.
  </li>
</ol>

## FAQ

**Безопасно ли это для генетических данных?**
Да. Приложение не передаёт ваши генетические данные разработчику и не
хранит их за пределами вашего устройства. Загрузка на сервер импутации
(Michigan Imputation Server для панели HRC, TOPMed Imputation Server /
BioData Catalyst для панели TOPMed) выполняется напрямую с вашего
компьютера, под вашей собственной учётной записью на этих сервисах.
Подробности — в [политике конфиденциальности](#privacy).

**Нужен ли Python?**
Нет, если вы используете готовый установщик — все зависимости уже
внутри. Python нужен только при запуске из исходного кода.

**Сколько это занимает по времени?**
Первый запуск для нового чипа — около двух часов, из которых
полтора часа — ожидание закачек референсного генома и донорских
образцов (эти файлы скачиваются один раз и переиспользуются дальше).
Повторный запуск для следующего теста на том же чипе — около 20 минут.

## Политика конфиденциальности

<span id="privacy"></span>
Полный текст — в файле
[PRIVACY.md](https://github.com/{{ site.repository }}/blob/main/PRIVACY.md)
в репозитории проекта.

---

<p><small>
© 2026 Поломошнов Кирилл Сергеевич ·
<a href="https://github.com/{{ site.repository }}/blob/main/EULA.md">EULA</a> ·
<a href="https://github.com/{{ site.repository }}/blob/main/LICENSE">Лицензия исходного кода</a> ·
<a href="https://github.com/{{ site.repository }}">Репозиторий на GitHub</a>
</small></p>
