"""
tests/conftest.py

Гарантирует, что КОРЕНЬ ПРОЕКТА (родитель папки tests/) присутствует в
sys.path — без этого `import adapters`, `import core`, `import template`,
`import main` внутри тестовых модулей падают с ModuleNotFoundError, если
pytest запущен не через `python -m pytest` из корня и без pythonpath в
pyproject.toml/pytest.ini.

pytest автоматически обнаруживает и выполняет conftest.py в каждой
родительской директории коллекции ДО импорта самих тестовых файлов — это
единая точка правды для sys.path, вместо того чтобы дублировать
sys.path.insert(...) в каждом отдельном test_*.py файле (как было сделано
временно в test_assembler_chrom_prefix.py — теперь этот костыль убран
оттуда, см. обновлённую версию файла).

⚠ Если в проекте уже ЕСТЬ conftest.py в tests/ (например, из более ранней
сессии, судя по упоминанию test_main_liftover_and_migration.py в истории
проекта) — этот файл нужно НЕ заменять слепо, а слить с существующим:
проверьте на дублирующиеся fixtures/настройки. Логика ниже идемпотентна
(проверка "not in sys.path" перед insert) и не создаст проблем при
повторном добавлении того же пути, но остальное содержимое существующего
conftest.py (если оно там уже есть — например, общие fixtures для
ChainLiftover/ReferenceGenome) нужно сохранить.
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
