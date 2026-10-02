"""Regression coverage for scheduler-visible pytest failures."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


RUNNER = Path(__file__).resolve().parents[1] / 'scripts/diagnostics/run_pytest_comparison.sh'
pytestmark = pytest.mark.skipif(not shutil.which('bash') or not shutil.which('git'),
                                reason='Comparison runner requires bash and git')


def make_tree(root, name, passing):
    tree = root / name
    (tree / 'tests').mkdir(parents=True)
    (tree / 'tests/test_example.py').write_text(f'def test_example():\n    assert {passing!r}\n')
    subprocess.run(['git', 'init', '-q', str(tree)], check=True)
    subprocess.run(['git', '-C', str(tree), 'add', 'tests/test_example.py'], check=True)
    subprocess.run(['git', '-C', str(tree), '-c', 'user.name=Test',
                    '-c', 'user.email=test@example.invalid', '-c', 'commit.gpgsign=false',
                    'commit', '-qm', 'Fixture'], check=True)
    return tree


def run_comparison(output, *trees):
    env = os.environ.copy()
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
    return subprocess.run(['bash', str(RUNNER), str(output), *map(str, trees)],
                          env=env, text=True, capture_output=True, timeout=60)


def test_failure_survives_a_later_passing_suite(tmp_path):
    bad = make_tree(tmp_path, 'failing tree', False)
    good = make_tree(tmp_path, 'passing tree', True)
    output = tmp_path / 'comparison results'
    result = run_comparison(output, bad, good)
    assert result.returncode == 1, result.stdout + result.stderr
    with (output / 'results.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [r['pytest_exit'] for r in rows] == ['1', '0']
    assert (output / 'tree_1.xml').is_file()
    assert (output / 'tree_2.xml').is_file()
    assert '1 failed' in (output / 'tree_1.log').read_text()
    assert '1 passed' in (output / 'tree_2.log').read_text()


def test_success_and_existing_evidence_are_preserved(tmp_path):
    tree = make_tree(tmp_path, 'passing tree', True)
    output = tmp_path / 'results'
    result = run_comparison(output, tree)
    assert result.returncode == 0, result.stdout + result.stderr
    original = {p.name: p.read_bytes() for p in output.iterdir()}
    assert run_comparison(output, tree).returncode != 0
    assert {p.name: p.read_bytes() for p in output.iterdir()} == original
