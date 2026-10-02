import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/diagnostics/summarize_same_edits.py'
    spec = importlib.util.spec_from_file_location('summarize_same_edits', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_holm_keeps_order_and_missing_values():
    summary = load()
    assert summary.holm([.01, .04, None, .03]) == pytest.approx([.03, .06, None, .06])


def test_mcnemar_pairs_by_sample_and_seed():
    summary = load()
    ok, bad = 'validated_descents', 'ts_force_unconverged'
    arrows = [dict(sample=0, seed=17, status=ok), dict(sample=0, seed=29, status=ok)]
    edits = [dict(sample=0, seed=17, status=ok), dict(sample=0, seed=29, status=bad),
             dict(sample=0, seed=1059, status=ok)]                     # unpaired extra control seed
    result = summary.mcnemar(arrows, edits)
    assert result['pairs'] == 2 and result['arrows_only'] == 1 and result['bond_edits_only'] == 0
    assert result['p'] == pytest.approx(1.)


def test_permutation_test_separated_clusters():
    summary = load()
    p = summary.permutation_p(['a']*3 + ['b']*3, [0]*3 + [1]*3, draws=2000)
    assert 0 < p < .2                                                   # 2 of 20 labelings separate
    assert summary.permutation_p(['a', 'b'], [0, 0]) is None
