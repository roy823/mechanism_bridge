import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/diagnostics/summarize_replays.py'
    spec = importlib.util.spec_from_file_location('summarize_replays', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows(group, variant, accepted, n, evaluations):
    return [dict(group=group, attempt=f'{group}/{k}', variant=variant, evaluations=evaluations,
                 status='validated_descents' if k < accepted else 'ts_force_unconverged',
                 historical_match=True) for k in range(n)]


def test_p1c_rule_needs_a_large_gain_without_small_losses():
    summary = load()
    data = []
    for variant, large, small, cost in (('legacy', 1, 5, 1000), ('p1c_handoff03', 4, 5, 1050),
                                        ('p1c_handoff10', 5, 3, 900), ('p1c_long_dimer', 1, 5, 1000)):
        data += rows('P1c_large', variant, large, 10, 1800) + rows('P1c_small', variant, small, 10, cost)
    result = summary.p1c_summary(data)
    # handoff10 gains more on large systems but loses 20 points on small ones.
    assert result['qualifying'] == ['p1c_handoff03'] and result['decision'] == 'p1c_handoff03'
    assert result['large']['legacy_historical_match'] == 10
