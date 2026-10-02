import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/diagnostics/seed_checks.py'
    spec = importlib.util.spec_from_file_location('seed_checks', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_signature_ignores_rng_seed_and_timing():
    checks = load()
    attempt = dict(proposal=dict(template_id='t', active_atoms=[1, 2], random_seed=1, generation_seconds=.1),
                   status='new_connection', evaluations=900, seconds=3.)
    other = dict(attempt, proposal=dict(attempt['proposal'], random_seed=2, generation_seconds=.2), seconds=4.)
    assert checks.signature(dict(attempts=[attempt])) == checks.signature(dict(attempts=[other]))
    assert checks.signature(dict(attempts=[attempt])) != checks.signature(dict(attempts=[dict(other, evaluations=901)]))


def test_required_seeds_noise_and_heterogeneity():
    checks = load()
    noisy = dict(arrows={'a': [0, 10], 'b': [10, 20]}, bond_edits={'a': [0, 0], 'b': [10, 10]})
    result = checks.required_seeds(noisy, 'arrows', 'bond_edits', systems=10)
    # var 50 + 0, tau^2 = 0, delta = 0.2*7.5: S = ceil(50/(1.5^2*10/2.8016^2)) = 18
    assert result['tau2'] == 0 and result['seeds_required'] == 18 and result['seeds_noise_only'] == 18
    spread = dict(arrows={'a': [0, 0], 'b': [20, 20]}, bond_edits={'a': [0, 0], 'b': [10, 10]})
    result = checks.required_seeds(spread, 'arrows', 'bond_edits', systems=10)
    assert result['heterogeneity_limited'] and result['seeds_required'] is None
    assert result['seeds_noise_only'] == 1


def test_required_seeds_result_is_json_serializable():
    import json
    checks = load()
    noisy = dict(arrows={'a': [0, 10], 'b': [10, 20]}, bond_edits={'a': [0, 0], 'b': [10, 10]})
    json.dumps(checks.required_seeds(noisy, 'arrows', 'bond_edits', systems=10))
