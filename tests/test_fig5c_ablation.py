import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/figures/fig5c_ablation.py'
    spec = importlib.util.spec_from_file_location('fig5c_ablation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ablations_are_tested_against_full_arrows_on_common_systems(tmp_path, monkeypatch):
    fig = load()
    rows = []
    arms = {'arrows': 5, 'arrows@FP-JCTC-1-ablation-no_chain_order': 3,
            'arrows@FP-JCTC-1-ablation-reversed_arrows': 5, 'bond_edits': 6}
    for arm, level in arms.items():
        for k, start in enumerate(['a', 'b', 'c', 'd']):
            if arm == 'bond_edits' and start == 'd':
                continue                      # 'd' lacks a bond_edits run, so it is left out
            rows += [dict(strategy=arm, start=start, seed=s, pairs_at_budget=level + k + s % 2) for s in (17, 29)]
    summary = tmp_path/'summary.json'
    summary.write_text(json.dumps(dict(budget=16000, rows=rows)))
    stats = dict(attempts=10, evaluations_per_attempt=700., connected_fraction=.5,
                 status_fractions=dict(new_connection=.4, ts_force_unconverged=.2))
    outcomes = tmp_path/'outcomes.json'
    outcomes.write_text(json.dumps(dict(per_strategy={a: stats for a in arms})))
    out = tmp_path/'fig5c'
    monkeypatch.setattr('sys.argv', ['f', str(summary), '--out', str(out), '--outcomes', str(outcomes), '--draws', '200'])
    fig.main()
    result = json.loads(out.with_suffix('.json').read_text())
    assert result['systems'] == ['a', 'b', 'c']
    chain = result['versus_arrows']['arrows@FP-JCTC-1-ablation-no_chain_order']
    assert chain['mean_difference'] == -2 and chain['losses'] == 3
    assert result['versus_arrows']['arrows@FP-JCTC-1-ablation-reversed_arrows']['p'] is None
    assert out.with_suffix('.png').exists() and result['attempt_outcomes']['arrows']['ts_force_unconverged'] == .2
    assert fig.holm([.01, None, .04]) == [.02, None, .04]
