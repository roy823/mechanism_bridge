import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_overview_draws_without_run_files(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('fig5b', ROOT/'scripts/figures/fig5b_same_edits.py')
    fig = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fig)
    arm = lambda acc, hit, n: dict(attempts=n, accepted=acc, intended=hit)
    summary = [dict(start=f's{k}', group=0, predicted_graph='CC', from_resonance=[k == 1, k == 1],
                    charged_atoms=dict(source=2, predicted=2 + 2*k),
                    per_arm=dict(arrows0=arm(3, 1, 12), arrows1=arm(0, 0, 12), bond_edits=arm(5, 2, 24)))
               for k in range(2)]
    path = tmp_path/'summary.json'
    path.write_text(json.dumps(summary))
    monkeypatch.setattr('sys.argv', ['f', str(path), '--overview', '--out', str(tmp_path/'ov')])
    fig.main()
    assert (tmp_path/'ov.png').exists() and (tmp_path/'ov.pdf').exists()
