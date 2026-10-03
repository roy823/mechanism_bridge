import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WATER = [[0., 0., 0.], [0.96, 0., 0.], [-0.24, 0.93, 0.]]


def load():
    path = ROOT/'scripts/diagnostics/root_preservation.py'
    spec = importlib.util.spec_from_file_location('root_preservation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_root_changes_are_flagged(tmp_path, monkeypatch):
    check = load()
    start = dict(id='w', atomic_numbers=[8, 1, 1], charge=0, positions_A=WATER)
    key = check.connectivity(start, WATER)[0]
    for k, positions in enumerate([WATER, [[0., 0., 0.], [0.96, 0., 0.], [-2.4, 9.3, 0.]]]):
        run = tmp_path/f'campaign{k}'/'runs'/'w'/'arrows'/'s17'/'w'/'arrows'
        run.mkdir(parents=True)
        (run/'network.json').write_text(json.dumps(dict(status='completed', start=dict(start, id=f'w{k}'),
                                                        nodes=[dict(positions_A=positions)])))
    refs = tmp_path/'refs.jsonl'
    refs.write_text('\n'.join(json.dumps(dict(id=f'w{k}', atomic_numbers=[8, 1, 1], reactant_bonds=[[0, 1], [0, 2]],
                                              reactant_key=key, admission='accepted')) for k in range(2)))
    out = tmp_path/'out.json'
    monkeypatch.setattr('sys.argv', ['r', str(tmp_path/'campaign0'), str(tmp_path/'campaign1'),
                                     '--references', str(refs), '--out', str(out)])
    check.main()
    result = json.loads(out.read_text())
    assert result['preserved'] == 1 and [r['start'] for r in result['not_preserved']] == ['w1']
    assert result['reactions'] == 2 and result['reactions_root_mapped'] == 1
