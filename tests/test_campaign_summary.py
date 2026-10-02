import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load():
    path = ROOT/'scripts/exploration/summarize_campaign.py'
    spec = importlib.util.spec_from_file_location('summarize_campaign', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def network():
    nodes = [dict(id=i, graph_smiles=g) for i, g in enumerate(['A', 'B', 'B', 'C', 'D'])]
    edges = [dict(id=0, attempt=0, nodes=[0, 1], kind='chemical'),
             dict(id=1, attempt=1, nodes=[1, 2], kind='conformational'),
             dict(id=2, attempt=2, nodes=[2, 3], kind='chemical'),
             dict(id=3, attempt=3, nodes=[4, 3], kind='chemical')]   # D joins via C
    attempts = [dict(evaluations=100), dict(evaluations=200), dict(evaluations=300), dict(evaluations=50)]
    return dict(nodes=nodes, edges=edges, attempts=attempts, initialization_evaluations=40)


def test_anytime_curve_counts_root_connected_pairs():
    summary = load()
    curve = summary.anytime(network())
    assert curve == [(140, 1), (340, 1), (640, 2), (690, 3)]
    assert summary.value_at(curve, 650) == 2 and summary.value_at(curve, 100) == 0


def test_paired_test_reports_wins():
    summary = load()
    per_system = dict(arrows={'a': 3., 'b': 2., 'c': 5.}, bond_edits={'a': 1., 'b': 2., 'c': 4.})
    result = summary.paired_test(per_system, 'arrows', 'bond_edits')
    assert result['systems'] == 3 and result['wins'] == 2 and result['losses'] == 0


def test_manifest_rows_have_no_carriage_return(tmp_path, monkeypatch):
    path = ROOT/'scripts/exploration/build_campaign_manifest.py'
    spec = importlib.util.spec_from_file_location('build_campaign_manifest', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    out = tmp_path/'campaign'
    monkeypatch.setattr('sys.argv', ['build', '--name', 't', '--starts', str(ROOT/'data/processed/t1x_test_starts.jsonl'),
                                     '--limit', '1', '--strategies', 'geometry', '--seeds', '17',
                                     '--protocol', str(ROOT/'configs/FP-JCTC-1-draft.json'), '--out', str(out)])
    builder.main()
    raw = (out/'manifest.tsv').read_bytes()
    assert b'\r' not in raw
    assert raw.decode().splitlines()[1].split('\t')[-1].endswith('/geometry/s17')
