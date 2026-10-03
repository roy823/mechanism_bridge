import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_root_check_subsets_from_saved_rows(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('scorer', ROOT/'scripts/exploration/score_reference_recovery.py')
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    rows = [dict(start=s, reference=s, strategy='arrows', seed=k, S1_mapped=hit, S1_unmapped=hit)
            for s, hit in (('a', True), ('b', False)) for k in range(3)]
    (tmp_path/'old.json').write_text(json.dumps(dict(rows=rows)))
    (tmp_path/'refs.jsonl').write_text('')
    (tmp_path/'rep.json').write_text(json.dumps(dict(representable_ids=['a', 'b'])))
    (tmp_path/'roots.json').write_text(json.dumps(dict(not_preserved=[dict(start='b')])))
    out = tmp_path/'new.json'
    monkeypatch.setattr('sys.argv', ['s', '--from-scores', str(tmp_path/'old.json'), '--references', str(tmp_path/'refs.jsonl'),
                                     '--representable', str(tmp_path/'rep.json'), '--root-check', str(tmp_path/'roots.json'),
                                     '--bootstrap', '200', '--out', str(out)])
    scorer.main()
    result = json.loads(out.read_text())
    assert result['subsets']['representable']['arrows']['S1_mapped_rate'] == .5
    assert result['subsets']['representable_root_preserved']['arrows']['S1_mapped_rate'] == 1.
    assert result['subsets']['all_root_preserved']['arrows']['starts'] == 1
