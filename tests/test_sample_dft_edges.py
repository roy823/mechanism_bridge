import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def campaign(root, name, protocol, strategy):
    camp = root/name
    run = camp/'runs'/f'sys/{strategy}/s17'
    net = run/'sys'/strategy
    (net/'attempt_000').mkdir(parents=True)
    (net/'attempt_000'/'result.json').write_text(json.dumps(dict(status='validated_descents')))
    network = dict(nodes=[dict(graph_smiles='A'), dict(graph_smiles='B')],
                   attempts=[dict(artifact='attempt_000/result.json')],
                   edges=[dict(id=0, attempt=0, nodes=[0, 1], kind='chemical', ts_energy_eV=-1.,
                               ts_positions_A=[[0., 0., 0.], [0., 0., 1.]])])
    (net/'network.json').write_text(json.dumps(network))
    (camp/'manifest.tsv').write_text('task_id\tstart_file\tstart_id\tstrategy\tseed\tpotential\toutdir\n'
                                     f'1\tx\tsys\t{strategy}\t17\taimnet2-rxn\tsys/{strategy}/s17\n')
    (camp/'campaign.json').write_text(json.dumps(dict(protocol=protocol)))
    return camp


def test_strata_label_non_main_protocols(tmp_path):
    path = ROOT/'scripts/qc/sample_dft_edges.py'
    spec = importlib.util.spec_from_file_location('sample_dft_edges', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    main = campaign(tmp_path, 'C1', 'configs/FP-JCTC-1.json', 'arrows')
    legacy = campaign(tmp_path, 'C7', 'configs/FP-JCTC-1-arrows-legacy.json', 'arrows')
    edges, kept = module.frame([main, legacy])
    assert [e['stratum'] for e in edges] == ['arrows', 'arrows@FP-JCTC-1-arrows-legacy']
    assert [e['strategy'] for e in edges] == ['arrows', 'arrows']      # run directories keep the strategy
    assert len(kept) == 2                                             # duplicates are removed within a stratum only
    again, kept = module.frame([main, main])
    assert len(again) == 2 and len(kept) == 1
