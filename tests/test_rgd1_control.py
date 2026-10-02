import importlib.util
import json
from pathlib import Path

import h5py
from rdkit import Chem
from rdkit.Chem import AllChem

from mechbridge.symbolic_library import parse_explicit

ROOT = Path(__file__).resolve().parents[1]
ENOL = '[C:1]([H:4])([H:5])=[C:2]([H:6])[O:3][H:7]'
KETO = '[C:1]([H:4])([H:5])([H:7])[C:2]([H:6])=[O:3]'
ETHANOL = '[C:1]([H:4])([H:5])([H:6])[C:2]([H:7])([H:8])[O:3][H:9]'
ALDEHYDE_H2 = '[C:1]([H:4])([H:5])([H:6])[C:2]([H:7])=[O:3].[H:8][H:9]'


def load():
    path = ROOT/'scripts/data/prepare_rgd1_control.py'
    spec = importlib.util.spec_from_file_location('prepare_rgd1_control', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(tmp_path, rows):
    with h5py.File(tmp_path/'rgd1.h5', 'w') as h5:
        for reaction, reactant, _, ts_energy in rows:
            mol = parse_explicit(reactant)
            order = sorted(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetAtomMapNum())
            mol = Chem.RenumberAtoms(mol, order)          # atom k-1 carries map k, as in RGD1
            assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
            AllChem.MMFFOptimizeMolecule(mol)
            group = h5.create_group(reaction)
            group['elements'] = [a.GetAtomicNum() for a in mol.GetAtoms()]
            group['RG'] = group['TSG'] = mol.GetConformer().GetPositions()
            group['R_E'], group['TS_E'], group['P_E'] = -1., ts_energy, -1.1
    (tmp_path/'rgd1.csv').write_text('reaction,reactant,product\n' +
                                     ''.join(f'{r},{a},{b}\n' for r, a, b, _ in rows))


def test_rgd1_control_frame_dedup_and_mapping(tmp_path, monkeypatch):
    prepare = load()
    write_inputs(tmp_path, [('MR_1', ENOL, KETO, -.9), ('MR_2', ENOL, KETO, -.8),       # TS conformers
                            ('MR_3', ETHANOL, ALDEHYDE_H2, -.7), ('MR_4', ETHANOL, ETHANOL, -.7)])
    out = {k: tmp_path/f'{k}.json' for k in ('starts', 'references', 'report')}
    monkeypatch.setattr('sys.argv', ['prep', '--csv', str(tmp_path/'rgd1.csv'), '--h5', str(tmp_path/'rgd1.h5'),
                                     '--size', '2', *[x for k, v in out.items() for x in (f'--{k}', str(v))]])
    prepare.main()
    report = json.loads(out['report'].read_text())
    assert report['frame_unique_pairs'] == 2
    assert report['csv_rejected'] == {'reactant and product share the same connectivity': 1}
    references = [json.loads(l) for l in out['references'].read_text().splitlines()]
    assert [r['id'] for r in references] == ['rgd1_MR_1', 'rgd1_MR_3']          # MR_2 has the higher TS
    assert [r['stratum'] for r in references] == ['b1f1', 'other']
    enol = references[0]
    assert [2, 6] in enol['reactant_bonds'] and [0, 6] in enol['product_bonds']   # O-H broken, C-H formed
    assert enol['product_graph_source'] == 'csv_mapped_smiles'
