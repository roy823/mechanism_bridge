"""Prepare sparse charged glycolysis roots for skeleton reconstruction."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from mechbridge.encounters import assemble_encounter
from mechbridge.event_graph import geometry_mol, graph_smiles, resonance_equivalent
from mechbridge.symbolic_library import ArrowLibrary, parse_explicit


SPECIES = {
    'g6p': 'O=CC(O)C(O)C(O)C(O)COP(=O)([O-])O',
    'f6p': 'OCC(=O)C(O)C(O)C(O)COP(=O)([O-])O',
    'fbp': 'O=P([O-])(O)OCC(=O)C(O)C(O)C(O)COP(=O)([O-])O',
    'gap': 'O=C[C@H](O)COP(=O)([O-])O',
    'dhap': 'O=C(COP(=O)([O-])O)CO',
    'pi': 'O=P([O-])(O)O',
    'water': 'O',
}

CASES = [
    dict(id='g6p', parts=['g6p'], segment='G6P <-> enediol <-> F6P', variants=2),
    dict(id='f6p', parts=['f6p'], segment='G6P <-> enediol <-> F6P', variants=2),
    dict(id='gap', parts=['gap'], segment='GAP <-> enediol <-> DHAP', variants=2),
    dict(id='dhap', parts=['dhap'], segment='GAP <-> enediol <-> DHAP', variants=2),
    dict(id='f6p_pi', parts=['f6p', 'pi'], segment='F6P + H2PO4- <-> FBP + H2O', variants=2),
    dict(id='fbp_water', parts=['fbp', 'water'], segment='F6P + H2PO4- <-> FBP + H2O', variants=2),
]


def embed(fragment, seed):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ids = list(AllChem.EmbedMultipleConfs(fragment, numConfs=12, params=params))
    optimized = AllChem.UFFOptimizeMoleculeConfs(fragment, maxIters=600, numThreads=1)
    good = [(cid, energy) for cid, (status, energy) in zip(ids, optimized) if status == 0]
    if not good:
        raise ValueError('No converged UFF conformer')
    cid, energy = min(good, key=lambda item: item[1])
    return np.asarray(fragment.GetConformer(cid).GetPositions()), float(energy)


def constitutional_smiles(mol):
    copy = Chem.RemoveHs(Chem.Mol(mol))
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy, isomericSmiles=False)


def main():
    RDLogger.DisableLog('rdApp.*')
    library = ArrowLibrary(ROOT / 'data/raw/synepd/polar.json')
    starts, definitions = [], []
    for case_index, case in enumerate(CASES):
        mol = parse_explicit('.'.join(SPECIES[name] for name in case['parts']))
        fragment_ids = Chem.GetMolFrags(mol)
        fragments = Chem.GetMolFrags(mol, asMols=True)
        numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        charge = Chem.GetFormalCharge(mol)
        base = np.zeros((mol.GetNumAtoms(), 3))
        energies = []
        for fragment_index, (atom_ids, fragment) in enumerate(zip(fragment_ids, fragments)):
            positions, energy = embed(fragment, 81000 + 100 * case_index + fragment_index)
            base[list(atom_ids)] = positions
            energies.append(energy)
        proposals = library.propose(mol, limit=64)
        definitions.append(dict(id=case['id'], parts=case['parts'], segment=case['segment'],
            canonical_smiles=graph_smiles(mol), formula=rdMolDescriptors.CalcMolFormula(mol),
            charge=charge, atoms=mol.GetNumAtoms(), variants=case['variants'],
            proposal_templates=sorted({p['template_id'] for p in proposals}),
            microstate='one negative charge per phosphate group; free Pi is H2PO4-',
            atom_conserved=True, charge_conserved=True, multiplicity=1,
            reference_product_or_TS_geometry_used=False))
        for variant in range(case['variants']):
            seed = 82000 + 100 * case_index + variant
            positions = (assemble_encounter(numbers, base, fragment_ids, seed)
                         if len(fragment_ids) == 2 else base - base.mean(0))
            perceived = geometry_mol(numbers, positions, charge)
            if (not resonance_equivalent(mol, perceived)
                    and constitutional_smiles(mol) != constitutional_smiles(perceived)):
                raise ValueError(f'Input perception mismatch: {case["id"]}_{variant}')
            starts.append(dict(id=f'{case["id"]}_c{variant}', atomic_numbers=numbers,
                positions_A=positions.tolist(), charge=charge, multiplicity=1,
                provenance=dict(system=case['id'], segment=case['segment'], parts=case['parts'],
                    input_smiles=graph_smiles(perceived), microstate=definitions[-1]['microstate'],
                    conformer_source='ETKDGv3+UFF', conformer_energies=energies,
                    random_seed=seed, reference_TS_or_product_geometry_used=False)))
    out = ROOT / 'data/processed/glycolysis_charged_sparse_starts.jsonl'
    out.write_text(''.join(json.dumps(row) + '\n' for row in starts), encoding='utf-8')
    definitions_out = ROOT / 'data/processed/glycolysis_charged_sparse_definitions.json'
    definitions_out.write_text(json.dumps(definitions, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(starts=len(starts), cases=len(CASES), charges=sorted({s['charge'] for s in starts}),
        templates={d['id']: d['proposal_templates'] for d in definitions}), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
