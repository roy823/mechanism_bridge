"""Prepare one-sided charged roots for nine glycolysis skeleton inventories."""
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
from mechbridge.symbolic_library import parse_explicit


SPECIES = {
    'glucose': 'O=CC(O)C(O)C(O)C(O)CO',
    'g6p': 'O=CC(O)C(O)C(O)C(O)COP(=O)([O-])O',
    'f6p': 'OCC(=O)C(O)C(O)C(O)COP(=O)([O-])O',
    'fbp': 'O=P([O-])(O)OCC(=O)C(O)C(O)C(O)COP(=O)([O-])O',
    'gap': 'O=C[C@H](O)COP(=O)([O-])O',
    'dhap': 'O=C(COP(=O)([O-])O)CO',
    'bpg13': 'O=C(OP(=O)([O-])O)[C@H](O)COP(=O)([O-])O',
    'pg3': 'O=C(O)[C@H](O)COP(=O)([O-])O',
    'pg2': 'O=C(O)[C@H](OP(=O)([O-])O)CO',
    'pep': 'C=C(OP(=O)([O-])O)C(=O)O',
    'pyruvate': 'CC(=O)C(=O)O',
    'pi': 'O=P([O-])(O)O',
    'water': 'O',
}

ROOTS = [
    dict(id='hexose_isomerization', parts=['g6p']),
    dict(id='fbp_cleavage', parts=['fbp']),
    dict(id='triose_isomerization', parts=['gap']),
    dict(id='phosphoglycerate_shift', parts=['pg3']),
    dict(id='phosphoglycerate_dehydration', parts=['pg2']),
    dict(id='pep_hydrolysis', parts=['pep', 'water']),
    dict(id='glucose_phosphorylation', parts=['glucose', 'pi']),
    dict(id='f6p_phosphorylation', parts=['fbp', 'water']),
    dict(id='bpg_hydrolysis', parts=['bpg13', 'water']),
]


def conformers(fragment, seed, count):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ids = list(AllChem.EmbedMultipleConfs(fragment, numConfs=20, params=params))
    optimized = AllChem.UFFOptimizeMoleculeConfs(fragment, maxIters=700, numThreads=1)
    ranked = sorted((float(energy), cid) for cid, (status, energy) in zip(ids, optimized)
                    if status == 0)
    if not ranked:
        raise ValueError('No converged UFF conformer')
    chosen = []
    for energy, cid in ranked:
        if not chosen or all(AllChem.GetConformerRMS(fragment, cid, old, prealigned=False) > .35
                             for _, old in chosen):
            chosen.append((energy, cid))
        if len(chosen) == count:
            break
    while len(chosen) < count:
        chosen.append(ranked[min(len(chosen), len(ranked)-1)])
    return [(np.asarray(fragment.GetConformer(cid).GetPositions()), energy)
            for energy, cid in chosen]


def constitutional_smiles(mol):
    copy = Chem.RemoveHs(Chem.Mol(mol))
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy, isomericSmiles=False)


def main():
    RDLogger.DisableLog('rdApp.*')
    starts, definitions = [], []
    for root_index, root in enumerate(ROOTS):
        mol = parse_explicit('.'.join(SPECIES[name] for name in root['parts']))
        fragment_ids = Chem.GetMolFrags(mol)
        fragments = Chem.GetMolFrags(mol, asMols=True)
        numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        charge = Chem.GetFormalCharge(mol)
        variants = []
        if len(fragments) == 1:
            for positions, energy in conformers(fragments[0], 91000 + root_index, 2):
                variants.append((positions - positions.mean(0), [energy], None))
        else:
            fragment_positions, energies = [], []
            base = np.zeros((mol.GetNumAtoms(), 3))
            for fragment_index, (atom_ids, fragment) in enumerate(zip(fragment_ids, fragments)):
                positions, energy = conformers(fragment, 91000 + 100 * root_index + fragment_index, 1)[0]
                base[list(atom_ids)] = positions
                energies.append(energy)
            for variant in range(2):
                seed = 92000 + 100 * root_index + variant
                variants.append((assemble_encounter(numbers, base, fragment_ids, seed), energies, seed))
        definitions.append(dict(id=root['id'], parts=root['parts'], input_smiles=graph_smiles(mol),
            formula=rdMolDescriptors.CalcMolFormula(mol), charge=charge, variants=len(variants),
            root_side_only=True, target_geometry_used=False,
            microstate='one negative charge per phosphate; free Pi is H2PO4-'))
        for variant, (positions, energies, encounter_seed) in enumerate(variants):
            perceived = geometry_mol(numbers, positions, charge)
            if (not resonance_equivalent(mol, perceived)
                    and constitutional_smiles(mol) != constitutional_smiles(perceived)):
                raise ValueError(f'Input perception mismatch: {root["id"]}_c{variant}')
            starts.append(dict(id=f'{root["id"]}_c{variant}', atomic_numbers=numbers,
                positions_A=positions.tolist(), charge=charge, multiplicity=1,
                provenance=dict(system=root['id'], parts=root['parts'],
                    input_smiles=graph_smiles(perceived), root_side_only=True,
                    conformer_source='ETKDGv3+UFF', conformer_energies=energies,
                    encounter_seed=encounter_seed, target_TS_or_product_geometry_used=False)))
    starts_path = ROOT / 'data/processed/glycolysis_sparse_endpoint_starts.jsonl'
    starts_path.write_text(''.join(json.dumps(row) + '\n' for row in starts), encoding='utf-8')
    definitions_path = ROOT / 'data/processed/glycolysis_sparse_endpoint_definitions.json'
    definitions_path.write_text(json.dumps(definitions, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(chemical_roots=len(ROOTS), starts=len(starts),
        charges=sorted({row['charge'] for row in starts}), definitions=definitions),
        ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
