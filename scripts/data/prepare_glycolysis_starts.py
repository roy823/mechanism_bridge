"""Prepare neutral, enzyme-free glycolysis reconstruction inventories.

The physical searches conserve every atom.  Water and phosphoric acid are
explicit surrogate reservoirs; ATP/ADP and NAD+/NADH are not silently removed.
All structures are reactant-side inputs generated without product or TS geometry.
"""
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
    'glucose_open': 'O=CC(O)C(O)C(O)C(O)CO',
    'g6p_open': 'O=CC(O)C(O)C(O)C(O)COP(=O)(O)O',
    'f6p_open': 'OCC(=O)C(O)C(O)C(O)COP(=O)(O)O',
    'fbp_open': 'O=P(O)(O)OCC(=O)C(O)C(O)C(O)COP(=O)(O)O',
    'gap': 'O=C[C@H](O)COP(=O)(O)O',
    'dhap': 'O=C(COP(=O)(O)O)CO',
    'bpg13': 'O=C(OP(=O)(O)O)[C@H](O)COP(=O)(O)O',
    'pg3': 'O=C(O)[C@H](O)COP(=O)(O)O',
    'pg2': 'O=C(O)[C@H](OP(=O)(O)O)CO',
    'pep': 'C=C(OP(=O)(O)O)C(=O)O',
    'pyruvic_acid': 'CC(=O)C(=O)O',
    'water': 'O',
    'phosphoric_acid': 'OP(=O)(O)O',
}


CASES = [
    dict(id='hmp_g6p', parts=['g6p_open'], segment='G6P <-> F6P', evidence='full_atom'),
    dict(id='hmp_f6p', parts=['f6p_open'], segment='G6P <-> F6P', evidence='full_atom'),
    dict(id='cleavage_fbp', parts=['fbp_open'], segment='FBP <-> GAP + DHAP', evidence='full_atom'),
    dict(id='cleavage_gap_dhap', parts=['gap', 'dhap'], segment='FBP <-> GAP + DHAP', evidence='full_atom'),
    dict(id='triose_gap', parts=['gap'], segment='GAP <-> DHAP', evidence='full_atom'),
    dict(id='triose_dhap', parts=['dhap'], segment='GAP <-> DHAP', evidence='full_atom'),
    dict(id='pg_shift_3pg', parts=['pg3'], segment='3PG <-> 2PG', evidence='full_atom'),
    dict(id='pg_shift_2pg', parts=['pg2'], segment='3PG <-> 2PG', evidence='full_atom'),
    dict(id='dehydration_2pg', parts=['pg2'], segment='2PG <-> PEP + H2O', evidence='full_atom'),
    dict(id='dehydration_pep_water', parts=['pep', 'water'], segment='2PG <-> PEP + H2O', evidence='full_atom'),
    dict(id='pep_hydrolysis_pep_water', parts=['pep', 'water'],
         segment='PEP + H2O <-> pyruvate + H3PO4', evidence='phosphate_surrogate'),
    dict(id='pep_hydrolysis_pyruvate_pi', parts=['pyruvic_acid', 'phosphoric_acid'],
         segment='PEP + H2O <-> pyruvate + H3PO4', evidence='phosphate_surrogate'),
    dict(id='glucose_phosphorylation_glucose_pi', parts=['glucose_open', 'phosphoric_acid'],
         segment='glucose + H3PO4 <-> G6P + H2O', evidence='phosphate_surrogate'),
    dict(id='glucose_phosphorylation_g6p_water', parts=['g6p_open', 'water'],
         segment='glucose + H3PO4 <-> G6P + H2O', evidence='phosphate_surrogate'),
    dict(id='f6p_phosphorylation_f6p_pi', parts=['f6p_open', 'phosphoric_acid'],
         segment='F6P + H3PO4 <-> FBP + H2O', evidence='phosphate_surrogate'),
    dict(id='f6p_phosphorylation_fbp_water', parts=['fbp_open', 'water'],
         segment='F6P + H3PO4 <-> FBP + H2O', evidence='phosphate_surrogate'),
    dict(id='bpg_hydrolysis_bpg_water', parts=['bpg13', 'water'],
         segment='1,3-BPG + H2O <-> 3PG + H3PO4', evidence='phosphate_surrogate'),
    dict(id='bpg_hydrolysis_3pg_pi', parts=['pg3', 'phosphoric_acid'],
         segment='1,3-BPG + H2O <-> 3PG + H3PO4', evidence='phosphate_surrogate'),
]


def embed_fragment(fragment, seed):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conformers = list(AllChem.EmbedMultipleConfs(fragment, numConfs=16, params=params))
    optimized = AllChem.UFFOptimizeMoleculeConfs(fragment, maxIters=700, numThreads=1)
    good = [(cid, energy) for cid, (status, energy) in zip(conformers, optimized) if status == 0]
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
    inventories = {}
    for case_index, case in enumerate(CASES):
        smiles = '.'.join(SPECIES[name] for name in case['parts'])
        mol = parse_explicit(smiles)
        fragment_ids = Chem.GetMolFrags(mol)
        fragments = Chem.GetMolFrags(mol, asMols=True)
        xyz = np.zeros((mol.GetNumAtoms(), 3))
        conformer_energies = []
        for fragment_index, (ids, fragment) in enumerate(zip(fragment_ids, fragments)):
            positions, energy = embed_fragment(fragment, 71000 + 100 * case_index + fragment_index)
            xyz[list(ids)] = positions
            conformer_energies.append(energy)
        orientations = 2 if len(fragment_ids) == 2 else 1
        numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        inventory_key = (tuple(sorted(numbers)), Chem.GetFormalCharge(mol))
        inventory_id = inventories.setdefault(inventory_key, f'inventory_{len(inventories):02d}')
        proposals = library.propose(mol)
        definition = dict(case, smiles=smiles, canonical_smiles=graph_smiles(mol),
            formula=rdMolDescriptors.CalcMolFormula(mol), inventory=inventory_id,
            atoms=mol.GetNumAtoms(), fragments=len(fragment_ids), proposal_count=len(proposals),
            proposal_templates=sorted({p['template_id'] for p in proposals}),
            proposal_products=sorted({p['predicted_graph'] for p in proposals}),
            microstate='fully protonated neutral phosphate surrogate',
            enzymes_present=False, ATP_ADP_present=False, NAD_NADH_present=False,
            reference_product_or_TS_geometry_used=False)
        definitions.append(definition)
        for orientation in range(orientations):
            seed = 72000 + 100 * case_index + orientation
            positions = (assemble_encounter(numbers, xyz, fragment_ids, seed)
                         if len(fragment_ids) == 2 else xyz - xyz.mean(0))
            perceived = geometry_mol(numbers, positions, 0)
            if (not resonance_equivalent(mol, perceived)
                    and constitutional_smiles(mol) != constitutional_smiles(perceived)):
                raise ValueError(f'Input perception mismatch: {case["id"]}')
            starts.append(dict(id=f'{case["id"]}_o{orientation}', atomic_numbers=numbers,
                positions_A=positions.tolist(), charge=0, multiplicity=1,
                provenance=dict(system=case['id'], segment=case['segment'],
                    evidence=case['evidence'], inventory=inventory_id,
                    input_smiles=graph_smiles(perceived), parts=case['parts'],
                    microstate=definition['microstate'], conformer_source='ETKDGv3+UFF',
                    conformer_energies=conformer_energies, random_seed=seed,
                    reference_TS_or_product_geometry_used=False)))
    output = ROOT / 'data/processed/glycolysis_starts.jsonl'
    output.write_text(''.join(json.dumps(start) + '\n' for start in starts), encoding='utf-8')
    (ROOT / 'data/processed/glycolysis_definitions.json').write_text(
        json.dumps(definitions, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(dict(starts=len(starts), cases=len(CASES), inventories=len(inventories),
        proposal_counts={row['id']: row['proposal_count'] for row in definitions}),
        ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
