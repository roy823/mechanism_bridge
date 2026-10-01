"""Build small reactant-only 3D starts for reaction classes shown by FlowER.

The molecules are minimal representatives of the published reaction classes, not
the exact proprietary NameRxn evaluation records. No product or TS geometry is
used to construct a start.
"""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from mechbridge.encounters import assemble_encounter
from mechbridge.event_graph import geometry_mol, graph_smiles, resonance_equivalent
from mechbridge.symbolic_library import ArrowLibrary, parse_explicit


CASES = [
    dict(id='flower_diels_alder_minimal', reactant='C=CC=C.C=C',
         reaction_class='Diels-Alder cycloaddition',
         expected_products=['C1=CCCCC1'], paper_location='Fig. 4d and Fig. S16'),
    dict(id='flower_prilezhaev_minimal', reactant='C=C.OOC=O',
         reaction_class='Prilezhaev epoxidation',
         expected_products=['C1CO1.O=CO'], paper_location='Fig. 4d and Fig. S16'),
    dict(id='flower_transamidation_minimal', reactant='CC(=O)N.CN',
         reaction_class='Transamidation',
         expected_products=['CNC(C)=O.N'], paper_location='Fig. 4 and Fig. S16'),
]


def embedded_fragment(fragment, seed, prefer_s_cis=False):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conformers = list(AllChem.EmbedMultipleConfs(fragment, numConfs=24, params=params))
    optimized = AllChem.UFFOptimizeMoleculeConfs(fragment, maxIters=500, numThreads=1)
    usable = [(cid, energy) for cid, (status, energy) in zip(conformers, optimized) if status == 0]
    if not usable:
        raise ValueError('No converged reactant conformer')
    policy = 'lowest_UFF_energy_reactant_conformer'
    if prefer_s_cis:
        query = Chem.MolFromSmarts('[C]=[C]-[C]=[C]')
        match = fragment.GetSubstructMatch(query)
        if match:
            # A reactant-only conformational prior exposes both diene termini.
            window = min(e for _, e in usable) + 5.0
            pool = [(cid, e) for cid, e in usable if e <= window]
            cid, _ = min(pool, key=lambda item: np.linalg.norm(
                fragment.GetConformer(item[0]).GetAtomPosition(match[0]) -
                fragment.GetConformer(item[0]).GetAtomPosition(match[3])))
            policy = 'lowest_terminal_distance_within_5_kcal_mol_UFF_window'
            return np.asarray(fragment.GetConformer(cid).GetPositions()), policy
    cid, _ = min(usable, key=lambda item: item[1])
    return np.asarray(fragment.GetConformer(cid).GetPositions()), policy


def main():
    RDLogger.DisableLog('rdApp.*')
    output = ROOT / 'data/processed/flower_case_starts.jsonl'
    case_file = ROOT / 'data/processed/flower_case_definitions.json'
    library = ArrowLibrary(ROOT / 'data/raw/synepd/polar.json')
    starts = []
    definitions = []
    for case_index, case in enumerate(CASES):
        mol = parse_explicit(case['reactant'])
        fragment_ids = Chem.GetMolFrags(mol)
        fragments = Chem.GetMolFrags(mol, asMols=True)
        xyz = np.zeros((mol.GetNumAtoms(), 3))
        conformer_policies = []
        for fragment_index, (ids, fragment) in enumerate(zip(fragment_ids, fragments)):
            prefer_s_cis = case['id'].startswith('flower_diels_alder') and fragment.GetNumHeavyAtoms() == 4
            positions, policy = embedded_fragment(
                fragment, 41000 + 100 * case_index + fragment_index, prefer_s_cis)
            xyz[list(ids)] = positions
            conformer_policies.append(policy)
        proposals = library.propose(mol)
        definitions.append(dict(case,
            role='minimal_representative_of_published_reaction_class',
            source='https://doi.org/10.1038/s41586-025-09426-9',
            exact_NameRxn_record_available=False,
            symbolic_proposals=[dict(template_id=p['template_id'],
                                     predicted_graph=p['predicted_graph']) for p in proposals]))
        for orientation in range(3):
            seed = 42000 + 100 * case_index + orientation
            numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
            encounter = assemble_encounter(numbers, xyz, fragment_ids, seed)
            perceived = geometry_mol(numbers, encounter, 0)
            if not resonance_equivalent(mol, perceived):
                raise ValueError(f'Input perception mismatch for {case["id"]}')
            starts.append(dict(
                id=f'{case["id"]}_o{orientation}',
                atomic_numbers=numbers, positions_A=encounter.tolist(),
                charge=0, multiplicity=1,
                provenance=dict(system=case['id'], input_smiles=graph_smiles(perceived),
                    reaction_class=case['reaction_class'], paper_location=case['paper_location'],
                    source='https://doi.org/10.1038/s41586-025-09426-9',
                    representative_scope='minimal_same-class_system; not exact paper evaluation molecule',
                    initial_fragment_atom_ids=[list(ids) for ids in fragment_ids],
                    conformer_policies=conformer_policies,
                    assembly='independent uniform SO(3) rotations, vdW gap 0.2 A',
                    random_seed=seed, reference_TS_or_product_geometry_used=False)))
    output.write_text(''.join(json.dumps(start) + '\n' for start in starts), encoding='utf-8')
    case_file.write_text(json.dumps(definitions, indent=2), encoding='utf-8')
    print(json.dumps(dict(starts=len(starts), cases=definitions), indent=2))


if __name__ == '__main__':
    main()
