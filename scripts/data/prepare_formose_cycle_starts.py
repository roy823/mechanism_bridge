"""Prepare reactant-only C2/C3/C4 starts spanning the canonical formose cycle."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.encounters import assemble_encounter
from mechbridge.event_graph import geometry_mol,graph_smiles,resonance_equivalent
from mechbridge.symbolic_library import ArrowLibrary,parse_explicit


CASES=[
    dict(id='c2_two_formaldehyde',smiles='C=O.C=O',inventory='C2H4O2',role='initiation: 2 FA -> GO'),
    dict(id='c2_glycolaldehyde',smiles='O=CCO',inventory='C2H4O2',role='glycolaldehyde product/catalyst'),
    dict(id='c2_enediol',smiles='O/C=C/O',inventory='C2H4O2',role='C2 enediol'),
    dict(id='c3_go_fa',smiles='O=CCO.C=O',inventory='C3H6O3',role='GO + FA -> glyceraldehyde'),
    dict(id='c3_glyceraldehyde_d',smiles='O=C[C@H](O)CO',inventory='C3H6O3',role='D-glyceraldehyde'),
    dict(id='c3_dha',smiles='O=C(CO)CO',inventory='C3H6O3',role='dihydroxyacetone'),
    dict(id='c3_enediol',smiles='O/C=C(/O)CO',inventory='C3H6O3',role='LdB-AvE C3 enediol'),
    dict(id='c4_ga_fa',smiles='O=C[C@H](O)CO.C=O',inventory='C4H8O4',role='glyceraldehyde + FA -> tetrose'),
    dict(id='c4_dha_fa',smiles='O=C(CO)CO.C=O',inventory='C4H8O4',role='DHA + FA -> ketotetrose'),
    dict(id='c4_enediol_fa',smiles='O/C=C(/O)CO.C=O',inventory='C4H8O4',role='C3 enediol + FA shortcut'),
    dict(id='c4_erythrose_d',smiles='O=C[C@H](O)[C@H](O)CO',inventory='C4H8O4',role='D-erythrose'),
    dict(id='c4_threose_d',smiles='O=C[C@@H](O)[C@H](O)CO',inventory='C4H8O4',role='D-threose'),
    dict(id='c4_erythrulose_d',smiles='O=C(CO)[C@H](O)CO',inventory='C4H8O4',role='D-erythrulose'),
    dict(id='c4_two_glycolaldehyde',smiles='O=CCO.O=CCO',inventory='C4H8O4',role='retro-aldol target: 2 GO'),
]


def embed_fragment(fragment,seed):
    params=AllChem.ETKDGv3();params.randomSeed=seed
    conformers=list(AllChem.EmbedMultipleConfs(fragment,numConfs=16,params=params))
    optimized=AllChem.UFFOptimizeMoleculeConfs(fragment,maxIters=500,numThreads=1)
    good=[(cid,energy) for cid,(status,energy) in zip(conformers,optimized) if status==0]
    if not good:raise ValueError('No converged reactant conformer')
    cid,_=min(good,key=lambda item:item[1])
    return np.asarray(fragment.GetConformer(cid).GetPositions())


def main():
    RDLogger.DisableLog('rdApp.*')
    library=ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    starts=[];definitions=[]
    for case_index,case in enumerate(CASES):
        mol=parse_explicit(case['smiles']);fragment_ids=Chem.GetMolFrags(mol)
        fragments=Chem.GetMolFrags(mol,asMols=True);xyz=np.zeros((mol.GetNumAtoms(),3))
        for fragment_index,(ids,fragment) in enumerate(zip(fragment_ids,fragments)):
            xyz[list(ids)]=embed_fragment(fragment,51000+100*case_index+fragment_index)
        orientations=2 if len(fragment_ids)==2 else 1
        proposals=library.propose(mol)
        definitions.append(dict(case,atoms=mol.GetNumAtoms(),fragments=len(fragment_ids),
            proposal_count=len(proposals),proposal_products=sorted({p['predicted_graph'] for p in proposals}),
            source='ReactionAtlas canonical formose cycle, arXiv:2606.30778',
            reference_product_or_TS_geometry_used=False))
        for orientation in range(orientations):
            seed=52000+100*case_index+orientation
            positions=(assemble_encounter([a.GetAtomicNum() for a in mol.GetAtoms()],xyz,fragment_ids,seed)
                       if len(fragment_ids)==2 else xyz-xyz.mean(0))
            numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
            perceived=geometry_mol(numbers,positions,0)
            if not resonance_equivalent(mol,perceived):
                raise ValueError(f'Input perception mismatch: {case["id"]}')
            starts.append(dict(id=f"{case['id']}_o{orientation}",atomic_numbers=numbers,
                positions_A=positions.tolist(),charge=0,multiplicity=1,
                provenance=dict(system=case['id'],inventory=case['inventory'],role=case['role'],
                    input_smiles=graph_smiles(perceived),source=definitions[-1]['source'],
                    random_seed=seed,reference_TS_or_product_geometry_used=False)))
    output=ROOT/'data/processed/formose_cycle_starts.jsonl'
    output.write_text(''.join(json.dumps(start)+'\n' for start in starts),encoding='utf-8')
    (ROOT/'data/processed/formose_cycle_definitions.json').write_text(
        json.dumps(definitions,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(starts=len(starts),definitions=definitions),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
