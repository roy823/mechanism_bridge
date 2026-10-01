"""Prepare neutral tetrose and 2GO conformer-cluster starts for cycle closure."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import RDLogger
from rdkit.Chem import AllChem

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.reaction_network import aligned_rmsd
from mechbridge.symbolic_library import parse_explicit


TETROSES=[('erythrose','O=C[C@H](O)[C@H](O)CO'),
          ('threose','O=C[C@@H](O)[C@H](O)CO')]
TWO_GO=ROOT/'reports/formose_cycle_v12/recovery/c4_two_go/c4_two_glycolaldehyde_o0/hybrid/network.json'


def tetrose_conformers(label,smiles,seed):
    mol=parse_explicit(smiles);params=AllChem.ETKDGv3();params.randomSeed=seed
    conformers=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=params))
    optimized=AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=700,numThreads=1)
    good=sorted(((cid,energy) for cid,(status,energy) in zip(conformers,optimized) if status==0),
                key=lambda item:item[1])
    selected=[]
    for cid,energy in good:
        positions=np.asarray(mol.GetConformer(cid).GetPositions())
        if all(aligned_rmsd(positions,np.asarray(old['positions_A']))>=.5 for old in selected):
            selected.append(dict(id=f'{label}_uff_c{len(selected)}',atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                positions_A=positions.tolist(),charge=0,multiplicity=1,
                provenance=dict(system='formose_closure',species=label,input_smiles=smiles,
                    conformer_source='ETKDGv3+UFF; heavy/all-atom aligned RMSD selection',UFF_energy=energy,
                    reference_TS_or_product_geometry_used=False)))
        if len(selected)==4:break
    return selected


def main():
    RDLogger.DisableLog('rdApp.*');starts=[]
    for index,(label,smiles) in enumerate(TETROSES):
        starts.extend(tetrose_conformers(label,smiles,61000+index))
    network=json.loads(TWO_GO.read_text(encoding='utf-8'))
    for cluster in network['conformer_clusters']:
        if cluster['graph_smiles']!='O=CCO.O=CCO':continue
        node=network['nodes'][cluster['representative_node']]
        starts.append(dict(id=f"two_go_ml_c{cluster['id']}",atomic_numbers=network['start']['atomic_numbers'],
            positions_A=node['positions_A'],charge=0,multiplicity=1,
            provenance=dict(system='formose_closure',species='2GO',input_smiles='O=CCO.O=CCO',
                conformer_source='AIMNet2-2025 physical minimum from prior 2GO search',
                source_network=TWO_GO.relative_to(ROOT).as_posix(),source_node=node['id'],
                reference_TS_or_product_geometry_used=False)))
    output=ROOT/'data/processed/formose_closure_starts.jsonl'
    output.write_text(''.join(json.dumps(start)+'\n' for start in starts),encoding='utf-8')
    summary=dict(starts=len(starts),by_species={name:sum(s['provenance']['species']==name for s in starts)
        for name in ('erythrose','threose','2GO')})
    (ROOT/'data/processed/formose_closure_starts_summary.json').write_text(
        json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
