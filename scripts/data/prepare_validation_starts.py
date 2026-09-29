"""Keep two independent v1 starts and add two preselected public symbol reactants."""
import json
from pathlib import Path
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import parse_explicit, ArrowLibrary
from mechbridge.event_graph import geometry_mol,graph_smiles


def main():
    RDLogger.DisableLog('rdApp.*')
    source=ROOT/'data/raw/synepd/polar.json'
    records=json.loads(source.read_text(encoding='utf-8'))['records']
    library=ArrowLibrary(source)
    starts=[r for r in map(json.loads,(ROOT/'data/processed/network_starts.jsonl').read_text().splitlines())
            if r['id'] in ('MR_8342_1','MR_460039_0')]
    for sid in (22,1741):
        record=next(r for r in records if r['id']==sid)
        mol=parse_explicit(record['rsmi'].split('>>')[0])
        params=AllChem.ETKDGv3();params.randomSeed=202609
        confs=list(AllChem.EmbedMultipleConfs(mol,numConfs=8,params=params))
        if not confs:raise RuntimeError('No reactant conformer generated')
        energies=AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=400,numThreads=1)
        converged=[i for i,(status,_) in enumerate(energies) if status==0]
        if not converged:raise RuntimeError('No converged UFF starting conformer')
        best=min(converged,key=lambda i:energies[i][1])
        xyz=mol.GetConformer(confs[best]).GetPositions()
        numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
        perceived=geometry_mol(numbers,xyz,0)
        if graph_smiles(perceived)!=graph_smiles(mol) or not library.propose(perceived):
            raise ValueError('Starting conformer or library match failed')
        starts.append(dict(id=f'SynEPD_{sid}',atomic_numbers=numbers,positions_A=xyz.tolist(),
            charge=0,multiplicity=1,provenance=dict(dataset='SynEPD',source_id=sid,
                geometry='RDKit ETKDGv3 then UFF: reactant-only initialization, not quantum geometry',
                selected_conformer=confs[best],embedding_seed=202609,conformers=len(confs),
                reactant_smiles=graph_smiles(mol))))
    path=ROOT/'data/processed/validation_starts.jsonl'
    path.write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    print(json.dumps([dict(id=s['id'],atoms=len(s['atomic_numbers']),provenance=s['provenance']) for s in starts],indent=2))


if __name__=='__main__':main()
