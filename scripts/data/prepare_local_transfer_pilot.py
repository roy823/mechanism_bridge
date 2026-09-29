"""Audit fixed reactant panel and prepare reactant-only local-transfer searches."""
import hashlib
import json
from pathlib import Path
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import ArrowLibrary,parse_explicit
from mechbridge.event_graph import graph_smiles,geometry_mol
from mechbridge.reaction_network import atomic_json

# Fixed before PES results; homolog/functional-group coverage, not unseen families.
PANEL = ['CC=O','CCC=O','CCCC=O','CC(C)=O','CCC(C)=O','CC(=O)CC(C)=O',
         'O=CCO','O=CC=O','CC(=O)O','CCC(=O)O','CC(N)=O','CCC(N)=O',
         'CC#N','CCC#N','C[N+](=O)[O-]','CC[N+](=O)[O-]',
         'O=C1CCC1','O=C1CCCC1','CO','CCO','CCC','O=Cc1ccccc1']
SEARCH_STARTS = [('propanal','CCC=O'),('acetone','CC(C)=O')]


def main():
    RDLogger.DisableLog('rdApp.*')
    dest=ROOT/'reports/local_transfer_v3'
    dest.mkdir(parents=True,exist_ok=True)
    if (dest/'coverage.json').exists():raise FileExistsError('Preserve completed coverage audit')
    source=ROOT/'data/raw/synepd/polar.json'
    library=ArrowLibrary(source)
    rows=[]
    for smiles in PANEL:
        mol=parse_explicit(smiles)
        proposals=library.propose(mol,limit=4096)
        rows.append(dict(smiles=graph_smiles(mol),source_graph_present=graph_smiles(mol) in library.source_graphs,
            proposals=len(proposals),distinct_products=len({p['predicted_graph'] for p in proposals}),
            transferred_proposals=sum(p['transferred'] for p in proposals),actions=proposals))
    atomic_json(dest/'coverage.json',dict(library_policy=library.policy,library_audit=dict(library.audit),
        library_rejections=library.rejections,
        library_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        compiled_patterns=len(library.templates),source_graphs=len(library.source_graphs),
        panel=rows,interpretation='Fixed engineering coverage panel; no physical validation or family holdout'))
    starts=[]
    for name,smiles in SEARCH_STARTS:
        mol=parse_explicit(smiles)
        params=AllChem.ETKDGv3();params.randomSeed=20260930
        confs=list(AllChem.EmbedMultipleConfs(mol,numConfs=8,params=params))
        energies=AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=400,numThreads=1)
        converged=[i for i,(status,_) in enumerate(energies) if status==0]
        if not converged:raise RuntimeError('No converged starting conformer')
        best=min(converged,key=lambda i:energies[i][1])
        xyz=mol.GetConformer(confs[best]).GetPositions()
        numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
        if graph_smiles(geometry_mol(numbers,xyz,0)) != graph_smiles(mol):
            raise ValueError('Geometry/graph disagreement')
        starts.append(dict(id=name,atomic_numbers=numbers,positions_A=xyz.tolist(),charge=0,multiplicity=1,
            provenance=dict(reactant_smiles=graph_smiles(mol),geometry='ETKDGv3 + UFF, reactant only',
                embedding_seed=20260930,selected_conformer=int(confs[best]),conformers=len(confs),
                reference_product_or_TS_used=False)))
    output=ROOT/'data/processed/local_transfer_starts.jsonl'
    if output.exists():raise FileExistsError(output)
    output.write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    print(json.dumps(dict(patterns=len(library.templates),audit=dict(library.audit),
        panel=[{k:v for k,v in r.items() if k!='actions'} for r in rows],starts=[s['id'] for s in starts]),indent=2))


if __name__=='__main__':main()
