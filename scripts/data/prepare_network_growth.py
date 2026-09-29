"""Prepare a fixed CHNO panel using only reactant structures and public arrows."""
import json
from pathlib import Path
import sys
from rdkit import RDLogger
from rdkit.Chem import AllChem

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import ArrowLibrary,parse_explicit
from mechbridge.event_graph import graph_smiles,geometry_mol

SYSTEMS=[('acetone','CC(C)=O'),('propanal','CCC=O'),('glycolaldehyde','O=CCO'),
         ('glyceraldehyde','O=CC(O)CO'),('acetylacetone','CC(=O)CC(C)=O'),
         ('cyclobutanone','O=C1CCC1'),('acetamide','CC(N)=O'),('nitroethane','CC[N+](=O)[O-]')]


def main():
    RDLogger.DisableLog('rdApp.*')
    path=ROOT/'data/processed/network_growth_starts.jsonl'
    if path.exists():raise FileExistsError(path)
    library=ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    starts=[]
    for name,smiles in SYSTEMS:
        mol=parse_explicit(smiles)
        params=AllChem.ETKDGv3();params.randomSeed=20261001
        confs=list(AllChem.EmbedMultipleConfs(mol,numConfs=8,params=params))
        energies=AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=400,numThreads=1)
        good=[i for i,(status,_) in enumerate(energies) if status==0]
        if not good:raise ValueError(f'No UFF conformer: {name}')
        best=min(good,key=lambda i:energies[i][1])
        xyz=mol.GetConformer(confs[best]).GetPositions()
        numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
        perceived=geometry_mol(numbers,xyz,0)
        proposals=library.propose(perceived)
        # Unspecified input stereocenters are realized by the generated geometry.
        starts.append(dict(id=name,atomic_numbers=numbers,positions_A=xyz.tolist(),charge=0,multiplicity=1,
            provenance=dict(input_smiles=smiles,realized_smiles=graph_smiles(perceived),
                geometry='ETKDGv3/UFF reactant-only ensemble',embedding_seed=20261001,
                selected_conformer=int(confs[best]),conformers=len(confs),symbolic_proposals=len(proposals),
                hypothesis_product_graphs=sorted({p['predicted_graph'] for p in proposals}),
                reference_TS_or_product_geometry_used=False)))
    path.write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    print(json.dumps([dict(id=s['id'],atoms=len(s['atomic_numbers']),proposals=s['provenance']['symbolic_proposals'])
                      for s in starts],indent=2))


if __name__=='__main__':main()
