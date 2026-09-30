"""Reactant-only encounters inspired by ReactionAtlas and actual Coley records."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.symbolic_library import ArrowLibrary,parse_explicit
from mechbridge.encounters import assemble_encounter
from mechbridge.event_graph import graph_smiles,geometry_mol,resonance_equivalent


def main():
    RDLogger.DisableLog('rdApp.*')
    destination=ROOT/'data/processed/bimolecular_starts.jsonl'
    if destination.exists():raise FileExistsError(destination)
    library=ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    records=list(csv.DictReader((ROOT/'data/raw/coley_dipolar/full_dataset.csv').open()))
    eligible=[]
    for row in records:
        mol=parse_explicit(row['rxn_smiles'].split('>>')[0])
        if (set(a.GetAtomicNum() for a in mol.GetAtoms())<={1,6,7,8}
            and len(Chem.GetMolFrags(mol))==2 and Chem.GetFormalCharge(mol)==0
            and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())):
            eligible.append((mol.GetNumAtoms(),int(row['rxn_id']),graph_smiles(mol),row))
    selected=[];seen=set()
    for _,_,smi,row in sorted(eligible):
        if smi not in seen:
            selected.append(row);seen.add(smi)
        if len(selected)==2:break
    systems=[('formaldehyde_water','C=O.O','ReactionAtlas seed pair'),
        ('formaldehyde_dimer','C=O.C=O','ReactionAtlas formose initiation'),
        ('formaldehyde_glycolaldehyde','C=O.O=CCO','ReactionAtlas formose building blocks'),
        ('formaldehyde_enediol','C=O.OC=CO','ReactionAtlas formose enediol context')]
    systems += [('coley_'+r['rxn_id'],r['rxn_smiles'].split('>>')[0],
                 'Coley dipolar dataset: reactants from record '+r['rxn_id']) for r in selected]
    starts=[];coverage=[]
    for index,(name,smi,source) in enumerate(systems):
        mol=parse_explicit(smi)
        for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
        fragments=Chem.GetMolFrags(mol)
        # Embed fragments independently; no product geometry is ever generated.
        x=np.zeros((mol.GetNumAtoms(),3))
        for ids,part in zip(fragments,Chem.GetMolFrags(mol,asMols=True)):
            params=AllChem.ETKDGv3();params.randomSeed=23000+index
            conformers=list(AllChem.EmbedMultipleConfs(part,numConfs=4,params=params))
            optimized=AllChem.UFFOptimizeMoleculeConfs(part,maxIters=500,numThreads=1)
            good=[i for i,(status,_) in enumerate(optimized) if status==0]
            if not good:raise ValueError('No converged reactant conformer: '+name)
            best=min(good,key=lambda i:optimized[i][1])
            x[list(ids)]=part.GetConformer(conformers[best]).GetPositions()
        numbers=[a.GetAtomicNum() for a in mol.GetAtoms()]
        proposals=library.propose(mol)
        coverage.append(dict(system=name,graph=graph_smiles(mol),atoms=len(numbers),
            proposals=len(proposals),intermolecular=sum(p['intermolecular'] for p in proposals),
            actions=[dict(name=p['name'],origin=p['origin'],product=p['predicted_graph']) for p in proposals]))
        for orientation in range(2):
            seed=27000+index*10+orientation
            xyz=assemble_encounter(numbers,x,fragments,seed)
            actual=geometry_mol(numbers,xyz,0)
            connectivity=Chem.Mol(actual);expected=Chem.Mol(mol)
            Chem.RemoveStereochemistry(connectivity);Chem.RemoveStereochemistry(expected)
            if not resonance_equivalent(expected,connectivity):raise ValueError('Encounter changed graph: '+name)
            starts.append(dict(id=f'{name}_o{orientation}',atomic_numbers=numbers,positions_A=xyz.tolist(),
                charge=0,multiplicity=1,provenance=dict(system=name,input_smiles=graph_smiles(actual),
                source=source,initial_fragment_atom_ids=[list(f) for f in fragments],
                source_reactant_smiles=graph_smiles(mol),
                resonance_representation_changed=graph_smiles(connectivity)!=graph_smiles(expected),
                assembly='independent uniform SO(3) rotations, vdW gap 0.2 A',random_seed=seed,
                reference_TS_or_product_geometry_used=False)))
    destination.write_text(''.join(json.dumps(s)+'\n' for s in starts),encoding='utf-8')
    report=ROOT/'reports/bimolecular_v5';report.mkdir(parents=True,exist_ok=True)
    (report/'coverage.json').write_text(json.dumps(dict(library_audit=dict(library.audit),
        templates=len(library.templates),systems=coverage),indent=2),encoding='utf-8')
    # Held out from the search inputs. Water-solvated free energies are not
    # comparable to gas-phase electronic barriers from AIMNet2-rxn.
    (report/'coley_held_out_reference.json').write_text(json.dumps(dict(selection_rule=
        'First two distinct neutral closed-shell CHNO reactant pairs sorted by explicit atom count then numeric reaction ID',
        eligible_rows=len(eligible),records=selected,
        used_for_proposals=False,used_for_seed_geometry=False),indent=2),encoding='utf-8')
    print(json.dumps(coverage,indent=2))


if __name__=='__main__':main()
