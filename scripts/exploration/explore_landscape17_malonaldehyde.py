"""ReactionAtlas-style PES loop on the Landscape17 malonaldehyde benchmark."""
from collections import Counter
from itertools import permutations,product
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.io import read,write
from ase.mep import NEB
from ase.optimize import BFGS,FIRE
from rdkit import Chem
from rdkit.Chem import Lipinski,rdMolTransforms

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol,graph_smiles
from mechbridge.physics import finite_hessian
from mechbridge.potentials import load_potential
from mechbridge.reaction_network import (CountedCalculator,SearchProtocol,aligned_rmsd,
    inspect_point,search_sella_connection,atomic_json)

START=ROOT/'reports/aimnet2025_transitionnet_v9/landscape17_malonaldehyde_pi/autonomous_start.jsonl'
REFERENCE=ROOT/'data/raw/landscape17/extracted/Landscape17/malonaldehyde/Conformations'
OUT=ROOT/'reports/aimnet2025_transitionnet_v9/landscape17_malonaldehyde_pes_loop'


def pi_rmsd(numbers,x,y):
    numbers=np.asarray(numbers);x=np.asarray(x);y=np.asarray(y)
    groups=[np.flatnonzero(numbers==z) for z in sorted(set(numbers))]
    best=float('inf')
    for choices in product(*(permutations(group.tolist()) for group in groups)):
        order=np.empty(len(numbers),dtype=int)
        for group,choice in zip(groups,choices):order[group]=choice
        best=min(best,aligned_rmsd(x,y[order]),aligned_rmsd(-x,y[order]))
    return best


def torsions(mol):
    result=[]
    for i,j in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
        left=[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx()!=j]
        right=[a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx()!=i]
        if left and right:
            left.sort(key=lambda k:(mol.GetAtomWithIdx(k).GetAtomicNum()==1,k))
            right.sort(key=lambda k:(mol.GetAtomWithIdx(k).GetAtomicNum()==1,k))
            result.append((left[0],i,j,right[0]))
    return result


def rotate(numbers,positions,mol,dihedrals,rng):
    work=Chem.Mol(mol);work.RemoveAllConformers();conf=Chem.Conformer(len(numbers))
    for i,p in enumerate(positions):conf.SetAtomPosition(i,tuple(map(float,p)))
    work.AddConformer(conf,assignId=True)
    count=int(rng.integers(1,len(dihedrals)+1))
    for index in rng.choice(len(dihedrals),size=count,replace=False):
        rdMolTransforms.SetDihedralDeg(work.GetConformer(),*dihedrals[index],
                                      float(rng.choice([60,120,180,240,300])))
    return np.asarray(work.GetConformer().GetPositions())


def main():
    if OUT.exists():raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    start=json.loads(START.read_text(encoding='utf-8'))
    numbers=np.asarray(start['atomic_numbers']);root_graph=start['provenance'].get('input_smiles','O=CCC=O')
    backend,potential=load_potential('aimnet2-2025',ROOT,'cpu',False)
    protocol=SearchProtocol(total_evaluations=18000,evaluations_per_attempt=18000,
        hessian_batch_size=32,fmax=.005,sella_steps=400,descent_steps=300)
    calculator=CountedCalculator(backend,120000);calculator.attempt_limit=120000
    nodes=[];representatives=[];symmetry_images=[];attempts=[]
    def add_minimum(atoms,stationary,source,trial):
        for node in nodes:
            rms=pi_rmsd(numbers,np.asarray(node['positions_A']),atoms.positions)
            if rms<.3 and abs(node['energy_eV']-stationary['energy_eV'])<.03:
                fixed=aligned_rmsd(np.asarray(node['positions_A']),atoms.positions)
                duplicate=any(s['node']==node['id'] and aligned_rmsd(np.asarray(s['positions_A']),
                              atoms.positions)<.1 for s in symmetry_images)
                if fixed>.3 and not duplicate and len(symmetry_images)<4:
                    symmetry_images.append(dict(node=node['id'],positions_A=atoms.positions.tolist(),
                        fixed_index_rmsd_A=fixed,pi_rmsd_A=rms,source=source,trial=trial))
                return node['id'],False
        node=dict(id=len(nodes),source=source,trial=trial,positions_A=atoms.positions.tolist(),
                  graph_smiles=graph_smiles(geometry_mol(numbers,atoms.positions,0)),**stationary)
        nodes.append(node);representatives.append(atoms.copy());write(OUT/f'minimum_{node["id"]:02d}.xyz',atoms,write_results=False)
        return node['id'],True
    root=Atoms(numbers=numbers,positions=start['positions_A'],calculator=calculator)
    with BFGS(root,maxstep=.1,logfile=str(OUT/'root.log')) as opt:opt.run(fmax=.003,steps=400)
    root_stat,_=inspect_point(root,protocol);root_graph=graph_smiles(geometry_mol(numbers,root.positions,0))
    add_minimum(root,root_stat,None,-1)
    mol=geometry_mol(numbers,root.positions,0);dihedrals=torsions(mol);rng=np.random.default_rng(20261001)
    for trial in range(32):
        source=trial%len(nodes);candidate=representatives[source].copy()
        candidate.positions=rotate(numbers,candidate.positions,geometry_mol(numbers,candidate.positions,0),dihedrals,rng)
        candidate.calc=calculator;before=calculator.calls
        try:
            with BFGS(candidate,maxstep=.1,logfile=str(OUT/f'basin_{trial:03d}.log')) as opt:
                converged=bool(opt.run(fmax=.003,steps=400))
            stationary,_=inspect_point(candidate,protocol)
            graph=graph_smiles(geometry_mol(numbers,candidate.positions,0))
            if converged and stationary['force_converged'] and stationary['imaginary_count']==0 and graph==root_graph:
                node,new=add_minimum(candidate,stationary,source,trial);status='new_minimum' if new else 'known_minimum'
            else:status='rejected_stationary_or_graph'
        except Exception as exc:
            status='calculation_failed';node=None
        attempts.append(dict(trial=trial,source=source,status=status,node=node,
                             evaluations=calculator.calls-before))
    # Connect each minimum to its three nearest fixed-index neighbours.  Symmetry
    # images are retained as candidates for degenerate self-edges.
    candidates=[dict(node=n['id'],positions=np.asarray(n['positions_A'])) for n in nodes]
    candidates += [dict(node=s['node'],positions=np.asarray(s['positions_A'])) for s in symmetry_images]
    pair_jobs=[(i,j) for i in range(len(nodes)) for j in range(i+1,len(nodes))]
    for index in range(len(nodes),len(candidates)):
        representative=candidates[index]['node']
        pair_jobs.append((representative,index))
    connections=[]
    for job,(a,b) in enumerate(pair_jobs):
        left,right=candidates[a],candidates[b]
        calc=CountedCalculator(backend,18000);calc.attempt_limit=18000
        images=[Atoms(numbers=numbers,positions=left['positions'])]
        images += [Atoms(numbers=numbers,positions=left['positions']) for _ in range(13)]
        images += [Atoms(numbers=numbers,positions=right['positions'])]
        for image in images:image.calc=calc
        neb=NEB(images,k=10.,climb=False,allow_shared_calculator=True,
                remove_rotation_and_translation=True,method='improvedtangent')
        neb.interpolate(method='idpp')
        with FIRE(neb,maxstep=.05,dt=.05,logfile=str(OUT/f'pair_{job:03d}_neb.log')) as opt:
            neb_ok=bool(opt.run(fmax=.05,steps=300))
        energies=[float(i.get_potential_energy()) for i in images]
        maxima=[i for i in range(1,len(images)-1) if energies[i]>=energies[i-1] and energies[i]>=energies[i+1]]
        maxima.sort(key=lambda i:energies[i],reverse=True)
        refinements=[];accepted=None
        for rank,index in enumerate(maxima[:4]):
            direction=images[index+1].positions-images[index-1].positions
            direction/=np.linalg.norm(direction)
            result=search_sella_connection(images[index].copy(),direction,calc,
                OUT/f'pair_{job:03d}_refine_{rank:02d}',protocol,0)
            assigned=[]
            for endpoint in result.get('endpoints',[]):
                rms=[pi_rmsd(numbers,np.asarray(n['positions_A']),np.asarray(endpoint['positions_A'])) for n in nodes]
                best=int(np.argmin(rms));assigned.append(best if rms[best]<.3 else None)
            target=(result['status']=='validated_descents' and sorted(assigned)==sorted([left['node'],right['node']]))
            refinements.append(dict(image=index,status=result['status'],assigned_nodes=assigned,
                                    target_pair=target,evaluations=result['evaluations']))
            if target:accepted=result;break
        connections.append(dict(job=job,candidate_indices=[a,b],target_nodes=[left['node'],right['node']],
            neb_converged=neb_ok,neb_energies_eV=energies,refinements=refinements,
            accepted=accepted is not None,evaluations=calc.calls))
    references=[read(p).positions for p in sorted(REFERENCE.glob('coordsmin_*.xyz'))]
    node_matches=[]
    for node in nodes:
        rms=[pi_rmsd(numbers,np.asarray(node['positions_A']),r) for r in references]
        best=int(np.argmin(rms));node_matches.append(dict(node=node['id'],reference_minimum=best,
            rmsd_A=rms[best],matched=rms[best]<.3,physical_graph=node['graph_smiles']==root_graph))
    summary=dict(method='torsion_basin_search_then_NEB_Sella_connectivity',potential=potential,
        protocol=protocol.__dict__,reference_TS_or_product_geometry_used=False,
        reference_minimum_used_as_initial_seed=True,root_graph=root_graph,dihedrals=dihedrals,
        basin_trials=len(attempts),basin_outcomes=dict(Counter(a['status'] for a in attempts)),
        evaluations=calculator.calls,nodes=nodes,symmetry_images=symmetry_images,
        node_reference_matches=node_matches,pair_jobs=len(pair_jobs),connections=connections,
        accepted_connections=sum(c['accepted'] for c in connections))
    atomic_json(OUT/'summary.json',summary)
    print(json.dumps(dict(nodes=len(nodes),symmetry_images=len(symmetry_images),
        matched_reference_minima=sum(n['matched'] for n in node_matches),
        accepted_connections=summary['accepted_connections']),indent=2))


if __name__=='__main__':main()
