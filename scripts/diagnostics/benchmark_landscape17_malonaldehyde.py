"""Reproduce the complete Landscape17 malonaldehyde KTN with AIMNet2-2025.

Reference TS coordinates initialize a separate upper-bound diagnostic.  They
are never mixed with autonomous discovery results.
"""
import hashlib
import json
import re
from pathlib import Path
import sys

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol,graph_smiles
from mechbridge.potentials import load_potential
from mechbridge.reaction_network import (CountedCalculator,SearchProtocol,aligned_rmsd,
    inspect_point,molecular_rmsd,search_sella_connection,atomic_json)


DATA=ROOT/'data/raw/landscape17/extracted/Landscape17/malonaldehyde'
OUT=ROOT/'reports/aimnet2025_transitionnet_v9/landscape17_malonaldehyde_pi'


def energy_from_xyz(path):
    line=path.read_text(encoding='utf-8').splitlines()[1]
    return float(re.search(r'energy=([-+0-9.eE]+)',line).group(1))


def closest_minimum(numbers,positions,minima):
    mol=geometry_mol(numbers,positions,0)
    values=[]
    for row in minima:
        matches=row['mol'].GetSubstructMatches(mol,uniquify=False,useChirality=False,maxMatches=512)
        values.append(min((min(aligned_rmsd(row['positions'][list(match)],positions),
                               aligned_rmsd(-row['positions'][list(match)],positions))
                           for match in matches),default=float('inf')))
    index=int(np.argmin(values))
    return index,float(values[index])


def main():
    if OUT.exists():raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    backend,potential=load_potential('aimnet2-2025',ROOT,'cpu',False)
    protocol=SearchProtocol(total_evaluations=14000,evaluations_per_attempt=14000,
        hessian_batch_size=32,fmax=.005,descent_steps=300,sella_steps=400)
    minima=[]
    for index,path in enumerate(sorted((DATA/'Conformations').glob('coordsmin_*.xyz'))):
        reference=read(path);numbers=reference.numbers.tolist()
        calc=CountedCalculator(backend,5000);calc.attempt_limit=5000
        atoms=reference.copy();atoms.calc=calc
        with BFGS(atoms,maxstep=.1,logfile=str(OUT/f'minimum_{index}.log'),
                  trajectory=str(OUT/f'minimum_{index}.traj')) as opt:
            converged=bool(opt.run(fmax=.003,steps=400))
        stationary,_=inspect_point(atoms,protocol)
        mol=geometry_mol(numbers,atoms.positions,0)
        row=dict(id=index,reference_file=path.relative_to(ROOT).as_posix(),
            reference_energy_eV=energy_from_xyz(path),optimizer_converged=converged,
            reference_to_model_rmsd_A=aligned_rmsd(reference.positions,atoms.positions),
            evaluations=calc.calls,graph_smiles=graph_smiles(mol),stationary=stationary,
            positions=atoms.positions.copy(),mol=mol)
        minima.append(row);write(OUT/f'minimum_{index}.xyz',atoms,write_results=False)
    ts_rows=[]
    paths=sorted((DATA/'Conformations').glob('coordsts_*.xyz'))
    for index,path in enumerate(paths):
        match=re.match(r'coordsts_(\d+)_(\d+)(?:_\d+)?\.xyz',path.name)
        pair=[int(match.group(1)),int(match.group(2))]
        reference=read(path);numbers=reference.numbers.tolist()
        calc=CountedCalculator(backend,protocol.total_evaluations);calc.attempt_limit=protocol.total_evaluations
        probe=reference.copy();probe.calc=calc
        initial,modes=inspect_point(probe,protocol)
        destination=OUT/f'ts_{index:02d}_{pair[0]}_{pair[1]}'
        result=search_sella_connection(probe,modes[0],calc,destination,protocol,0)
        assigned=[]
        for endpoint in result.get('endpoints',[]):
            assigned.append(closest_minimum(numbers,np.asarray(endpoint['positions_A']),minima))
        observed=[a[0] for a in assigned if a[1] < .3]
        reference_pair_recovered=(result['status']=='validated_descents' and
                                  sorted(observed)==sorted(pair))
        ts_rows.append(dict(id=index,reference_file=path.relative_to(ROOT).as_posix(),
            reference_pair=pair,reference_energy_eV=energy_from_xyz(path),
            reference_barriers_eV=[energy_from_xyz(path)-minima[i]['reference_energy_eV'] for i in pair],
            model_at_reference=initial,model_refinement_status=result['status'],
            model_ts_rmsd_A=(aligned_rmsd(reference.positions,np.asarray(result['ts_positions_A']))
                             if 'ts_positions_A' in result else None),
            model_endpoint_assignments=[dict(minimum=i,rmsd_A=r) for i,r in assigned],
            reference_pair_recovered=reference_pair_recovered,evaluations=calc.calls,
            result_file=(destination/'result.json').relative_to(ROOT).as_posix()))
    clean=[]
    for row in minima:
        clean.append({k:v for k,v in row.items() if k not in ('positions','mol')})
    archive=ROOT/'data/raw/landscape17/Landscape17.zip'
    summary=dict(scope='reference_TS_assisted_complete_KTN_reproduction_upper_bound',
        autonomous_discovery=False,reference_product_or_TS_used=True,
        reference=dict(dataset='Landscape17',doi='10.6084/m9.figshare.29949230.v1',
            method='omegaB97x/6-31G(d)',license='CC BY 4.0',minima=len(minima),
            transition_states=len(ts_rows),archive_md5=hashlib.md5(archive.read_bytes()).hexdigest()),
        potential=potential,protocol=protocol.__dict__,minima=clean,transition_states=ts_rows,
        recovered_transition_states=sum(r['reference_pair_recovered'] for r in ts_rows),
        complete_network_recovered=all(r['reference_pair_recovered'] for r in ts_rows))
    atomic_json(OUT/'summary.json',summary)
    start=dict(id='landscape17_malonaldehyde_min0',atomic_numbers=numbers,
        positions_A=minima[0]['positions'].tolist(),charge=0,multiplicity=1,
        provenance=dict(dataset='Landscape17',reference_minimum=0,
            reference_TS_or_product_geometry_used=False,
            note='A reference minimum initializes the separate autonomous search; no other KTN data enter'))
    (OUT/'autonomous_start.jsonl').write_text(json.dumps(start)+'\n',encoding='utf-8')
    print(json.dumps(dict(recovered=summary['recovered_transition_states'],total=len(ts_rows),
        complete=summary['complete_network_recovered']),indent=2))


if __name__=='__main__':main()
