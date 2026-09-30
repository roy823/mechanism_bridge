"""Separate tighter-relaxation diagnosis; never modifies primary network edges."""
import json
import argparse
from pathlib import Path
import sys
import time
import hashlib
import zipfile
from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from rdkit import RDLogger
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.physics import analyze_stationary
from mechbridge.event_graph import geometry_mol,graph_smiles
from mechbridge.reaction_network import CountedCalculator,molecular_rmsd,atomic_json
import numpy as np


def main():
    RDLogger.DisableLog('rdApp.*')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    root=ROOT/'reports/bimolecular_v5';out=args.outdir
    if out.exists():raise FileExistsError(out)
    out.mkdir()
    summary=json.loads((root/'summary.json').read_text())
    selected=[];seen=set()
    for e in summary['events']:
        if (e['classification']=='intermolecular_heavy_atom_bond' and e['strategy']=='arrows'
            and any(e['input_pair_endpoint_match']) and e['system'] not in seen):
            selected.append(e);seen.add(e['system'])
    sources=list((ROOT/'src/mechbridge').glob('*.py'))+[Path(__file__)]
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for p in sources:archive.write(p,p.relative_to(ROOT).as_posix())
    atomic_json(out/'manifest.json',dict(selection='First arrows heavy-atom event with input-pair endpoint per system; posthoc basin diagnosis',
        cases=selected,fmax_eV_A=.001,steps=500,optimizer='ASE BFGS',physical_nodes_merged=False,
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))
    backend=ReactionPotential(ROOT/'models/aimnet2-rxn',threads=1)
    results=[]
    for event in selected:
        n=json.loads((ROOT/event['network_file']).read_text());numbers=n['start']['atomic_numbers']
        endpoint=event['nodes'][event['input_pair_endpoint_match'].index(True)]
        folder=out/event['system'];folder.mkdir()
        calc=CountedCalculator(backend,10000);relaxed=[];started=time.perf_counter()
        for label,index in [('initial',0),('observed',endpoint)]:
            atoms=Atoms(numbers=numbers,positions=n['nodes'][index]['positions_A']);atoms.calc=calc
            with BFGS(atoms,maxstep=.1,logfile=str(folder/(label+'.log')),
                               trajectory=str(folder/(label+'.traj'))) as opt:
                opt.run(fmax=.001,steps=500)
            info=analyze_stationary(atoms,.001,batch_forces=calc.batch_forces,batch_size=32);info.pop('modes')
            info.update(positions_A=atoms.positions.tolist(),graph_smiles=graph_smiles(geometry_mol(numbers,atoms.positions,0)))
            relaxed.append(info);write(folder/(label+'.xyz'),atoms,write_results=False)
        mols=[geometry_mol(numbers,r['positions_A'],0) for r in relaxed]
        rmsd=molecular_rmsd(mols[0],np.array(relaxed[0]['positions_A']),mols[1],np.array(relaxed[1]['positions_A']))
        result=dict(system=event['system'],source=event,refined=relaxed,
            rmsd_A=rmsd if np.isfinite(rmsd) else None,
            energy_difference_eV=abs(relaxed[0]['energy_eV']-relaxed[1]['energy_eV']),
            same_refined_minimum_by_registry_criteria=bool(rmsd<.15 and abs(relaxed[0]['energy_eV']-relaxed[1]['energy_eV'])<.03
                and all(r['force_converged'] and r['imaginary_count']==0 for r in relaxed)),
            evaluations=calc.calls,seconds=time.perf_counter()-started,primary_network_modified=False)
        atomic_json(folder/'result.json',result);results.append(result)
        atomic_json(out/'summary.json',results)
        print(json.dumps({k:result[k] for k in ['system','rmsd_A','energy_difference_eV','same_refined_minimum_by_registry_criteria','evaluations']}),flush=True)


if __name__=='__main__':main()
