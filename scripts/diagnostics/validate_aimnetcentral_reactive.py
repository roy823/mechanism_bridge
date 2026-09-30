"""Reference-assisted saddle transfer across AIMNetCentral families; not discovery."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from ase.io import read
from rdkit import RDLogger
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnetcentral_backend import AIMNetCentralPotential,FAMILIES
from mechbridge.event_graph import geometry_mol,resonance_equivalent
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,search_connection,atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=FAMILIES,required=True)
    p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
    if a.outdir.exists():raise FileExistsError(a.outdir)
    a.outdir.mkdir(parents=True);RDLogger.DisableLog('rdApp.*')
    backend=AIMNetCentralPotential(a.family,ROOT/'models/aimnetcentral',device='cpu',threads=2)
    protocol=SearchProtocol(seed_policy='DFT_TS_and_negative_mode_oracle_transfer_v8',
        max_attempts=1,total_evaluations=2400,evaluations_per_attempt=2400,ts_steps=220,
        descent_steps=350,fmax=.005,hessian_batch_size=32,mode_displacement=.15)
    cases=[('acetone_tautomerization',ROOT/'reports/local_transfer_v3/qc/acetone'),
           ('formaldehyde_hydration',ROOT/'reports/bimolecular_v5/qc/formaldehyde_water'),
           ('formaldehyde_dimerization',ROOT/'reports/bimolecular_v5/qc/formaldehyde_dimer_refined')]
    results=[];started=time.perf_counter()
    for name,folder in cases:
        q=json.loads((folder/'verification.json').read_text());atoms=read(folder/'ts.xyz')
        data=np.load(folder/'ts_modes.npz');direction=np.asarray(data['modes'][0])
        calculator=CountedCalculator(backend,protocol.total_evaluations)
        result=search_connection(atoms,direction,calculator,a.outdir/name,protocol)
        actual=[]
        for endpoint in result.get('endpoints',[]):
            try:actual.append(geometry_mol(atoms.numbers,endpoint['positions_A'],q['charge']))
            except ValueError:pass
        expected=[geometry_mol(atoms.numbers,read(folder/f'minimum_{d}.xyz').positions,q['charge'])
                  for d in ('forward','reverse')]
        preserved=(len(actual)==2 and any(all(resonance_equivalent(expected[i],actual[j])
            for i,j in enumerate(order)) for order in ((0,1),(1,0))))
        row=dict(case=name,family=a.family,status=result['status'],evaluations=result['evaluations'],
            seconds=result['seconds'],DFT_reference_pair_preserved=preserved,
            endpoint_graphs=[e['graph_smiles'] for e in result.get('endpoints',[])],
            ts=result.get('ts'),barriers_eV=[e['barrier_eV'] for e in result.get('endpoints',[])],
            reference_assisted=True,reference_product_geometry_used=False,
            reference_TS_and_negative_mode_used=True)
        results.append(row);atomic_json(a.outdir/'results.json',results)
        print(json.dumps({k:row[k] for k in ['case','family','status','evaluations','DFT_reference_pair_preserved']}),flush=True)
    atomic_json(a.outdir/'manifest.json',dict(family=a.family,model=backend.describe(),
        protocol=protocol.__dict__,elapsed_seconds=time.perf_counter()-started,
        evidence='DFT TS plus DFT negative mode supplied; model-PES refinement and descents, not autonomous discovery'))


if __name__=='__main__':main()
