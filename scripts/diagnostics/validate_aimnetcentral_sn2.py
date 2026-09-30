"""Run the official AIMNetCentral Cl-/CH3Cl worked example with strict GPA checks."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from ase import Atoms
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnetcentral_backend import AIMNetCentralPotential
from mechbridge.event_graph import geometry_mol,bond_orders
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,search_connection,atomic_json
from mechbridge.search_seeds import internal_direction


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--family',choices=['aimnet2','aimnet2-2025','aimnet2-nse'],required=True)
    p.add_argument('--outdir',type=Path,required=True);a=p.parse_args()
    if a.outdir.exists():raise FileExistsError(a.outdir)
    a.outdir.mkdir(parents=True)
    atoms=Atoms(numbers=[6,17,17,1,1,1],positions=[
        [0.,0.,0.],[2.35,0.,0.],[-2.35,0.,0.],
        [0.,1.02,.15],[0.,-.51,.885],[0.,-.51,-.885]])
    direction=np.zeros((6,3));direction[1,0]=1.;direction[2,0]=1.
    direction=internal_direction(atoms.positions,direction)
    backend=AIMNetCentralPotential(a.family,ROOT/'models/aimnetcentral',charge=-1,
        multiplicity=1,device='cpu',threads=2)
    protocol=SearchProtocol(seed_policy='AIMNetCentral_documented_symmetric_SN2_TS_guess_v8',
        total_evaluations=3200,evaluations_per_attempt=3200,ts_steps=300,descent_steps=400,
        fmax=.005,hessian_batch_size=32,mode_displacement=.15)
    calculator=CountedCalculator(backend,protocol.total_evaluations)
    result=search_connection(atoms,direction,calculator,a.outdir/'search',protocol,charge=-1)
    bonded=[]
    for endpoint in result.get('endpoints',[]):
        mol=geometry_mol(atoms.numbers,endpoint['positions_A'],-1)
        bonded.append([i for i in (1,2) if tuple(sorted((0,i))) in bond_orders(mol)])
    passed=(result['status']=='validated_descents' and len(bonded)==2 and
            {tuple(v) for v in bonded}=={(1,),(2,)})
    report=dict(family=a.family,status=result['status'],evaluations=result['evaluations'],
        seconds=result['seconds'],strict_self_exchange_passed=passed,bonded_chlorine_indices=bonded,
        endpoint_graphs=[e['graph_smiles'] for e in result.get('endpoints',[])],
        barriers_eV=[e['barrier_eV'] for e in result.get('endpoints',[])],
        source='https://isayevlab.github.io/aimnetcentral/advanced/reaction_paths/',
        source_input='Documented symmetric SN2 TS guess and reaction direction; reference-assisted, not discovery',
        reference_product_geometry_used=False)
    atomic_json(a.outdir/'result.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
