"""Oracle-seed diagnostic, explicitly separate from reactant-only discovery."""
from dataclasses import asdict
import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
from ase.io import read
from rdkit import RDLogger
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.aimnet_backend import ReactionPotential
from mechbridge.reaction_network import CountedCalculator,SearchProtocol,search_connection,inspect_point,atomic_json
from mechbridge.search_seeds import internal_direction
from mechbridge.symbolic_library import parse_explicit
from mechbridge.event_graph import geometry_mol,resonance_equivalent


def main():
    RDLogger.DisableLog('rdApp.*')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--fmax',type=float,default=.005)
    args=parser.parse_args()
    root=ROOT/'reports/bimolecular_v5';out=args.outdir
    if out.exists():raise FileExistsError(out)
    out.mkdir()
    reference=json.loads((root/'coley_held_out_reference.json').read_text())
    protocol=SearchProtocol(seed_policy='author_reference_TS_and_vibration_oracle',
        total_evaluations=6000,evaluations_per_attempt=6000,
        ts_steps=300,descent_steps=500,fmax=args.fmax,hessian_batch_size=32)
    sources=list((ROOT/'src/mechbridge').glob('*.py'))+[Path(__file__)]
    with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for f in sources:archive.write(f,f.relative_to(ROOT).as_posix())
    atomic_json(out/'manifest.json',dict(reference_TS_used=True,reference_vibration_used=True,
        purpose='Oracle-seed PES/optimizer diagnostic, not reactant-only discovery',protocol=asdict(protocol),
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))
    backend=ReactionPotential(ROOT/'models/aimnet2-rxn',threads=2);results=[]
    for record in reference['records']:
        folder=ROOT/'data/raw/coley_dipolar/selected_profiles'/record['rxn_id']
        paths=[p for p in folder.glob('TS_*.xyz') if p.name!='TS_imag_mode.xyz']
        if len(paths)!=1:raise ValueError('Ambiguous reference TS')
        atoms=read(paths[0]);frames=read(folder/'TS_imag_mode.xyz',index=':')
        if not np.array_equal(atoms.numbers,frames[0].numbers):raise ValueError('Vibration atom order mismatch')
        direction=internal_direction(atoms.positions,frames[1].positions-frames[0].positions)
        calculator=CountedCalculator(backend,6000);atoms.calc=calculator
        initial,_=inspect_point(atoms,protocol)
        result=search_connection(atoms,direction,calculator,out/record['rxn_id'],protocol)
        expected=[parse_explicit(s) for s in record['rxn_smiles'].split('>>')]
        actual=[geometry_mol(atoms.numbers,e['positions_A'],0) for e in result['endpoints']]
        pair_match=(len(actual)==2 and any(all(resonance_equivalent(expected[i],actual[j])
            for i,j in enumerate(order)) for order in ((0,1),(1,0))))
        results.append(dict(rxn_id=record['rxn_id'],reference_TS_used=True,initial_MLIP_stationarity=initial,
            outcome=result['status'],reference_pair_preserved=pair_match,
            endpoint_graphs=[e['graph_smiles'] for e in result['endpoints']],
            total_evaluations=calculator.calls,search_evaluations=result['evaluations'],
            source_TS=paths[0].relative_to(ROOT).as_posix(),
            source_TS_sha256=hashlib.sha256(paths[0].read_bytes()).hexdigest()))
        atomic_json(out/'summary.json',results)
        print(json.dumps(results[-1]),flush=True)


if __name__=='__main__':main()
