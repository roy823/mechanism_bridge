#!/usr/bin/env python3
"""Test whether automated arrows survive a rigid rotation of the same physical path."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.spatial.transform import Rotation
from ase.io import read,write

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.electronic import compute_ibo
from mechbridge.event_graph import geometry_mol,arrow_hypotheses,contract_atom_relays


def signature(arrows):
    return sorted((tuple(a['source']),tuple(a['sink'])) for a in arrows)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--event',required=True)
    p.add_argument('--outdir',type=Path)
    args=p.parse_args(); source=ROOT/'reports/feasibility'/args.event
    report=json.loads((source/'verification.json').read_text())
    baseline=json.loads((source/'ibo/electronic_report.json').read_text())
    frames=read(source/'electronic_frames.xyz',index=':')
    rotation=Rotation.from_euler('xyz',[.37,.61,-.2]).as_matrix()
    for atoms in frames: atoms.positions=atoms.positions@rotation.T+np.array([2.,-1.,.5])
    output=args.outdir if args.outdir is not None else source/'rigid_rotation_check'
    output.mkdir(parents=True,exist_ok=True)
    if (output/'result.json').exists(): raise FileExistsError(output)
    write(output/'frames.xyz',frames)
    electronic=compute_ibo(frames,output/'ibo',report['charge'],report['multiplicity'],
                           report['method'],report['basis'],ordered_path=True,threads=2)
    rm=geometry_mol(frames[0].numbers,frames[0].positions,report['charge'])
    pm=geometry_mol(frames[-1].numbers,frames[-1].positions,report['charge'])
    initial=np.load(output/'ibo/frame_0000.npz')['iao_atom_populations_e']
    final=np.load(output/f'ibo/frame_{len(frames)-1:04d}.npz')['iao_atom_populations_e']
    symbolic=arrow_hypotheses(initial,final,rm,pm)
    normalized_before=contract_atom_relays(report['symbolic']['arrows'])
    normalized_after=contract_atom_relays(symbolic['arrows'])
    delta=max(abs(a['energy_hartree']-b['energy_hartree']) for a,b in zip(baseline['frames'],electronic['frames']))*27.211386245988
    result={'event_id':args.event,'transformation':'same rigid rotation and translation applied to every frame',
            'max_energy_difference_eV':delta,
            'raw_arrow_set_unchanged':signature(report['symbolic']['arrows'])==signature(symbolic['arrows']),
            'baseline_arrow_count':len(report['symbolic']['arrows']),
            'rotated_arrow_hypotheses':symbolic,
            'atom_relay_normalization_before':normalized_before,
            'atom_relay_normalization_after':normalized_after,
            'normalized_flow_unchanged':normalized_before['normalized_flows']==normalized_after['normalized_flows'],
            'caveat':'Finite DFT grids and orbital localization/tracking can contribute to changes; no independent arrow truth is supplied'}
    (output/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__': main()
