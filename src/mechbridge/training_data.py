"""Lazy physical training input. Frames remain grouped by source reaction.

This is a data loader, not a trained ML potential or an arrow predictor.
"""
import numpy as np

def transition1x_frames(path, split="train", max_frames=None, chunk_size=128):
    import h5py
    count=0
    with h5py.File(path,"r") as f:
        for formula in f[split]:
            for reaction,g in f[split][formula].items():
                numbers=np.asarray(g['atomic_numbers'][()])
                n=len(g['positions'])
                for start in range(0,n,chunk_size):
                    positions=g['positions'][start:start+chunk_size]
                    energies=g['wB97x_6-31G(d).energy'][start:start+chunk_size]
                    forces=g['wB97x_6-31G(d).forces'][start:start+chunk_size]
                    for j,(x,e,force) in enumerate(zip(positions,energies,forces)):
                        if max_frames is not None and count>=max_frames:return
                        yield {'event_id':f'transition1x:{formula}/{reaction}',
                               'frame_index':start+j,'atomic_numbers':numbers,
                               'positions_A':x,'energy_eV':float(np.asarray(e).reshape(-1)[0]),
                               'forces_eV_A':force,'source_split':split,
                               'path_kind':'NEB_optimization_sample_not_IRC'}
                        count+=1
