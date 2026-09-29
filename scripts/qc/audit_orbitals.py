#!/usr/bin/env python3
"""Audit orbital gauge ambiguity using cross-geometry overlaps, not new SCF labels."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from ase.io import read
from pyscf import gto,lib

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.event_graph import geometry_mol,arrow_hypotheses


def signature(result):
    return tuple(sorted((tuple(a["source"]),tuple(a["sink"])) for a in result["arrows"]))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--event",required=True)
    a=p.parse_args(); lib.num_threads(1)
    folder=ROOT/"reports/feasibility"/a.event
    meta=json.loads((folder/"ibo/electronic_report.json").read_text())
    frames=read(folder/"electronic_frames.xyz",index=":")
    records=[]; previous=None
    for i,atoms in enumerate(frames):
        data=np.load(folder/f"ibo/frame_{i:04d}.npz")
        mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
                  charge=meta["charge"],spin=0,basis=meta["basis"],unit="Angstrom",verbose=0)
        coeff=data["ibo_coeff_ao"]
        if previous is not None:
            pmol,pc=previous
            overlap=pc.T@gto.intor_cross("int1e_ovlp",pmol,mol)@coeff
            records.append({"frame":i,"minimum_individual_squared_overlap":float(np.diag(overlap**2).min()),
                            "minimum_occupied_subspace_squared_overlap":float(np.linalg.svd(overlap,compute_uv=False).min()**2)})
        previous=(mol,coeff)
    weights=overlap**2
    rows,cols=linear_sum_assignment(-weights)
    best=float(weights[rows,cols].sum())
    source_pop=np.load(folder/"ibo/frame_0000.npz")["iao_atom_populations_e"]
    target_pop=np.load(folder/f"ibo/frame_{len(frames)-1:04d}.npz")["iao_atom_populations_e"]
    rm=geometry_mol(frames[0].numbers,frames[0].positions,meta["charge"])
    pm=geometry_mol(frames[-1].numbers,frames[-1].positions,meta["charge"])
    candidates=[]; seen=set()
    for forbidden in [None,*zip(rows,cols)]:
        cost=-weights.copy()
        if forbidden is not None: cost[forbidden]=1e6
        rr,cc=linear_sum_assignment(cost)
        gap=best-float(weights[rr,cc].sum())
        if gap>.05: continue
        hypothesis=arrow_hypotheses(source_pop,target_pop[:,cc[np.argsort(rr)]],rm,pm)
        sig=signature(hypothesis)
        if sig in seen: continue
        seen.add(sig)
        candidates.append({"total_squared_overlap_score_gap":gap,**hypothesis})
    result={"event_id":a.event,"adjacent_frames":records,
            "last_step_alternative_hypotheses":candidates,
            "alternatives_are_exhaustive":False,"overlap_score_gap_cutoff":.05,
            "review_required":True,
            "meaning":"Individual IBO identity can be ambiguous even when occupied subspace is stable; candidates are not independent chemical truth"}
    (folder/"orbital_ambiguity.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({"event_id":a.event,"last_step":records[-1],"alternative_arrow_sets":len(candidates)},indent=2))


if __name__=="__main__": main()
