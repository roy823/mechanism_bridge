#!/usr/bin/env python3
"""Blind, held-out graph-proxy TS seeds versus endpoint interpolation on DFT.

Reference TS is used only after optimization for recovery assessment. This tests
net-edit guidance; it does not test the extra value of complete curved arrows.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.backends import PySCFCalculator
from mechbridge.verification import stationary


def align(x, reference):
    xc=x-x.mean(axis=0); rc=reference-reference.mean(axis=0)
    u,_,vt=np.linalg.svd(xc.T@rc)
    d=np.eye(3); d[-1,-1]=np.linalg.det(u@vt)
    return xc@(u@d@vt)+reference.mean(axis=0)


def distance_embedding(matrix):
    n=len(matrix); center=np.eye(n)-np.ones((n,n))/n
    gram=-.5*center@matrix**2@center
    values,vectors=np.linalg.eigh(gram)
    return vectors[:,-3:]*np.sqrt(np.maximum(values[-3:],0))


class BudgetExceeded(RuntimeError): pass


class BudgetCalculator(PySCFCalculator):
    def __init__(self,budget,**kwargs):
        super().__init__(**kwargs); self.budget=budget
    def calculate(self,*args,**kwargs):
        if self.evaluation_count>=self.budget:
            raise BudgetExceeded("Predeclared gradient budget exhausted")
        return super().calculate(*args,**kwargs)


def main():
    from sella import Sella
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--event",required=True)
    p.add_argument("--budget",type=int,default=60)
    p.add_argument("--outdir",type=Path,default=ROOT/"reports/seed_recovery")
    a=p.parse_args()
    record=next(r for r in map(json.loads,(ROOT/"data/processed/feasibility_events.jsonl").read_text().splitlines()) if r["event_id"]==a.event)
    reference_dir=ROOT/"reports/feasibility"/a.event
    reference_report=json.loads((reference_dir/"verification.json").read_text())
    if not reference_report["physical_event_verified"]:
        raise ValueError("Independent reference connection must be verified first")
    predictions=np.load(ROOT/"reports/graph_bridge_benchmark/heldout_predictions.npz")
    mask=predictions["event_ids"]==a.event
    if not mask.any(): raise ValueError("Missing held-out model prediction")
    n=len(record["atomic_numbers"]); distances=np.zeros((n,n))
    for (i,j),value in zip(predictions["atom_pairs"][mask],predictions["forward_R_graph_edits"][mask]):
        distances[i,j]=distances[j,i]=value
    r=np.asarray(record["positions_A"]["reactant"])
    product=align(np.asarray(record["positions_A"]["product"]),r)
    seeds={"endpoint_interpolation":(r+product)/2,
           "heldout_graph_model":align(distance_embedding(distances),r)}
    out=a.outdir/a.event; out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name,positions in seeds.items():
        folder=out/name; folder.mkdir(exist_ok=True)
        if (folder/"result.json").exists(): raise FileExistsError(folder)
        atoms=Atoms(numbers=record["atomic_numbers"],positions=positions)
        calc=BudgetCalculator(a.budget,charge=record["charge"],multiplicity=record["multiplicity"],
                              method=reference_report["method"],basis=reference_report["basis"],threads=2)
        atoms.calc=calc; write(folder/"initial.xyz",atoms)
        result={"event_id":a.event,"seed_method":name,"gradient_budget":a.budget,
                "reference_ts_used_for_seed":False,"status":"started"}
        start=time.time()
        try:
            with Sella(atoms,internal=False,order=1,logfile=str(folder/"optimization.log"),
                       trajectory=str(folder/"optimization.traj")) as opt:
                result["optimizer_converged"]=bool(opt.run(fmax=.02,steps=100))
            check,_,_=stationary(atoms,calc)
            result["stationary"]=check
            reference=read(reference_dir/"ts.xyz")
            rms=float(np.sqrt(np.mean(np.sum((align(atoms.positions,reference.positions)-reference.positions)**2,axis=1))))
            result["reference_ts_RMSD_A"]=rms
            result["reference_energy_difference_eV"]=abs(check["energy_eV"]-reference_report["ts"]["energy_eV"])
            result["recovered_reference_saddle"]=(check["force_converged"] and check["imaginary_count"]==1 and
                rms<.05 and result["reference_energy_difference_eV"]<.01)
            result["status"]="completed"
        except Exception as exc:
            result["status"]="budget_exhausted" if isinstance(exc,BudgetExceeded) else "failed"
            result["error"]=f"{type(exc).__name__}: {exc}"
        result["gradient_evaluations"]=calc.evaluation_count
        result["seconds"]=time.time()-start
        write(folder/"final.xyz",atoms,write_results=False)
        (folder/"result.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
        rows.append(result); print(json.dumps(result),flush=True)
    (out/"comparison.json").write_text(json.dumps({"results":rows,
        "limitations":["One-event search diagnostic, not a generalization claim",
        "Model is conditioned on graph edits, not independent curved-arrow labels",
        "Matching a previously verified saddle is a recovery metric, not a fresh IRC verification"]},indent=2),encoding="utf-8")


if __name__=="__main__": main()
