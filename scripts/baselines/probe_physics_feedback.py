#!/usr/bin/env python3
"""Small new-DFT residual update probe; NOT a multi-round active-learning claim."""
import json
from pathlib import Path
import numpy as np
from ase.io import read
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import mean_absolute_error
import joblib

ROOT=Path(__file__).resolve().parents[2]


def main():
    source=np.load(ROOT/"reports/graph_bridge_benchmark/heldout_predictions.npz")
    features=[]; predicted=[]; observed=[]; changes=[]; groups=[]; ids=[]
    for path in sorted((ROOT/"reports/feasibility").glob("*/verification.json")):
        record=json.loads(path.read_text())
        if not (record.get("physical_event_verified") and record.get("expected_endpoint_match")): continue
        name=record["event_id"]
        mask=source["event_ids"]==name
        atoms=read(path.parent/"ts.xyz")
        for (i,j),old,edit,group in zip(source["atom_pairs"][mask],source["forward_R_graph_edits"][mask],
                                       source["target_edits"][mask],source["formula_groups"][mask]):
            z=sorted([int(atoms.numbers[i]),int(atoms.numbers[j])])
            features.append([old,edit,z[0],z[1]])
            predicted.append(old)
            observed.append(float(np.linalg.norm(atoms.positions[i]-atoms.positions[j])))
            changes.append(edit!=0); groups.append(group); ids.append(name)
    folder=ROOT/"reports/physics_feedback_probe"; folder.mkdir(exist_ok=True)
    if len(set(groups))<3:
        result={"status":"insufficient_independent_formula_groups","groups":len(set(groups)),
                "claim":"No feedback benefit can be assessed from fewer than three groups in this probe"}
    else:
        x=np.asarray(features); base=np.asarray(predicted); y=np.asarray(observed)
        groups=np.asarray(groups); ids=np.asarray(ids); changed=np.asarray(changes)
        corrected=np.zeros(len(y))
        for train,test in LeaveOneGroupOut().split(x,y,groups):
            update=make_pipeline(StandardScaler(),Ridge(alpha=10.))
            update.fit(x[train],y[train]-base[train])
            corrected[test]=base[test]+update.predict(x[test])
        events=[]
        for event in np.unique(ids):
            mask=(ids==event)&changed
            events.append({"event_id":event,"changed_pair_MAE_before_A":mean_absolute_error(y[mask],base[mask]),
                           "changed_pair_MAE_after_A":mean_absolute_error(y[mask],corrected[mask])})
        result={"status":"completed","source_model_label_level":"RGD1 B3LYP-D3/TZVP geometry",
                "new_evidence_level":"new restricted wb97x/6-31g(d) verified TS geometry",
                "task":"residual calibration to new quantum evidence, not relabeling source data",
                "groups":len(set(groups)),"events":events,
                "event_mean_MAE_before_A":float(np.mean([e["changed_pair_MAE_before_A"] for e in events])),
                "event_mean_MAE_after_A":float(np.mean([e["changed_pair_MAE_after_A"] for e in events])),
                "candidate_promoted":False,
                "model_artifact_status":"experimental_residual_candidate; original graph models unchanged",
                "validation":"Leave-one-formula-out new-QC updates; base predictions also excluded each target formula",
                "limitations":["Very small nonrandom pilot; no statistical/general active-learning conclusion",
                               "One residual-update experiment, not repeated acquisition rounds",
                               "No inverse curved-arrow model update tested",
                               "Distance accuracy does not certify generated structures or TS connectivity"]}
        result["validation_improved"] = result["event_mean_MAE_after_A"] < result["event_mean_MAE_before_A"]
        final=make_pipeline(StandardScaler(),Ridge(alpha=10.))
        final.fit(x,y-base)
        joblib.dump(final,folder/"qc_residual_update.joblib")
        np.savez_compressed(folder/"heldout_corrections.npz",event_ids=ids,formula_groups=groups,
                            before_A=base,after_A=corrected,target_A=y,changed_pairs=changed)
    (folder/"result.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()
