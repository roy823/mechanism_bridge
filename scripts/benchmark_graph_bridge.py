#!/usr/bin/env python3
"""Grouped graph-level baselines; explicitly NOT an arrow-mechanism benchmark.

Forward: reactant features / + net graph edits -> source TS pair distances.
Inverse: reactant features / + source TS geometry -> signed net bond changes.
All conformers and reverse events of a formula stay on the same fold.
"""
import json
from pathlib import Path
import sys
import importlib.metadata
import joblib
from collections import Counter
import numpy as np
from rdkit import Chem
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error, f1_score, precision_recall_fscore_support

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.event_graph import geometry_mol, bond_orders


def main():
    records=[json.loads(l) for l in (ROOT/"data/processed/feasibility_pool.jsonl").read_text().splitlines()]
    xr=[]; xs=[]; xp=[]; ydistance=[]; yedit=[]; groups=[]; ids=[]; pairs=[]
    for record in records:
        z=np.asarray(record["atomic_numbers"])
        pos={k:np.asarray(v) for k,v in record["positions_A"].items()}
        rm=geometry_mol(z,pos["reactant"],record["charge"])
        pm=geometry_mol(z,pos["product"],record["charge"])
        rb,pb=bond_orders(rm),bond_orders(pm)
        topo=Chem.GetDistanceMatrix(rm)
        formula=str(sorted(Counter(z.tolist()).items()))
        for i in range(len(z)):
            for j in range(i+1,len(z)):
                r=rb.get((i,j),0.); p=pb.get((i,j),0.)
                dr=float(np.linalg.norm(pos["reactant"][i]-pos["reactant"][j]))
                dt=float(np.linalg.norm(pos["ts"][i]-pos["ts"][j]))
                atom_features=sorted([(int(z[k]),rm.GetAtomWithIdx(k).GetDegree(),
                                       rm.GetAtomWithIdx(k).GetFormalCharge()) for k in (i,j)])
                base=[*atom_features[0],*atom_features[1],r,dr,float(topo[i,j]),len(z)]
                changes=[pb.get(b,0)-rb.get(b,0) for b in rb.keys()|pb.keys() if i in b or j in b]
                xr.append(base)
                xs.append(base+[p,p-r,sum(max(c,0) for c in changes),sum(max(-c,0) for c in changes)])
                xp.append(base+[dt,dt-dr])
                ydistance.append(dt); yedit.append(int(round(p-r)))
                groups.append(formula); ids.append(record["event_id"]); pairs.append([i,j])
    xr,xs,xp=map(np.asarray,(xr,xs,xp)); yd=np.asarray(ydistance); ye=np.asarray(yedit)
    groups=np.asarray(groups); ids=np.asarray(ids)
    outputs={name:np.zeros(len(yd)) for name in ["forward_R","forward_R_graph_edits","inverse_R","inverse_R_TS"]}
    folds=np.full(len(yd),-1)
    for fold,(train,test) in enumerate(GroupKFold(n_splits=5).split(xr,ye,groups)):
        folds[test]=fold
        for name,x in [("forward_R",xr),("forward_R_graph_edits",xs)]:
            model=RandomForestRegressor(n_estimators=120,max_depth=12,min_samples_leaf=2,
                                        random_state=17,n_jobs=2)
            model.fit(x[train],yd[train]); outputs[name][test]=model.predict(x[test])
        for name,x in [("inverse_R",xr),("inverse_R_TS",xp)]:
            model=RandomForestClassifier(n_estimators=120,max_depth=12,min_samples_leaf=2,
                                         class_weight="balanced_subsample",random_state=17,n_jobs=2)
            model.fit(x[train],ye[train]); outputs[name][test]=model.predict(x[test])
    metrics={}
    for name,pred in outputs.items():
        if name.startswith("forward"):
            metrics[name]={"pair_distance_MAE_A":mean_absolute_error(yd,pred),
                           "changed_pair_distance_MAE_A":mean_absolute_error(yd[ye!=0],pred[ye!=0])}
        else:
            precision,recall,f1,_=precision_recall_fscore_support(ye!=0,pred!=0,average="binary",zero_division=0)
            metrics[name]={"signed_edit_macro_F1":f1_score(ye,pred,average="macro",zero_division=0),
                           "changed_pair_precision":precision,"changed_pair_recall":recall,
                           "changed_pair_F1":f1,
                           "exact_event_edit_fraction":float(np.mean([np.all(ye[ids==k]==pred[ids==k]) for k in np.unique(ids)]))}
    # Event-group bootstrap avoids treating atom pairs as independent reactions.
    unique=np.unique(groups); differences=[]
    for group in unique:
        mask=(groups==group)&(ye!=0)
        if mask.any():
            differences.append(float(np.mean(np.abs(yd[mask]-outputs["forward_R"][mask]))-
                                     np.mean(np.abs(yd[mask]-outputs["forward_R_graph_edits"][mask]))))
    rng=np.random.default_rng(17)
    bootstrap=np.mean(rng.choice(differences,size=(2000,len(differences)),replace=True),axis=1)
    report={"scope":"graph-level feasibility proxy, not a curved-arrow or validated-path benchmark",
            "runtime_versions":{k:importlib.metadata.version(k) for k in ["numpy","scipy","rdkit","scikit-learn"]},
            "events":len(records),"formula_groups":len(unique),"atom_pair_rows":len(yd),
            "split":"5-fold GroupKFold by atomic composition; no reverse/conformer cross-fold leakage",
            "metrics":metrics,"forward_group_mean_MAE_reduction_A":float(np.mean(differences)),
            "forward_group_bootstrap_95_percent_CI_A":np.quantile(bootstrap,[.025,.975]).tolist(),
            "limitations":["Forward symbolic features are net bond edits, not full electron arrows",
                           "R+edits has product graph information unavailable to R-only; this does NOT establish benefit beyond R/P graphs",
                           "Inverse reconstructs net graph edits from R and TS; it is NOT electronic arrow decoding",
                           "TS labels are source B3LYP-D3/TZVP structures; most are not independently reverified here",
                           "No new-physics feedback-learning claim is made"]}
    dest=ROOT/"reports/graph_bridge_benchmark"; dest.mkdir(exist_ok=True)
    (dest/"metrics.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    np.savez_compressed(dest/"heldout_predictions.npz",event_ids=ids,formula_groups=groups,
                        atom_pairs=pairs,fold=folds,target_ts_distances_A=yd,target_edits=ye,**outputs)
    forward=RandomForestRegressor(n_estimators=120,max_depth=12,min_samples_leaf=2,random_state=17,n_jobs=2)
    inverse=RandomForestClassifier(n_estimators=120,max_depth=12,min_samples_leaf=2,
                                   class_weight="balanced_subsample",random_state=17,n_jobs=2)
    forward.fit(xs,yd); inverse.fit(xp,ye)
    joblib.dump({"forward":forward,"inverse":inverse,
                 "scope":"graph proxies only; metrics come from separate held-out predictions"},
                 dest/"graph_proxy_models.joblib",compress=3)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    axes[0].bar(["Reactant", "+ graph edits"],
                [metrics[k]["changed_pair_distance_MAE_A"] for k in ["forward_R","forward_R_graph_edits"]],
                color=["#718096","#2b6cb0"])
    axes[0].set_ylabel("Changed-pair TS distance MAE (Angstrom)")
    axes[0].set_title("Graph -> TS geometry proxy (lower is better)")
    axes[1].bar(["Reactant", "+ TS geometry"],
                [metrics[k]["changed_pair_F1"] for k in ["inverse_R","inverse_R_TS"]],
                color=["#718096","#2b6cb0"])
    axes[1].set_ylim(0,1); axes[1].set_ylabel("Changed-pair F1")
    axes[1].set_title("TS geometry -> graph edits (higher is better)")
    fig.suptitle("5-fold composition-grouped validation; graph proxies, not curved-arrow mechanisms",fontsize=10)
    fig.tight_layout(); fig.savefig(dest/"baseline_results.png",dpi=180); plt.close(fig)
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
