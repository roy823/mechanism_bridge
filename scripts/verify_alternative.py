#!/usr/bin/env python3
"""Investigate an unintended model-generated saddle as a separate candidate event."""
import argparse
import json
from pathlib import Path
import sys
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.verification import verify_event


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent-event",required=True)
    a=p.parse_args()
    source=ROOT/"reports/seed_recovery"/a.parent_event/"heldout_graph_model"
    trial=json.loads((source/"result.json").read_text())
    if trial.get("recovered_reference_saddle") or trial.get("stationary",{}).get("imaginary_count")!=1:
        raise ValueError("Requires an unintended first-order saddle")
    rec=next(r for r in map(json.loads,(ROOT/"data/processed/feasibility_events.jsonl").read_text().splitlines()) if r["event_id"]==a.parent_event)
    rec["event_id"] += "_model_alternative"
    rec["positions_A"]["ts"] = read(source/"final.xyz").positions.tolist()
    rec["candidate_origin"] = "heldout_graph_model_then_DFT_search"
    rec["parent_event_id"] = a.parent_event
    rec["requested_endpoints_are_only_a_hypothesis"] = True
    folder=ROOT/"reports/feasibility"/rec["event_id"]
    if (folder/"verification.json").exists(): raise FileExistsError(folder)
    result=verify_event(rec,folder,threads=2)
    print(json.dumps({k:result.get(k) for k in ["event_id","status","event_outcome","physical_event_verified","expected_endpoint_match","observed_reactant_smiles","observed_product_smiles","error"]}),flush=True)


if __name__=="__main__": main()
