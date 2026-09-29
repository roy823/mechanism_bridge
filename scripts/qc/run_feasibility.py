#!/usr/bin/env python3
"""Run real RGD1 event verification in the isolated QC environment."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from mechbridge.verification import verify_event


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limit",type=int,default=3)
    p.add_argument("--start",type=int,default=0)
    p.add_argument("--threads",type=int,default=2)
    p.add_argument("--outdir",type=Path,default=ROOT/"reports/feasibility")
    a=p.parse_args()
    records=[json.loads(l) for l in (ROOT/"data/processed/feasibility_events.jsonl").read_text().splitlines()]
    for rec in records[a.start:a.start+a.limit]:
        out=a.outdir/rec["event_id"]
        if (out/"verification.json").exists():
            raise FileExistsError(f"Choose a new run output: {out}")
        print(json.dumps({"event_id":rec["event_id"],"status":"starting"}),flush=True)
        result=verify_event(rec,out,threads=a.threads)
        print(json.dumps({k:result[k] for k in ["event_id","status","physical_event_verified","elapsed_seconds"]}),flush=True)


if __name__=="__main__": main()
