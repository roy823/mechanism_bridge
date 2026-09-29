#!/usr/bin/env python3
"""Select real RGD1 small, single-fragment, closed-shell graph-consistent events."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import h5py
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mechbridge.event_graph import geometry_mol, graph_smiles, pair_slots
from mechbridge.chemistry import canonical, graph_pair_key


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-atoms", type=int, default=8)
    p.add_argument("--limit", type=int, default=12)
    a = p.parse_args()
    RDLogger.DisableLog("rdApp.warning")
    source = ROOT / "data/raw/rgd1_zenodo/RGD1_allrxns.h5"
    candidates = []
    counts = Counter()
    with h5py.File(source, "r") as f:
        for key, g in f.items():
            if g["elements"].shape[0] > a.max_atoms:
                continue
            counts["small_records"] += 1
            r, q = g["Rsmiles"][()].decode(), g["Psmiles"][()].decode()
            if "." in r or "." in q:
                counts["multifragment_excluded"] += 1
                continue
            try:
                parsed = [Chem.MolFromSmiles(s) for s in (r, q)]
                if any(m is None for m in parsed):
                    raise ValueError("Source SMILES parse failed")
                if any(sum(x.GetFormalCharge() for x in m.GetAtoms()) or
                       any(x.GetNumRadicalElectrons() for x in m.GetAtoms()) for m in parsed):
                    counts["charge_or_radical_excluded"] += 1
                    continue
                if canonical(r) == canonical(q):
                    counts["same_graph_excluded"] += 1
                    continue
                numbers = g["elements"][()].tolist()
                positions = {out:g[src][()].tolist() for out,src in
                             [("reactant","RG"),("ts","TSG"),("product","PG")]}
                rm = geometry_mol(numbers, positions["reactant"], 0)
                pm = geometry_mol(numbers, positions["product"], 0)
                if graph_smiles(rm) != canonical(r) or graph_smiles(pm) != canonical(q):
                    counts["geometry_source_graph_disagreement"] += 1
                    continue
                if len(pair_slots(rm)) != len(pair_slots(pm)):
                    raise ValueError("Electron count mismatch")
                candidates.append({"event_id":key,"atomic_numbers":numbers,
                    "positions_A":positions,"reactant_smiles":r,"product_smiles":q,
                    "mapped_reactant":Chem.MolToSmiles(rm),"mapped_product":Chem.MolToSmiles(pm),
                    "graph_pair_key":graph_pair_key(r,q),"charge":0,"multiplicity":1,
                    "state_basis":"explicit closed-shell singlet pilot assumption; SCF checks required",
                    "source_doi":"10.5281/zenodo.7860446","source_group":key,
                    "source_sha256":"ed125b4cb1eac9af670a7cae8b9d29c88a9f8be24d0206a027f0f2c035f8f268",
                    "geometry_graph_agreement":True,"verified_pair":False})
            except (ValueError,RuntimeError) as exc:
                counts["graph_perception_failed"] += 1
    seen = set(); selected = []
    for rec in sorted(candidates,key=lambda x:(len(x["atomic_numbers"]),x["event_id"])):
        if rec["graph_pair_key"] in seen:
            continue
        seen.add(rec["graph_pair_key"])
        if len(selected) < a.limit:
            selected.append(rec)
    dest = ROOT / "data/processed/feasibility_events.jsonl"
    dest.write_text("".join(json.dumps(x)+"\n" for x in selected),encoding="utf-8")
    (ROOT / "data/processed/feasibility_pool.jsonl").write_text(
        "".join(json.dumps(x)+"\n" for x in candidates),encoding="utf-8")
    report = {"selection":"smallest atom count then source ID, no QC outcomes used",
              "max_atoms":a.max_atoms,"counts":dict(counts),"eligible":len(candidates),
              "distinct_graph_pairs":len(seen),"selected":[r["event_id"] for r in selected]}
    (ROOT / "reports/feasibility_selection.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
