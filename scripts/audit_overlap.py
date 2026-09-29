#!/usr/bin/env python3
"""Audit original FlowER splits against RGD1 mapped CSV; hits are NOT verified pairs."""
import csv
from collections import Counter, defaultdict
from functools import lru_cache
import io
import json
from pathlib import Path
import re
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mechbridge.chemistry import audit_reaction, canonical

ATOM = re.compile(r"\[(?:\d+)?([A-Z][a-z]?|[bcnops])[^\]]*\]")


def domain(smiles):
    # Only apply the fast filter to completely bracketed atom-mapped SMILES.
    # Unexpected notation is parsed normally instead of silently excluded.
    if re.search(r"[A-Za-z]", ATOM.sub("", smiles)):
        return None
    atoms = ATOM.findall(smiles)
    return frozenset(a.upper() for a in atoms), sum(a.upper() != "H" for a in atoms)


@lru_cache(maxsize=50000)
def graph(smiles):
    return canonical(smiles)


def pair(r, p):
    return tuple(sorted((graph(r), graph(p))))


def main():
    index = defaultdict(list)
    max_heavy = 0
    allowed = set()
    stats = Counter()
    with (ROOT / "data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv").open(
            encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            stats["rgd1_csv_rows"] += 1
            for s in (row["reactant"], row["product"]):
                d = domain(s)
                if d is None:
                    raise ValueError("RGD1 CSV domain filter assumptions require review")
                allowed.update(d[0])
                max_heavy = max(max_heavy, d[1])
            try:
                key = pair(row["reactant"], row["product"])
            except ValueError:
                stats["rgd1_parse_errors"] += 1
                continue
            if key[0] != key[1]:
                index[key].append(row["reaction"])
            if stats["rgd1_csv_rows"] % 25000 == 0:
                print(json.dumps(dict(stats)), flush=True)
    stats["rgd1_distinct_nonidentity_graph_pairs"] = len(index)
    graph.cache_clear()
    seen = set()
    hits_path = ROOT / "reports/full_overlap_candidates.jsonl"
    with zipfile.ZipFile(ROOT / "data/raw/flower/data.zip") as z, hits_path.open(
            "w", encoding="utf-8") as output:
        for split in ("train", "val", "test"):
            member = f"data/flower_dataset/{split}.txt"
            with z.open(member) as source:
                for line_no, line in enumerate(io.TextIOWrapper(source, encoding="utf-8"), 1):
                    stats["flower_rows"] += 1
                    stats[f"flower_{split}_rows"] += 1
                    reaction, sequence = line.strip().strip("`").rsplit("|", 1)
                    r, p = reaction.split(">>")
                    if r == p:
                        stats["identical_smiles_skipped"] += 1
                        continue
                    d = domain(r)
                    if d is not None and (not d[0].issubset(allowed) or d[1] > max_heavy):
                        stats["outside_rgd1_element_or_size_scope"] += 1
                        continue
                    stats["flower_rows_canonicalized"] += 1
                    try:
                        key = pair(r, p)
                    except ValueError:
                        stats["flower_parse_errors"] += 1
                        continue
                    ids = index.get(key)
                    if not ids:
                        continue
                    audit = audit_reaction(r, p)
                    if not (audit["composition_conserved"] and audit["charge_conserved"]):
                        stats["unbalanced_hits_excluded"] += 1
                        continue
                    stats["candidate_flower_rows"] += 1
                    stats["candidate_rgd1_associations"] += len(ids)
                    seen.add(key)
                    output.write(json.dumps({"member": member, "line": line_no,
                                             "sequence": sequence, "reaction": reaction,
                                             "rgd1_ids": ids, "canonical_pair": key,
                                             "verified_pair": False,
                                             "status": "csv_graph_candidate_only"}) + "\n")
            print(json.dumps(dict(stats)), flush=True)
    report = {"variant": "flower_dataset", "splits": ["train", "val", "test"],
              "rgd1_reference": "Zenodo mapped CSV; HDF5 identity not checked here",
              "domain_filter": {"elements": sorted(allowed), "maximum_heavy_atoms": max_heavy},
              **dict(stats), "distinct_candidate_graph_pairs": len(seen),
              "verified_pairs": 0,
              "limitations": ["CSV graph matches are only retrieval candidates",
                              "No geometry atom mapping, spin/environment, IRC or arrow verification",
                              "Old/new FlowER variants are not merged"]}
    (ROOT / "reports/full_overlap_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
