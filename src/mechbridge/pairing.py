"""Exact whole-system graph retrieval. Hits are candidates, never verified pairs."""
import hashlib
from .chemistry import audit_reaction, canonical

def endpoint_smiles(record):
    s = record["symbolic"]
    return s.get("reactant_smiles"), s.get("product_smiles")

def pair_candidates(symbolic, physical):
    index = {}
    for p in physical:
        r, q = endpoint_smiles(p)
        if not r or not q:
            continue
        try:
            audit = audit_reaction(r, q)
        except ValueError:
            continue
        if audit["composition_conserved"] and audit["charge_conserved"]:
            index.setdefault(audit["graph_pair_key"], []).append(p)
    for s in symbolic:
        # Overall reactions / ungrouped arrow sequences are not elementary event labels.
        if s["kind"] != "symbolic_elementary_step":
            continue
        r, q = endpoint_smiles(s)
        try:
            audit = audit_reaction(r, q)
        except ValueError:
            continue
        if not (audit["composition_conserved"] and audit["charge_conserved"]):
            continue
        for p in index.get(audit["graph_pair_key"], []):
            pr, _ = endpoint_smiles(p)
            yield {"symbolic_id": s["event_id"], "physical_id": p["event_id"],
                   "graph_pair_key": audit["graph_pair_key"],
                   "direction": "forward" if canonical(r) == canonical(pr) else "reverse",
                   "status": "whole_system_graph_candidate",
                   "verified_pair": False, "atom_correspondence": None,
                   "remaining_checks": ["joint_atom_mapping", "charge_spin_environment", "geometry_atom_order",
                                        "same_level_saddle_connectivity", "arrow_evidence"]}

def split_records(records, seed=17):
    """Union by same pair AND mechanism sequence, then deterministic component split.

    Prevents reverse/conformer/same-path leakage, but not all scaffold-family leakage.
    """
    records = list(records)
    parent = list(range(len(records)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    seen = {}
    for i, rec in enumerate(records):
        keys = ["record:" + rec["event_id"]]
        r, p = endpoint_smiles(rec)
        if r and p:
            try:
                keys.append("graph:" + audit_reaction(r, p)["graph_pair_key"])
            except ValueError:
                pass
        sequence = rec["symbolic"].get("sequence_id")
        if sequence is not None:
            keys.append("seq:" + rec["source"] + ":" + str(sequence))
        for key in keys:
            if key in seen:
                parent[find(i)] = find(seen[key])
            else:
                seen[key] = i
    groups = {}
    for i, rec in enumerate(records):
        groups.setdefault(find(i), []).append(rec["event_id"])
    for i, rec in enumerate(records):
        group = min(groups[find(i)])
        score = int(hashlib.sha256(f"{seed}:{group}".encode()).hexdigest()[:8], 16) / 2**32
        yield {"event_id": rec["event_id"], "component_id": group,
               "split": "train" if score < .8 else "validation" if score < .9 else "test"}
