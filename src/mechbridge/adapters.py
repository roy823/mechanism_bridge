"""Streaming source adapters. Preserve claims; never manufacture missing labels."""
import ast
import csv
import itertools
from pathlib import Path
import numpy as np
from .io import event, sha256
from .chemistry import audit_reaction, split_reaction, net_changes

def provenance(path, doi=None):
    return {"file": Path(path).name, "sha256": sha256(path), "doi": doi,
            "importer": "mechanism-bridge/0.1.0"}

def attach_graph(record, reaction):
    r, p, agents = split_reaction(reaction)
    record["symbolic"].update(reactant_smiles=r, product_smiles=p, agents_smiles=agents)
    try:
        a = audit_reaction(r, p)
        record["graph_audit"] = a
        if a["charge_conserved"]:
            record["system"]["charge"] = a["reactant"]["charge"]
        if a["mapped_atom_identity_conserved"]:
            record["symbolic"]["net_bond_changes"] = net_changes(r, p)
    except ValueError as exc:
        record["graph_audit"] = {"parse_error": str(exc)}
    return record

def flow_lines(lines, source, prov, limit=None, sequence_namespace=""):
    for i, line in enumerate(itertools.islice(lines, limit)):
        line = line.strip().strip('`')
        if not line:
            continue
        reaction, sequence = line.rsplit("|", 1)
        source_id = f"{sequence_namespace}{prov.get('zip_member', prov.get('file', 'inline'))}:{i}"
        rec = event(source, source_id, "symbolic_elementary_step", {**prov, "line": i + 1})
        attach_graph(rec, reaction)
        rec["symbolic"].update(sequence_id=sequence_namespace + sequence,
                               arrow_pairs=None, label_origin="upstream_mechanistic_dataset",
                               arrow_grouping="not_provided", identical_endpoint=(reaction.split(">>")[0] == reaction.split(">>")[-1]))
        rec["validation"]["arrows"] = "not_explicit_in_source"
        yield rec

def flower(path, limit=None):
    prov = provenance(path, "10.6084/m9.figshare.28359407.v3")
    with open(path, encoding="utf-8-sig") as f:
        yield from flow_lines(f, "flower", prov, limit, Path(path).parent.name + ":")

def arrow_sites(value):
    def site(x):
        if isinstance(x, int) and not isinstance(x, bool):
            return {"kind": "atom", "atom_maps": [x]}
        if isinstance(x, float) and np.isfinite(x):
            return {"kind": "source_specific_decimal_site", "source_token": str(x),
                    "resolved_atom_map": None}
        if isinstance(x, (list, tuple)) and len(x) == 2:
            return {"kind": "bond", "sites": [site(i) for i in x]}
        raise ValueError(f"Unsupported arrow site: {x!r}")
    if not isinstance(value, (list, tuple)):
        raise ValueError("Expected flat source-sink list")
    result = []
    for pair in value:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError("Malformed arrow pair; record retained for review")
        result.append({"source": site(pair[0]), "sink": site(pair[1]), "electrons": 2})
    return result

def mech_csv(path, reaction_column, arrows_column=None, id_column=None, limit=None):
    prov = provenance(path)
    with open(path, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for col in [reaction_column, arrows_column, id_column]:
            if col and col not in reader.fieldnames:
                raise ValueError(f"Column {col!r} absent. Actual columns: {reader.fieldnames}")
        for i, row in enumerate(itertools.islice(reader, limit)):
            source = "mech_uspto" if arrows_column else "uspto_source"
            rec = event(source, row[id_column] if id_column else i, "symbolic_overall_reaction",
                        {**prov, "row": i + 2})
            attach_graph(rec, row[reaction_column])
            rec["symbolic"]["upstream_row"] = row
            if arrows_column:
                try:
                    arrows = arrow_sites(ast.literal_eval(row[arrows_column]))
                    rec["symbolic"].update(arrow_pairs=arrows, arrow_grouping="unknown",
                                           label_origin="expert_template_generated")
                    rec["validation"]["arrows"] = "source_annotation_unverified"
                except (ValueError, SyntaxError, TypeError) as exc:
                    rec["symbolic"]["arrow_parse_error"] = str(exc)
            yield rec

def _scalar(group, key):
    value = np.asarray(group[key][()])
    if value.size != 1 or not np.isfinite(value).all():
        raise ValueError(f"Expected finite scalar for {key}")
    return float(value.reshape(-1)[0])

def _text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)

def rgd1(path, limit=None, mapped_csv=None):
    import h5py
    from .chemistry import canonical
    prov = provenance(path, "10.6084/m9.figshare.21066901.v6")
    mapped = {}
    if mapped_csv:
        with open(mapped_csv, newline="", encoding="utf-8-sig") as handle:
            mapped = {row["reaction"]: row for row in csv.DictReader(handle)}
        prov["mapped_csv"] = provenance(mapped_csv)
    with h5py.File(path, "r") as f:
        for key in itertools.islice(f.keys(), limit):
            g = f[key]
            rec = event("rgd1", key, "physical_event_candidate", {**prov, "hdf5_group": key})
            attach_graph(rec, _text(g["Rsmiles"][()]) + ">>" + _text(g["Psmiles"][()]))
            if key in mapped:
                row = mapped[key]
                original = dict(rec["symbolic"])
                csv_pair = [canonical(row["reactant"]), canonical(row["product"])]
                h5_pair = [canonical(original["reactant_smiles"]), canonical(original["product_smiles"])]
                forward = csv_pair == h5_pair
                reverse = csv_pair[::-1] == h5_pair
                agrees = forward or reverse
                rec["provenance"]["mapped_csv_endpoint_graphs_agree"] = agrees
                if agrees:
                    rec["provenance"]["mapped_csv_direction"] = "forward" if forward else "reverse"
                    mapped_r, mapped_p = (row["reactant"], row["product"]) if forward else (row["product"], row["reactant"])
                    attach_graph(rec, mapped_r + ">>" + mapped_p)
                    rec["symbolic"]["hdf5_unmapped_endpoints"] = original
                else:
                    rec["symbolic"]["mapped_csv_conflict"] = row
            numbers = np.asarray(g["elements"][()]).astype(int)
            geometries = {}
            for output, field in [("reactant", "RG"), ("ts", "TSG"), ("product", "PG")]:
                coordinates = np.asarray(g[field][()])
                if coordinates.shape != (len(numbers), 3) or not np.isfinite(coordinates).all():
                    raise ValueError(f"Bad {key}/{field} shape or coordinates")
                geometries[output] = coordinates.tolist()
            energies = {k: _scalar(g, k) for k in ["R_E", "P_E", "TS_E", "R_H", "P_H", "TS_H", "R_F", "P_F", "TS_F"] if k in g}
            rec["physical"] = {"atomic_numbers": numbers.tolist(), "positions_A": geometries,
                                "source_energies_hartree": energies,
                                "energy_geometry_correspondence": "R/P energies may be sums of optimized fragments; RG/PG are unoptimized",
                                "method": "B3LYP-D3/TZVP", "environment": "gas_phase",
                                "endpoint_geometry_status": "unoptimized_in_this_file",
                                "geometry_atom_order": "source_order_not_linked_to_smiles_maps",
                                "source_connectivity_claim": "xTB_IRC_and_DFT_classifier_assisted",
                                "irc": None}
            rec["system"]["environment"] = "gas_phase"
            yield rec

def transition1x(path, limit=None, split="data"):
    import h5py
    prov = provenance(path, "10.6084/m9.figshare.19614657.v4")
    count = 0
    with h5py.File(path, "r") as f:
        for formula in f[split]:
            for name, g in f[split][formula].items():
                if limit is not None and count >= limit:
                    return
                rec = event("transition1x", formula + "/" + name, "physical_event_candidate",
                            {**prov, "hdf5_group": g.name, "source_split": split})
                numbers = np.asarray(g["atomic_numbers"][()]).astype(int).reshape(-1)
                physical = {"atomic_numbers": numbers.tolist(), "positions_A": {}, "source_energies_eV": {},
                            "source_forces_eV_A": {},
                            "method": "wB97X/6-31G(d)", "environment": "gas_phase", "irc": None,
                            "trajectory_status": "NEB_optimization_samples_not_ordered_IRC",
                            "source_connectivity_claim": "NEB_endpoint_conditioned_not_independent_IRC"}
                for src, dst in [("reactant", "reactant"), ("transition_state", "ts"), ("product", "product")]:
                    sg = g[src]
                    xyz = np.asarray(sg["positions"][()]).reshape(-1, len(numbers), 3)
                    if len(xyz) != 1 or not np.isfinite(xyz).all():
                        raise ValueError("Expected exactly one finite final structure per endpoint subgroup")
                    physical["positions_A"][dst] = xyz[0].tolist()
                    physical["source_energies_eV"][dst] = _scalar(sg, "wB97x_6-31G(d).energy")
                    if "wB97x_6-31G(d).forces" in sg:
                        force = np.asarray(sg["wB97x_6-31G(d).forces"][()]).reshape(-1,len(numbers),3)
                        if len(force) != 1 or not np.isfinite(force).all():
                            raise ValueError("Invalid endpoint force data")
                        physical["source_forces_eV_A"][dst] = force[0].tolist()
                rec["physical"] = physical
                rec["system"]["environment"] = "gas_phase"
                yield rec
                count += 1
