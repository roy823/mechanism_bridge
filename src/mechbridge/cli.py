import argparse
import json
from pathlib import Path
from . import adapters
from .io import read_jsonl, write_jsonl

def main():
    p = argparse.ArgumentParser(prog="mechbridge")
    s = p.add_subparsers(dest="command", required=True)
    ingest = s.add_parser("ingest", help="Import source records without upgrading verification claims")
    ingest.add_argument("dataset", choices=["mech-csv", "flower", "rgd1", "transition1x"])
    ingest.add_argument("input", type=Path); ingest.add_argument("output", type=Path)
    ingest.add_argument("--limit", type=int); ingest.add_argument("--reaction-column", default="reaction")
    ingest.add_argument("--arrows-column"); ingest.add_argument("--id-column")
    ingest.add_argument("--mapped-csv", type=Path)
    pair = s.add_parser("pair"); pair.add_argument("symbolic"); pair.add_argument("physical"); pair.add_argument("output")
    split = s.add_parser("split"); split.add_argument("input"); split.add_argument("output"); split.add_argument("--seed", type=int, default=17)
    extract = s.add_parser("export-xyz"); extract.add_argument("input"); extract.add_argument("event_id"); extract.add_argument("outdir", type=Path)
    ibo = s.add_parser("ibo", help="Compute IBO features; no arrows automatically invented")
    ibo.add_argument("xyz"); ibo.add_argument("outdir", type=Path)
    ibo.add_argument("--ordered-path", action="store_true")
    verify = s.add_parser("verify", help="Hessian and two-sided downhill verification (not IRC)")
    verify.add_argument("xyz"); verify.add_argument("outdir", type=Path)
    verify.add_argument("--dimer", action="store_true")
    verify.add_argument("--fmax", type=float, default=0.02)
    verify.add_argument("--steps", type=int, default=300)
    verify.add_argument("--calculator", help="Custom ASE calculator factory as module:function")
    verify.add_argument("--calculator-config", type=Path)
    for cmd in [ibo, verify]:
        cmd.add_argument("--charge", type=int, required=True)
        cmd.add_argument("--multiplicity", type=int, required=True)
        cmd.add_argument("--method", default="wb97x")
        cmd.add_argument("--basis", default="6-31g(d)")
        cmd.add_argument("--threads", type=int, default=1)
    a = p.parse_args()
    if a.command == "ingest":
        if a.limit is not None and a.limit < 1:
            p.error("limit must be positive")
        if a.dataset == "mech-csv":
            records = adapters.mech_csv(a.input, a.reaction_column, a.arrows_column, a.id_column, a.limit)
        elif a.dataset == "rgd1":
            records = adapters.rgd1(a.input, a.limit, a.mapped_csv)
        else:
            records = getattr(adapters, a.dataset)(a.input, a.limit)
        print(json.dumps({"records": write_jsonl(a.output, records), "output": str(a.output)}))
    elif a.command == "pair":
        from .pairing import pair_candidates
        n = write_jsonl(a.output, pair_candidates(read_jsonl(a.symbolic), read_jsonl(a.physical)))
        print(json.dumps({"graph_candidates": n, "physically_verified_pairs": 0}))
    elif a.command == "split":
        from .pairing import split_records
        print(write_jsonl(a.output, split_records(read_jsonl(a.input), a.seed)))
    elif a.command == "export-xyz":
        from ase import Atoms
        from ase.io import write
        record = next((r for r in read_jsonl(a.input) if r["event_id"] == a.event_id), None)
        if record is None:
            raise ValueError("Event not found")
        a.outdir.mkdir(parents=True, exist_ok=True)
        for kind, xyz in record["physical"]["positions_A"].items():
            atoms = Atoms(numbers=record["physical"]["atomic_numbers"], positions=xyz)
            atoms.info.update(event_id=record["event_id"], source_status="unverified")
            write(a.outdir / (kind + ".xyz"), atoms)
    elif a.command == "ibo":
        from ase.io import read
        from .electronic import compute_ibo
        result = compute_ibo(read(a.xyz, index=":"), a.outdir, a.charge, a.multiplicity,
                             a.method, a.basis, a.ordered_path, a.threads)
        print(json.dumps({"frames": len(result["frames"]), "arrow_status": result["arrow_status"]}))
    elif a.command == "verify":
        from ase.io import read, write
        from .backends import PySCFCalculator, load_factory
        from .physics import verify_saddle, dimer_refine
        a.outdir.mkdir(parents=True, exist_ok=True)
        if a.calculator:
            config = json.loads(a.calculator_config.read_text()) if a.calculator_config else {}
            config.update(charge=a.charge, multiplicity=a.multiplicity)
            factory = load_factory(a.calculator, config)
            backend = {"factory": a.calculator, "config": config}
        else:
            factory = lambda: PySCFCalculator(a.charge, a.multiplicity, a.method, a.basis, a.threads)
            backend = {"factory": "PySCFCalculator", "method": a.method, "basis": a.basis}
        atoms = read(a.xyz)
        atoms.calc = factory()
        dimer_converged = dimer_refine(atoms, a.outdir, a.fmax, a.steps) if a.dimer else None
        result = verify_saddle(atoms, factory, a.outdir, a.fmax, a.steps)
        result.update(backend=backend, charge=a.charge, multiplicity=a.multiplicity,
                      dimer_optimizer_converged=dimer_converged, input_file=Path(a.xyz).name)
        (a.outdir / "verification.json").write_text(json.dumps(result, indent=2))
        print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
