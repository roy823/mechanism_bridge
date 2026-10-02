"""Real-data pilot: same-PES TS refinement, Sella IRC and endpoint identity checks."""
import json
from pathlib import Path
import time
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from .backends import PySCFCalculator
from .physics import vibrational_analysis
from .event_graph import geometry_mol, graph_smiles, arrow_hypotheses
from .chemistry import canonical
from .electronic import compute_ibo


def classify_connection(ts, endpoints, expected_smiles):
    """A valid unexpected connection is new evidence, never a failed reaction."""
    actual = [e["graph_smiles"] for e in endpoints]
    matched = sorted(actual) == sorted(expected_smiles)
    distinct = len(actual) == 2 and actual[0] != actual[1]
    verified = (len(endpoints) == 2 and distinct and ts["force_converged"] and
                ts["imaginary_count"] == 1 and all(
                    e["irc_converged"] and e["force_converged"] and
                    e["imaginary_count"] == 0 and e["barrier_electronic_eV"] >= 0
                    for e in endpoints))
    return {"expected_endpoint_match":matched,"distinct_endpoint_graphs":distinct,
            "physical_event_verified":verified,
            "event_outcome":("intended_valid_event" if matched else "alternative_valid_event")
                            if verified else "unresolved"}


def stationary(atoms, calculator, fmax=0.02):
    atoms.calc = calculator
    force = float(np.linalg.norm(atoms.get_forces(),axis=1).max())
    hessian = calculator.hessian(atoms)
    vib = vibrational_analysis(atoms.positions, atoms.get_masses(), hessian)
    modes = vib.pop("modes")
    return {"force_max_eV_A":force,"force_converged":force<=fmax,
            "energy_eV":float(atoms.get_potential_energy()),**vib}, hessian, modes


def verify_event(record, outdir, method="wb97x", basis="6-31g(d)", threads=2,
                 ts_steps=100, irc_steps=120, frames=17, ts_optimizer_fmax=0.01, device="cpu",
                 polish_steps=100):
    """Validate a supplied candidate on DFT; frames=0 omits orbital annotation."""
    if frames == 1 or frames < 0:
        raise ValueError('Use frames=0 for physical verification, or frames>=2 for orbital analysis')
    if not 0 < ts_optimizer_fmax <= .02:
        raise ValueError('TS optimizer target must be positive and no looser than the physical force gate')
    from sella import Sella, IRC
    outdir = Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    started = time.time()
    result = {"event_id":record["event_id"],"source":record,
              "method":method,"basis":basis,"charge":record["charge"],
              "multiplicity":record["multiplicity"],"environment":"gas_phase",
              "protocol":{"threads":threads,"device":device,"ts_max_steps":ts_steps,"irc_max_steps":irc_steps,
                          "endpoint_polish_max_steps":polish_steps,
                          "irc_dx_A_sqrt_amu":0.08,"ts_minimum_fmax_eV_A":0.02,
                          "ts_optimizer_fmax_eV_A":ts_optimizer_fmax,
                          "irc_fmax_eV_A":0.05,"irc_inner_fmax_eV_A":0.02,
                          "scf_initial_guess":"independent_atomic_guess",
                          "electronic_state_stability_test":"not_performed"},
              "physical_event_verified":False,"verified_pair":False,
              "symbolic_status":"not_computed","status":"started"}
    def save():
        result["elapsed_seconds"] = time.time()-started
        result['gradient_evaluations'] = calc.evaluation_count
        temp = outdir / "verification.json.tmp"
        temp.write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
        temp.replace(outdir / "verification.json")
    calc = PySCFCalculator(record["charge"],record["multiplicity"],method,basis,threads,device=device)
    try:
        atoms = Atoms(numbers=record["atomic_numbers"],positions=record["positions_A"]["ts"])
        atoms.calc = calc
        result["status"] = "refining_source_ts"; save()
        with Sella(atoms,order=1,internal=False,logfile=str(outdir/"ts.log"),
                   trajectory=str(outdir/"ts.traj")) as opt:
            result["ts_optimizer_converged"] = bool(opt.run(fmax=ts_optimizer_fmax,steps=ts_steps))
        ts, hessian, modes = stationary(atoms,calc)
        result["ts"] = ts
        np.savez_compressed(outdir/"ts_modes.npz",hessian_eV_A2=hessian,modes=modes)
        write(outdir/"ts.xyz",atoms)
        if not ts["force_converged"] or ts["imaginary_count"] != 1:
            result["status"] = "unresolved_ts"; save(); return result
        result["status"] = "integrating_irc"; save()
        ts_atoms = atoms.copy()
        branches = []; minima = []
        for direction in ("forward","reverse"):
            branch = ts_atoms.copy(); branch.calc = calc
            with IRC(branch,dx=0.08,ninner_iter=20,keep_going=False,
                     logfile=str(outdir/f"irc_{direction}.log"),
                     trajectory=str(outdir/f"irc_{direction}.traj")) as irc:
                converged = bool(irc.run(fmax=0.05,fmax_inner=0.02,
                                         steps=irc_steps,direction=direction))
            if not converged:
                result["status"] = "unresolved_irc"; save(); return result
            path = read(outdir/f"irc_{direction}.traj",index=":")
            with BFGS(branch,logfile=str(outdir/f"minimum_{direction}.log")) as opt:
                opt.run(fmax=0.02,steps=polish_steps)    # historical cap 100
            end, _, _ = stationary(branch,calc)
            end["irc_converged"] = converged
            end["graph_smiles"] = graph_smiles(geometry_mol(branch.numbers,branch.positions,record["charge"]))
            end["barrier_electronic_eV"] = ts["energy_eV"] - end["energy_eV"]
            write(outdir/f"minimum_{direction}.xyz",branch)
            branches.append(path); minima.append((branch.copy(),end))
            result["completed_irc_branches"] = [e for _,e in minima]
            save()
        expected = [canonical(record["reactant_smiles"]),canonical(record["product_smiles"])]
        actual = [end["graph_smiles"] for _,end in minima]
        result["endpoints"] = [end for _,end in minima]
        result.update(classify_connection(ts,result["endpoints"],expected))
        verified = result["physical_event_verified"]
        result["is_IRC"] = True
        if not verified:
            result["status"] = "connection_not_certified"; save(); return result
        if result["expected_endpoint_match"] and actual[0] != expected[0]:
            branches.reverse(); minima.reverse()
        result["observed_reactant_smiles"] = minima[0][1]["graph_smiles"]
        result["observed_product_smiles"] = minima[1][1]["graph_smiles"]
        path = [minima[0][0]] + branches[0][::-1] + [ts_atoms] + branches[1] + [minima[1][0]]
        write(outdir/"ordered_path.xyz",path)
        if frames == 0:
            result['status'] = 'physical_event_verified'
            result['symbolic_status'] = 'not_requested'
            save()
            return result
        indices = np.unique(np.linspace(0,len(path)-1,min(frames,len(path))).astype(int))
        sampled = [path[i] for i in indices]
        write(outdir/"electronic_frames.xyz",sampled)
        result["status"] = "computing_electronic_path"; save()
        electronic = compute_ibo(sampled,outdir/"ibo",record["charge"],record["multiplicity"],
                                 method,basis,ordered_path=True,threads=threads)
        rm = geometry_mol(sampled[0].numbers,sampled[0].positions,record["charge"])
        pm = geometry_mol(sampled[-1].numbers,sampled[-1].positions,record["charge"])
        pop_r = np.load(outdir/"ibo/frame_0000.npz")["iao_atom_populations_e"]
        pop_p = np.load(outdir/f"ibo/frame_{len(sampled)-1:04d}.npz")["iao_atom_populations_e"]
        symbolic = arrow_hypotheses(pop_r,pop_p,rm,pm)
        symbolic["tracking_low_overlap"] = any(
            f["tracking"].get("low_overlap_flag",False) for f in electronic["frames"])
        result["symbolic"] = symbolic
        result["symbolic_status"] = "automated_hypothesis_requires_independent_review"
        result["electronic_path_kind"] = "sampled_IRC_with_relaxed_endpoints"
        result["electronic_frame_indices"] = indices.tolist()
        result["status"] = "physical_event_verified_with_symbolic_hypotheses"
    except Exception as exc:
        import traceback
        result["status"] = "calculation_failed"
        result["error"] = f"{type(exc).__name__}: {exc}"
        (outdir/"error.log").write_text(traceback.format_exc(),encoding="utf-8")
    result["gradient_evaluations"] = calc.evaluation_count
    save()
    return result
