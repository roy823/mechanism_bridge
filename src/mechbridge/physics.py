"""Finite-difference curvature and two-sided minimization on a supplied ASE PES.

No chemistry rules. A downhill minimization is explicitly NOT an IRC integrator.
"""
import numpy as np
from scipy.linalg import null_space
from ase.optimize import BFGS

FREQUENCY_FACTOR = 521.470898  # sqrt(eV / Angstrom^2 / amu) -> cm^-1

def finite_hessian(atoms, step=0.005):
    if step <= 0 or atoms.constraints or np.any(atoms.pbc):
        raise ValueError("Requires positive displacement and unconstrained nonperiodic molecule")
    x = atoms.get_positions().copy()
    h = np.empty((x.size, x.size))
    try:
        for j in range(x.size):
            delta = np.zeros(x.size); delta[j] = step
            atoms.set_positions(x + delta.reshape(x.shape)); fp = atoms.get_forces().ravel()
            atoms.set_positions(x - delta.reshape(x.shape)); fm = atoms.get_forces().ravel()
            h[:, j] = -(fp - fm) / (2 * step)
    finally:
        atoms.set_positions(x)
    if not np.isfinite(h).all():
        raise ValueError("Nonfinite Hessian")
    return (h + h.T) / 2

def vibrational_analysis(positions, masses, hessian, imaginary_threshold_cm=30.):
    x, masses, h = np.asarray(positions), np.asarray(masses), np.asarray(hessian)
    n = len(masses)
    if x.shape != (n, 3) or h.shape != (3*n, 3*n) or np.any(masses <= 0):
        raise ValueError("Invalid geometry, mass or Hessian dimensions")
    if not (np.isfinite(x).all() and np.isfinite(h).all()):
        raise ValueError("Nonfinite input")
    root = np.sqrt(masses)[:, None]
    centered = x - np.average(x, axis=0, weights=masses)
    rigid = []
    for axis in np.eye(3):
        rigid.append((np.ones_like(x) * axis * root).ravel())
        rigid.append((np.cross(centered, axis) * root).ravel())
    rigid = np.array(rigid).T
    # Linear molecules correctly have five rigid modes, nonlinear molecules six.
    q = null_space(rigid.T, rcond=1e-10)
    m = np.repeat(np.sqrt(masses), 3)
    hw = h / m[:, None] / m[None, :]
    values, vectors = np.linalg.eigh(q.T @ hw @ q)
    freq = np.sign(values) * np.sqrt(np.abs(values)) * FREQUENCY_FACTOR
    modes = ((q @ vectors) / m[:, None]).T.reshape(-1, n, 3)
    modes /= np.maximum(np.linalg.norm(modes.reshape(len(values), -1), axis=1)[:, None, None], 1e-30)
    return {"frequencies_cm-1": freq.tolist(), "imaginary_count": int(np.sum(freq < -imaginary_threshold_cm)),
            "imaginary_threshold_cm-1": imaginary_threshold_cm,
            "rigid_mode_count": int(3*n - q.shape[1]), "modes": modes}

def analyze_stationary(atoms, fmax=0.02, step=0.005, imaginary_threshold_cm=30.):
    force = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    vibration = vibrational_analysis(atoms.positions, atoms.get_masses(), finite_hessian(atoms, step), imaginary_threshold_cm)
    return {"force_max_eV_A": force, "force_converged": force <= fmax,
            "energy_eV": float(atoms.get_potential_energy()), **vibration}

def verify_saddle(atoms, calculator_factory, outdir, fmax=0.02, steps=300,
                  displacement=0.10, hessian_step=0.005):
    from pathlib import Path
    from ase.io import write
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    atoms.calc = calculator_factory()
    ts = analyze_stationary(atoms, fmax, hessian_step)
    modes = ts.pop("modes")
    result = {"ts": ts, "method": "negative_mode_displacement_then_BFGS",
              "is_IRC": False, "expected_endpoint_match": "not_checked",
              "status": "unresolved", "endpoints": []}
    if not ts["force_converged"] or ts["imaginary_count"] != 1:
        result["status"] = "not_converged_index1_at_this_geometry"
        return result
    write(outdir / "ts.xyz", atoms)
    minima = []
    for sign in [-1, 1]:
        a = atoms.copy(); a.calc = calculator_factory()
        a.positions += sign * displacement * modes[0]
        stem = "minus" if sign == -1 else "plus"
        opt = BFGS(a, logfile=str(outdir / f"{stem}.log"), trajectory=str(outdir / f"{stem}.traj"))
        opt.run(fmax=fmax, steps=steps)
        end = analyze_stationary(a, fmax, hessian_step); end.pop("modes")
        end["geometry_file"] = stem + ".xyz"
        end["barrier_electronic_eV"] = ts["energy_eV"] - end["energy_eV"]
        result["endpoints"].append(end); minima.append(a.positions.copy())
        write(outdir / (stem + ".xyz"), a)
    # Rotation/translation-invariant, fixed atom-index distance comparison.
    d0 = np.linalg.norm(minima[0][:, None] - minima[0][None, :], axis=-1)
    d1 = np.linalg.norm(minima[1][:, None] - minima[1][None, :], axis=-1)
    result["endpoint_distance_matrix_rms_A"] = float(np.sqrt(np.mean((d0 - d1)**2)))
    valid = all(e["force_converged"] and e["imaginary_count"] == 0 and e["barrier_electronic_eV"] >= -1e-5 for e in result["endpoints"])
    result["status"] = "index1_with_two_relaxed_minima" if valid else "descent_unresolved"
    result["distinct_basins"] = "not_certified_geometry_metric_only"
    return result

def dimer_refine(atoms, outdir, fmax=0.02, steps=300, seed=17):
    from pathlib import Path
    from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    with DimerControl(initial_eigenmode_method="gauss", displacement_method="gauss",
                      logfile=str(outdir / "dimer.log")) as control:
        minimum = MinModeAtoms(atoms, control=control, random_seed=seed)
        with MinModeTranslate(minimum, logfile=str(outdir / "dimer_opt.log"),
                              trajectory=str(outdir / "dimer.traj")) as opt:
            converged = opt.run(fmax=fmax, steps=steps)
    return bool(converged)
