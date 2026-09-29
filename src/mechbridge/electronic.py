"""Compute IBO/IAO populations; optional overlap tracking on user-supplied ordered paths.

No hardcoded conversion from bond differences to arrow labels.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

def compute_ibo(frames, outdir, charge, multiplicity, method="wb97x", basis="6-31g(d)",
                ordered_path=False, threads=1):
    from pyscf import gto, scf, dft, lib
    from pyscf.lo import iao, ibo, orth
    if multiplicity != 1:
        raise ValueError("IBO prototype supports closed-shell singlets only")
    frames = list(frames)
    if not frames:
        raise ValueError("No frames")
    from .backends import configure_pyscf
    configure_pyscf(threads)
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    numbers = frames[0].numbers.tolist()
    previous = None; report = []
    for i, atoms in enumerate(frames):
        if atoms.numbers.tolist() != numbers or np.any(atoms.pbc):
            raise ValueError("Same atom order and nonperiodic geometry required")
        m = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions.tolist())),
                  charge=charge, spin=0, basis=basis, unit="Angstrom", verbose=0)
        mf = scf.RHF(m) if method.upper() == "HF" else dft.RKS(m, xc=method)
        mf.conv_tol = 1e-10; mf.max_cycle = 150
        energy = mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"SCF failed at frame {i}")
        s = m.intor_symmetric("int1e_ovlp")
        occupied = mf.mo_coeff[:, mf.mo_occ > 0]
        iaos = orth.vec_lowdin(iao.iao(m, occupied), s)
        # PM with exponent four in the IAO representation allows convergence auditing.
        localizer = ibo.PM(m, occupied, iaos, s, exponent=4)
        localizer.exponent = 4
        localizer.conv_tol = 1e-9
        localizer.conv_tol_grad = 1e-5
        convergence = {"converged": False}
        def callback(state):
            convergence["converged"] = bool(state.get("conv", False))
        coeff = localizer.kernel(callback=callback)
        localization_gradient = float(np.linalg.norm(localizer.get_grad()))
        if not convergence["converged"] or localization_gradient > 1e-5:
            raise RuntimeError(f"IBO localization failed at frame {i}")
        tracking = {"performed": False}
        if ordered_path and previous is not None:
            pmol, pc = previous
            overlap = pc.T @ gto.intor_cross("int1e_ovlp", pmol, m) @ coeff
            rows, cols = linear_sum_assignment(-np.abs(overlap)**2)
            permutation = cols[np.argsort(rows)]
            coeff = coeff[:, permutation]
            signed = np.diag(pc.T @ gto.intor_cross("int1e_ovlp", pmol, m) @ coeff)
            coeff *= np.where(signed < 0, -1, 1)[None, :]
            tracking = {"performed": True, "permutation": permutation.tolist(),
                        "squared_overlaps": (signed**2).tolist(),
                        "low_overlap_flag": bool(np.min(signed**2) < 0.5)}
        projection = iaos.T @ s @ coeff
        reference = iao.reference_mol(m)
        populations = np.stack([2 * np.sum(projection[p0:p1]**2, axis=0)
                                for _, _, p0, p1 in reference.aoslice_by_atom()])
        error = float(np.max(np.abs(coeff.T @ s @ coeff - np.eye(coeff.shape[1]))))
        np.savez_compressed(outdir / f"frame_{i:04d}.npz", positions_A=atoms.positions,
                            atomic_numbers=atoms.numbers, ibo_coeff_ao=coeff, overlap_ao=s,
                            iao_atom_populations_e=populations)
        report.append({"frame": i, "energy_hartree": float(energy), "scf_converged": True,
                       "localization_converged": convergence["converged"],
                       "localization_gradient_norm": localization_gradient,
                       "orthogonality_error": error, "tracking": tracking,
                       "orbital_population_sum_e": populations.sum(axis=0).tolist()})
        previous = (m, coeff)
    result = {"method": method, "basis": basis, "charge": charge, "multiplicity": multiplicity,
              "environment": "gas_phase", "ordered_path_asserted_by_user": ordered_path,
              "orbital_representation": "Pipek-Mezey exponent 4 in IAO basis (PySCF ibo.PM)",
              "scf_stability_test": "not_performed", "arrows": None,
              "arrow_status": "requires_annotation_or_trained_decoder", "frames": report}
    (outdir / "electronic_report.json").write_text(json.dumps(result, indent=2))
    return result
