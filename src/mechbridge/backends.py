"""ASE calculator contracts; PySCF backend is for explicit small closed-shell jobs."""
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

def configure_pyscf(threads):
    from pyscf import lib
    lib.num_threads(threads)
    try:
        lib.current_memory()
    except FileNotFoundError:
        # Some hosted containers lack /proc/<pid>/statm. Only memory telemetry
        # is replaced; energies, gradients and convergence criteria are unchanged.
        import resource
        def memory_without_proc():
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024., 0.
        lib.current_memory = memory_without_proc
        lib.misc.current_memory = memory_without_proc

class PySCFCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, charge, multiplicity, method="wb97x", basis="6-31g(d)", threads=1,
                 reuse_density=False, **kwargs):
        super().__init__(**kwargs)
        if multiplicity != 1:
            raise ValueError("v0.1 PySCF backend only supports closed-shell singlets")
        self.charge, self.multiplicity = charge, multiplicity
        self.method, self.basis, self.threads = method, basis, threads
        self.evaluation_count = 0
        self.mean_field = None
        self.reuse_density = reuse_density

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        from pyscf import gto, scf, dft, lib
        import numpy as np
        if self.atoms.constraints or np.any(self.atoms.pbc):
            raise ValueError("Only unconstrained gas-phase molecules supported")
        configure_pyscf(self.threads)
        m = gto.M(atom=list(zip(self.atoms.get_chemical_symbols(), self.atoms.positions.tolist())),
                  unit="Angstrom", charge=self.charge, spin=0, basis=self.basis, verbose=0)
        mf = scf.RHF(m) if self.method.upper() == "HF" else dft.RKS(m, xc=self.method)
        mf.conv_tol = 1e-10; mf.max_cycle = 150
        dm0 = None
        if (self.reuse_density and self.mean_field is not None and m.nelectron == self.mean_field.mol.nelectron
                and m.nao_nr() == self.mean_field.mol.nao_nr()
                and np.array_equal(m.atom_charges(), self.mean_field.mol.atom_charges())):
            dm0 = self.mean_field.make_rdm1()
        energy = mf.kernel(dm0=dm0)
        if not mf.converged:
            raise RuntimeError("SCF did not converge; do not use these labels")
        gradient = mf.nuc_grad_method().kernel()
        self.mean_field = mf
        self.evaluation_count += 1
        self.results = {"energy": float(energy * Hartree), "forces": -gradient * Hartree / Bohr}

    def hessian(self, atoms):
        """Analytic Cartesian Hessian on the same converged PES, in eV/A^2."""
        import numpy as np
        self.get_potential_energy(atoms)
        raw = np.asarray(self.mean_field.Hessian().kernel())
        result = raw.transpose(0, 2, 1, 3).reshape(3 * len(atoms), 3 * len(atoms))
        result = result * Hartree / Bohr**2
        if not np.isfinite(result).all():
            raise RuntimeError("Nonfinite analytic Hessian")
        return (result + result.T) / 2

def load_factory(spec, config):
    """User-provided MLIP module:function must return a fresh ASE Calculator.

    Its charge/spin/element/domain checks are the plugin author's responsibility.
    The project never silently substitutes a neutral or lower-level model.
    """
    import importlib
    module, function = spec.split(":", 1)
    fn = getattr(importlib.import_module(module), function)
    return lambda: fn(**config)
