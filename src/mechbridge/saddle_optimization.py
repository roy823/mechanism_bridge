"""ASE Dimer with a convergence norm matching the physical stationary-point audit."""
import numpy as np
from ase.mep import MinModeTranslate


class StationaryDimerTranslate(MinModeTranslate):
    """Keep ASE steps/rotations; stop on original per-atom force, not projected force.

    Mode reflection preserves the global force norm but not each atom's norm.
    Using projected per-atom maxima can stop just before the physical fmax test.
    """
    def gradient_converged(self, gradient):
        force=self.dimeratoms.get_forces(real=True)
        return bool(np.linalg.norm(force,axis=1).max()<=self.fmax and
                    self.dimeratoms.get_curvature()<0.)
