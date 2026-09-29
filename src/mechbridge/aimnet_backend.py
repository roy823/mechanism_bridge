"""Pinned AIMNet2-rxn inference through the official calculator (no retraining)."""
import numpy as np
import json
from pathlib import Path
from ase.calculators.calculator import Calculator, all_changes


class ReactionPotential(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, model_path, charge=0, multiplicity=1, threads=2):
        super().__init__()
        if charge != 0 or multiplicity != 1:
            raise ValueError('AIMNet2-rxn pilot requires charge=0 and multiplicity=1')
        import torch
        import warp
        warp.config.kernel_cache_dir = str(Path(model_path).resolve().parents[1]/'.cache/warp')
        from aimnet.calculators import AIMNet2Calculator
        torch.set_num_threads(threads)
        config = json.loads((Path(model_path)/'config.json').read_text())
        if config['needs_dispersion'] or config['implemented_species'] != [1,6,7,8]:
            raise ValueError('Unexpected AIMNet2-rxn artifact configuration')
        # Explicitly honor the pinned model card: VV10 is in training labels.
        # aimnet 0.2 family defaults otherwise add D3 despite this artifact flag.
        self.predict = AIMNet2Calculator(str(model_path), device='cpu', compile_model=False,
                                         needs_dispersion=False, needs_coulomb=True)

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        a = self.atoms
        if a.constraints or np.any(a.pbc) or not set(a.numbers) <= {1,6,7,8}:
            raise ValueError('Only unconstrained nonperiodic CHNO systems are supported')
        if sum(a.numbers) % 2:
            raise ValueError('Odd electron count is outside the closed-shell pilot')
        import torch
        # Includes nvalchemiops helper functions decorated with torch.compile.
        with torch.compiler.set_stance('force_eager'):
            r = self.predict(dict(coord=a.positions, numbers=a.numbers, charge=0.), forces=True)
        self.results = dict(energy=float(r['energy'].detach().cpu().reshape(-1)[0]),
                            forces=r['forces'].detach().cpu().numpy().reshape((-1,3)))

    def evaluate_many(self, numbers, positions):
        """Independent same-composition systems via the official batch API."""
        numbers=np.asarray(numbers)
        positions=np.asarray(positions)
        if (positions.ndim!=3 or positions.shape[1:]!=(len(numbers),3)
            or not len(positions) or not np.isfinite(positions).all()
            or not set(numbers)<={1,6,7,8} or sum(numbers)%2):
            raise ValueError('Batch requires finite neutral closed-shell CHNO geometries')
        import torch
        with torch.compiler.set_stance('force_eager'):
            r=self.predict(dict(coord=positions,
                numbers=np.broadcast_to(numbers,(len(positions),len(numbers))).copy(),
                charge=np.zeros(len(positions))),forces=True)
        return dict(energy=r['energy'].detach().cpu().numpy().reshape(-1),
                    forces=r['forces'].detach().cpu().numpy().reshape(positions.shape))
