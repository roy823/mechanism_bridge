"""Official AIMNetCentral model families behind one evidence-recorded interface."""
import json
import os
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes


FAMILIES = {
    "aimnet2": dict(file="aimnet2_wb97m_d3_0.pt", open_shell=False),
    "aimnet2-2025": dict(file="aimnet2_2025_b973c_d3_0.pt", open_shell=False),
    "aimnet2-nse": dict(file="aimnet2nse_wb97m_0.pt", open_shell=True),
    "aimnet2-rxn": dict(file="aimnet2_rxn_0.pt", open_shell=False),
}


class AIMNetCentralPotential(Calculator):
    """Energy/force adapter that preserves official model metadata defaults."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, family, model_dir, charge=0, multiplicity=1, device="cpu",
                 compile_model=False, threads=2, needs_dispersion=None, needs_coulomb=None):
        super().__init__()
        if family not in FAMILIES:
            raise ValueError(f"Unknown AIMNetCentral family: {family}")
        import torch
        import warp
        from aimnet.calculators import AIMNet2Calculator

        self.family = family
        self.model_path = Path(model_dir) / FAMILIES[family]["file"]
        if not self.model_path.is_file():
            raise FileNotFoundError(self.model_path)
        metadata = torch.load(self.model_path, map_location="cpu", weights_only=False)
        self.supported_species = {int(z) for z in metadata["implemented_species"]}
        self.open_shell = FAMILIES[family]["open_shell"]
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.device = str(device)
        self.compile_model = bool(compile_model)
        self.needs_dispersion_override = needs_dispersion
        self.needs_coulomb_override = needs_coulomb
        warp.config.kernel_cache_dir = os.environ.get('MECHBRIDGE_WARP_CACHE',
            str(Path(model_dir).resolve().parents[1]/'.cache/warp'))
        if self.device == "cpu":
            torch.set_num_threads(threads)
        # None means: honor the model's embedded metadata. Broad models need
        # external D3; the reaction model does not. Never impose one setting.
        self.predict = AIMNet2Calculator(str(self.model_path), device=self.device,
            compile_model=self.compile_model, needs_dispersion=needs_dispersion,
            needs_coulomb=needs_coulomb)
        self.metadata = dict(cutoff=float(metadata["cutoff"]),
            implemented_species=sorted(self.supported_species),
            needs_dispersion=bool(metadata["needs_dispersion"]),
            needs_coulomb=bool(metadata["needs_coulomb"]),
            family=metadata.get("family", family))
        self.validate_system([], self.charge, self.multiplicity, allow_empty=True)

    def validate_system(self, numbers, charge, multiplicity, allow_empty=False):
        numbers = np.asarray(numbers, dtype=int)
        if not allow_empty and not len(numbers):
            raise ValueError("A nonempty molecular system is required")
        unsupported = sorted(set(numbers) - self.supported_species)
        if unsupported:
            raise ValueError(f"Elements outside {self.family}: {unsupported}")
        if float(charge) != int(charge) or float(multiplicity) != int(multiplicity):
            raise ValueError("Integer charge and spin multiplicity are required")
        if int(multiplicity) < 1:
            raise ValueError("Spin multiplicity must be positive")
        if not self.open_shell and int(multiplicity) != 1:
            raise ValueError(f"{self.family} adapter only admits closed-shell singlets")
        if self.family == "aimnet2-rxn" and int(charge) != 0:
            raise ValueError("AIMNet2-rxn only supports net-neutral systems")
        electrons=int(numbers.sum())-int(charge)
        unpaired=int(multiplicity)-1
        if len(numbers) and (electrons < unpaired or (electrons-unpaired)%2):
            raise ValueError("Charge and multiplicity are inconsistent with electron parity")
        self.charge, self.multiplicity = int(charge), int(multiplicity)

    def _input(self, numbers, positions):
        positions = np.asarray(positions)
        batch = positions.ndim == 3
        data = dict(coord=positions, numbers=(np.broadcast_to(numbers,
            (len(positions), len(numbers))).copy() if batch else np.asarray(numbers)),
            charge=(np.full(len(positions), self.charge, dtype=float) if batch else float(self.charge)))
        if self.open_shell:
            data["mult"] = (np.full(len(positions), self.multiplicity, dtype=float)
                            if batch else float(self.multiplicity))
        return data

    def evaluate_many(self, numbers, positions):
        numbers = np.asarray(numbers, dtype=int);positions = np.asarray(positions, dtype=float)
        if (positions.ndim != 3 or positions.shape[1:] != (len(numbers), 3)
            or not len(positions) or not np.isfinite(positions).all()):
            raise ValueError("Batch requires finite equal-composition geometries")
        self.validate_system(numbers, self.charge, self.multiplicity)
        import torch
        if self.compile_model:
            result = self.predict(self._input(numbers, positions), forces=True)
        else:
            with torch.compiler.set_stance("force_eager"):
                result = self.predict(self._input(numbers, positions), forces=True)
        energy = result["energy"].detach().cpu().numpy().reshape(-1)
        forces = result["forces"].detach().cpu().numpy().reshape(positions.shape)
        if not np.isfinite(energy).all() or not np.isfinite(forces).all():
            raise ValueError("Nonfinite AIMNetCentral prediction")
        return {"energy": energy, "forces": forces}

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if self.atoms.constraints or np.any(self.atoms.pbc):
            raise ValueError("Reaction exploration requires unconstrained nonperiodic molecules")
        result = self.evaluate_many(self.atoms.numbers, self.atoms.positions[None])
        self.results = {"energy": float(result["energy"][0]),
                        "forces": result["forces"][0].copy()}

    def describe(self):
        result=dict(self.metadata)
        result.update(family=self.family,device=self.device,compile_model=self.compile_model,
                      model_path=str(self.model_path),
                      needs_dispersion_override=self.needs_dispersion_override,
                      needs_coulomb_override=self.needs_coulomb_override)
        return json.loads(json.dumps(result))
