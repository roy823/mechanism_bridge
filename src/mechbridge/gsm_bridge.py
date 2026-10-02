"""Budgeted ASE calculator for the external single-ended GSM baseline (pyGSM).

pyGSM builds its ASE calculator from an import path and JSON keyword arguments
inside its own process. This class loads the same pinned MLIP, bills every
evaluation exactly like CountedCalculator, writes the running count to
counter_file after every call (so the count survives pyGSM stopping on the
budget), and raises BudgetExceeded at the limit.
"""
from pathlib import Path

from ase.calculators.calculator import Calculator, all_changes

from .potentials import load_potential
from .reaction_network import CountedCalculator


class BudgetedPotential(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, potential, root, limit, counter_file, charge=0, multiplicity=1, threads=2):
        super().__init__()
        backend, _ = load_potential(potential, root, threads=int(threads))
        self._validator = getattr(backend, 'validate_system', None)
        self._state = (int(charge), int(multiplicity))
        self._validated = False
        self.counted = CountedCalculator(backend, int(limit))
        self.counter_file = Path(counter_file)
        self.counter_file.write_text('0')

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if not self._validated and self._validator is not None:
            self._validator(self.atoms.numbers, *self._state)
            self._validated = True
        try:
            self.counted.calculate(self.atoms, properties, system_changes)
        finally:
            self.counter_file.write_text(str(self.counted.calls))
        self.results = dict(energy=self.counted.results['energy'],
                            forces=self.counted.results['forces'].copy())
