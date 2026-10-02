"""Geometry-indexed Lewis graphs and explicitly constrained two-electron hypotheses."""
from collections import Counter
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from scipy.optimize import linear_sum_assignment


def geometry_mol(numbers, positions, charge):
    table = Chem.GetPeriodicTable()
    xyz = str(len(numbers)) + "\ngeometry\n" + "\n".join(
        f"{table.GetElementSymbol(int(z))} {x:.10f} {y:.10f} {zz:.10f}"
        for z, (x, y, zz) in zip(numbers, positions))
    mol = Chem.MolFromXYZBlock(xyz)
    rdDetermineBonds.DetermineBonds(mol, charge=int(charge), allowChargedFragments=True,
                                    embedChiral=True)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Open-shell Lewis graph outside this pilot")
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i + 1)
    return mol


def graph_smiles(mol):
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(Chem.RemoveHs(copy), isomericSmiles=True)


def resonance_equivalent(left, right):
    """Recognize enumerated Lewis resonance alternatives without moving nuclei.

    A positive result is evidence of equivalence; a truncated enumeration can
    miss alternatives. Stereochemistry and explicit hydrogen locations survive.
    """
    target=graph_smiles(right)
    if graph_smiles(left)==target:return True
    flags=Chem.UNCONSTRAINED_ANIONS|Chem.UNCONSTRAINED_CATIONS|Chem.ALLOW_CHARGE_SEPARATION
    # RDKit can fail to canonicalize an exotic resonance form ("Invariant
    # Violation", a RuntimeError); such a form cannot prove equivalence.
    try:
        structures=Chem.ResonanceMolSupplier(left,flags=flags,maxStructs=256)
    except RuntimeError:
        return False
    for m in structures:
        if m is None:continue
        try:
            if graph_smiles(m)==target:return True
        except RuntimeError:
            continue
    return False


def bond_orders(mol):
    return {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))):
            b.GetBondTypeAsDouble() for b in mol.GetBonds()}


def pair_slots(mol):
    """One token per occupied electron pair, including stationary core pairs.

    Atom tokens collect core and nonbonding pairs; bond tokens have integral
    bond-order multiplicity. Aromatic/fractional and open-shell graphs rejected.
    """
    slots = []
    bonds = bond_orders(mol)
    for bond, order in sorted(bonds.items()):
        if order != int(order):
            raise ValueError("Fractional/aromatic bond orders outside this pilot")
        slots.extend([bond] * int(order))
    for i, atom in enumerate(mol.GetAtoms()):
        nonbond_e = atom.GetAtomicNum() - atom.GetFormalCharge() - sum(
            order for bond, order in bonds.items() if i in bond)
        if nonbond_e < 0 or nonbond_e % 2 or atom.GetNumRadicalElectrons():
            raise ValueError("Lewis graph has non-paired electrons")
        slots.extend([(i,)] * int(nonbond_e / 2))
    return slots


def assign_orbitals(populations, mol):
    """Assign IBO populations to Lewis sites with exact site capacities.

    This is a graph-constrained hypothesis, not an independent arrow annotation.
    """
    pop = np.asarray(populations, dtype=float)
    slots = pair_slots(mol)
    if pop.shape != (mol.GetNumAtoms(), len(slots)):
        raise ValueError("Electron-pair/site count mismatch")
    prototypes = np.zeros_like(pop)
    for j, site in enumerate(slots):
        prototypes[list(site), j] = 2.0 / len(site)
    cost = ((pop[:, :, None] - prototypes[:, None, :]) ** 2).sum(axis=0)
    rows, columns = linear_sum_assignment(cost)
    assignment = [None] * len(slots)
    for row, column in zip(rows, columns):
        assignment[row] = slots[column]
    return assignment, float(cost[rows, columns].sum())


def arrow_hypotheses(pop_r, pop_p, mol_r, mol_p):
    sources, cost_r = assign_orbitals(pop_r, mol_r)
    sinks, cost_p = assign_orbitals(pop_p, mol_p)
    if len(sources) != len(sinks):
        raise ValueError("Electron number changed")
    arrows = [{"orbital": i, "source": list(a), "sink": list(b), "electrons": 2}
              for i, (a, b) in enumerate(zip(sources, sinks)) if a != b]
    before, after = Counter(pair_slots(mol_r)), Counter(pair_slots(mol_p))
    transported = before.copy()
    for arrow in arrows:
        transported[tuple(arrow["source"])] -= 1
        transported[tuple(arrow["sink"])] += 1
    conserved = all(transported[k] == after[k] for k in transported.keys() | after.keys())
    return {"arrows": arrows, "electron_pair_bookkeeping_passed": conserved,
            "atom_index_base": 0, "atom_map_convention": "atom_map = coordinate_index + 1",
            "assignment_cost": cost_r + cost_p,
            "label_origin": "Lewis_capacity_constrained_tracked_IBO_hypothesis",
            "independently_annotated": False, "unique_mechanism_certified": False}


def contract_atom_relays(arrows):
    """A restricted algebraic convention, not a chemical equivalence certificate.

    Contract only a balanced atom site with exactly one incoming and outgoing
    flow of equal multiplicity. Keep the original relay evidence. Do not merge
    different donor/sink pairings or remove closed cycles.
    """
    flows=Counter()
    for arrow in arrows:
        if arrow["electrons"] != 2:
            raise ValueError("This normalization only handles electron pairs")
        flows[(tuple(arrow["source"]),tuple(arrow["sink"]))]+=1
    contractions=[]
    while True:
        changed=False
        nodes=sorted({n for edge in flows for n in edge if len(n)==1})
        for node in nodes:
            incoming=[(a,n) for (a,b),n in flows.items() if b==node and a!=node and n>0]
            outgoing=[(b,n) for (a,b),n in flows.items() if a==node and b!=node and n>0]
            if len(incoming)!=1 or len(outgoing)!=1: continue
            (source,nin),(sink,nout)=incoming[0],outgoing[0]
            if nin!=nout or source==sink: continue
            del flows[(source,node)]; del flows[(node,sink)]
            flows[(source,sink)]+=nin
            contractions.append({"source":list(source),"relay_atom":list(node),
                                 "sink":list(sink),"pair_count":nin})
            changed=True; break
        if not changed: break
    normalized=[{"source":list(a),"sink":list(b),"pair_count":n}
                for (a,b),n in sorted(flows.items()) if n]
    return {"normalized_flows":normalized,"contracted_relays":contractions,
            "chemical_equivalence_certified":False}
