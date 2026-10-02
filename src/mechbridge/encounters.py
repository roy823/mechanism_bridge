"""Rigid fragment assembly; no reference products, TSs, or energy labels."""
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import brentq
from ase.data import vdw_radii, covalent_radii
from rdkit import Chem


def assemble_encounter(numbers, positions, fragments, random_seed):
    if len(fragments) != 2:
        raise ValueError('Encounter assembly requires exactly two fragments')
    rng = np.random.default_rng(random_seed)
    x = np.array(positions, dtype=float, copy=True)
    a, b = [np.array(f, dtype=int) for f in fragments]
    for ids in (a,b):
        x[ids] = (x[ids]-x[ids].mean(0)) @ Rotation.random(random_state=rng).as_matrix()
    direction = rng.normal(size=3); direction /= np.linalg.norm(direction)
    radii = vdw_radii[np.asarray(numbers)]
    if not np.isfinite(radii).all():
        raise ValueError('Missing van der Waals radius')
    def gap(distance):
        d = np.linalg.norm(x[a,None]-x[b][None,:]-distance*direction,axis=-1)
        return float(np.min(d-radii[a,None]-radii[b][None,:])-.2)
    upper = 2*np.max(np.linalg.norm(x,axis=1))+8.
    if gap(0) >= 0:
        raise ValueError('Overlapping-center encounter unexpectedly has no contact')
    distance = brentq(gap,0,upper)
    x[b] += distance*direction
    return x-x.mean(0)


def orient_reactive_encounter(numbers, positions, mol, edits, random_seed):
    """Sample rigid placements exposing proposed cross-fragment forming bonds.

    The 64 inexpensive geometric trials do not evaluate the PES. This rotation
    and translation is recorded separately from the matched deformation norm.
    """
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) < 2:
        return np.array(positions,copy=True), None
    component = {atom: index for index, fragment in enumerate(fragments) for atom in fragment}
    pairs = []
    for e in edits:
        i,j = e['atoms']
        if e['before']==0 and component[i] != component[j]:
            pairs.append((i,j))
    if not pairs:
        return np.array(positions,copy=True), None
    active_components = tuple(sorted((component[pairs[0][0]], component[pairs[0][1]])))
    pairs = [pair for pair in pairs
             if tuple(sorted((component[pair[0]], component[pair[1]]))) == active_components]
    ids_a, ids_b = [np.array(fragments[index],dtype=int) for index in active_components]
    aset = set(ids_a)
    pairs = [(i,j) if i in aset else (j,i) for i,j in pairs]
    numbers = np.asarray(numbers)
    pairs.sort(key=lambda p: (numbers[p[0]]==1 or numbers[p[1]]==1,p))
    rng = np.random.default_rng(random_seed)
    original = np.asarray(positions)
    best = _best_rigid_orientation(numbers, original, fragments, active_components, pairs, rng)
    return best[1],dict(policy='64_rigid_orientations_cross_bond_distance_and_clash_score',
        score=float(best[0]),selected_trial=best[2],cross_forming_pairs=pairs,
        active_components=list(active_components),spectator_components=len(fragments)-2,
        contact_distance_A=2.8,reference_product_or_TS_used=False,
        displacement_from_source_A=float(np.linalg.norm(best[1]-original)))


def orient_control_encounter(numbers, positions, mol, allowed_atoms, random_seed):
    """Proposal-free rigid placement for the geometry and center_random controls.

    Same 64-trial pool, 2.8 A contact and clash score as the symbolic encounter,
    but the contact pair is drawn uniformly from cross-fragment pairs of
    `allowed_atoms` (heavy-atom pairs first) instead of a proposed forming bond.
    """
    fragments = Chem.GetMolFrags(mol)
    if len(fragments) < 2:
        return np.array(positions,copy=True), None
    component = {atom: index for index, fragment in enumerate(fragments) for atom in fragment}
    numbers = np.asarray(numbers)
    allowed = sorted(allowed_atoms)
    cross = [(i,j) for i in allowed for j in allowed if i<j and component[i]!=component[j]]
    heavy = [p for p in cross if numbers[p[0]]>1 and numbers[p[1]]>1]
    candidates = heavy or cross
    if not candidates:
        return np.array(positions,copy=True), None
    rng = np.random.default_rng(random_seed)
    i,j = candidates[int(rng.integers(len(candidates)))]
    active_components = tuple(sorted((component[i], component[j])))
    pair = (i,j) if component[i]==active_components[0] else (j,i)
    original = np.asarray(positions)
    best = _best_rigid_orientation(numbers, original, fragments, active_components, [pair], rng)
    return best[1],dict(policy='64_rigid_orientations_random_cross_pair_contact_and_clash_score',
        score=float(best[0]),selected_trial=best[2],contact_pair=list(map(int,pair)),
        contact_pair_pool='heavy_atom_cross_pairs' if heavy else 'all_cross_pairs',
        candidate_pairs=len(candidates),active_components=list(active_components),
        spectator_components=len(fragments)-2,contact_distance_A=2.8,
        proposal_bond_used=False,reference_product_or_TS_used=False,
        displacement_from_source_A=float(np.linalg.norm(best[1]-original)))


def _best_rigid_orientation(numbers, original, fragments, active_components, pairs, rng):
    """Score 64 rigid placements that put pairs[0] at 2.8 A; no PES evaluation."""
    ids_a, ids_b = [np.array(fragments[index],dtype=int) for index in active_components]
    anchor_a,anchor_b = pairs[0]
    active_ids = np.concatenate((ids_a,ids_b))
    active_center = original[active_ids].mean(0)
    best = None
    for trial in range(64):
        y = original.copy()
        for ids,anchor in ((ids_a,anchor_a),(ids_b,anchor_b)):
            y[ids] = (y[ids]-y[anchor]) @ Rotation.random(random_state=rng).as_matrix()
        y[ids_b] += np.array([0.,0.,2.8])
        y[active_ids] += active_center-y[active_ids].mean(0)
        d = np.linalg.norm(y[ids_a,None]-y[ids_b][None,:],axis=-1)
        floor = 1.15*(covalent_radii[numbers[ids_a,None]]+covalent_radii[numbers[ids_b]][None,:])
        score = sum((np.linalg.norm(y[i]-y[j])-2.8)**2 for i,j in pairs)
        score += 50*np.maximum(floor-d,0).sum()**2
        spectators = np.array([i for index,fragment in enumerate(fragments)
                               if index not in active_components for i in fragment],dtype=int)
        if len(spectators):
            d_other=np.linalg.norm(y[active_ids,None]-y[spectators][None,:],axis=-1)
            floor_other=(covalent_radii[numbers[active_ids,None]]+
                         covalent_radii[numbers[spectators]][None,:])
            score += 50*np.maximum(floor_other-d_other,0).sum()**2
        if best is None or score < best[0]:
            best = score,y-y.mean(0),trial
    return best
