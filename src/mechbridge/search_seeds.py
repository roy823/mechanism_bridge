"""Geometry proposals and optional electron-source/sink coupling, without PES bias."""
import numpy as np
from scipy.optimize import least_squares
from ase.data import covalent_radii
from .event_graph import bond_orders


def internal_direction(x, direction):
    """Project translations and infinitesimal rotations out of Cartesian proposals."""
    c = x - x.mean(axis=0)
    rigid = np.stack([np.tile(v, (len(x), 1)).ravel() for v in np.eye(3)] +
                     [np.cross(c, v).ravel() for v in np.eye(3)], axis=1)
    d = direction.ravel()
    d = d - rigid @ np.linalg.lstsq(rigid, d, rcond=1e-10)[0]
    norm = np.linalg.norm(d)
    if norm < 1e-10:
        raise ValueError('Degenerate seed direction')
    return (d / norm).reshape(x.shape)


def make_seed(atoms, mol, strategy, proposal, sample, random_seed):
    x = atoms.positions.copy()
    rng = np.random.default_rng(random_seed)
    fraction = (0.35, 0.55, 0.75)[sample % 3]
    amplitude = (0.6, 1.0, 1.4)[sample % 3]
    if strategy in ('geometry', 'center_random'):
        # Internal random pair-distance perturbations, independent of arrows/target.
        active = set(range(len(x)))
        if strategy == 'center_random':
            if proposal is None:
                raise ValueError('Center control requires proposed active atom identities')
            active = {i for edit in proposal['edits'] for i in edit['atoms']}
        pairs = [(i, j) for i in sorted(active) for j in sorted(active) if j<i and np.linalg.norm(x[i]-x[j]) < 3.5]
        if not pairs:
            raise ValueError('No local pairs for geometry seed')
        d = np.zeros_like(x)
        for k in rng.choice(len(pairs), size=min(3, len(pairs)), replace=False):
            i, j = pairs[k]
            v = x[i] - x[j]
            v /= max(np.linalg.norm(v), 1e-8)
            v *= rng.choice([-1., 1.])
            d[i] += v
            d[j] -= v
        d = internal_direction(x, d)
        return x + amplitude * d, d, dict(fraction=fraction, template_id=None,
            displacement_norm_A=amplitude, active_atoms=sorted(active),
            information='active_atom_ids_only' if strategy=='center_random' else 'geometry_only')
    if strategy not in ('bond_edits', 'arrows') or proposal is None:
        raise ValueError('Symbolic strategy needs an applicable proposal')
    edits = proposal['edits']
    old_bonds = bond_orders(mol)
    changed = {tuple(e['atoms']) for e in edits}
    ij = np.array([e['atoms'] for e in edits])
    r0 = np.linalg.norm(x[ij[:, 0]] - x[ij[:, 1]], axis=1)
    target = []
    for e, distance in zip(edits, r0):
        i, j = e['atoms']
        if e['after'] == 0:
            target.append(distance + 1.1)
        else:
            length = covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]]
            target.append(length * {1: 1., 2: .88, 3: .80}[e['after']])
    delta = np.array(target) - r0
    delta = np.where(np.abs(delta) < .08, np.where(delta >= 0, .08, -.08), delta)
    unchanged = [(i, j) for i, j in old_bonds if (i, j) not in changed]
    links = []
    index = {tuple(e['atoms']): k for k, e in enumerate(edits)}
    if strategy == 'arrows':
        for arrow in proposal['arrows']:
            a, b = tuple(sorted(arrow['source'])), tuple(sorted(arrow['sink']))
            if a in index and b in index:
                links.append((index[a], index[b]))
    # This synchrony term is an explicit proposal heuristic, not a physical law.
    def residual(flat):
        y = flat.reshape(x.shape)
        distances = np.linalg.norm(y[ij[:, 0]] - y[ij[:, 1]], axis=1)
        progress = (distances - r0) / delta
        out = list(2. * (distances - r0 - fraction * delta))
        out.extend((np.linalg.norm(y[i]-y[j]) - np.linalg.norm(x[i]-x[j])) for i,j in unchanged)
        out.extend((.12 * (y-x)).ravel())
        out.extend(0.7 * (progress[a] - progress[b]) for a,b in links)
        for i in range(len(y)):
            for j in range(i):
                floor = .6 * (covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]])
                out.append(3 * max(0., floor - np.linalg.norm(y[i]-y[j])))
        return out
    initial = x + rng.normal(0., .015, x.shape)
    fit = least_squares(residual, initial.ravel(), max_nfev=100)
    y = fit.x.reshape(x.shape)
    direction = internal_direction(x, y-x)
    return x + amplitude * direction, direction, dict(fraction=fraction, template_id=proposal['template_id'],
                            displacement_norm_A=amplitude,
                            unnormalized_fit_displacement_A=float(np.linalg.norm(y-x)),
                            source_sink_couplings=links, geometric_fit_cost=float(fit.cost),
                            predicted_graph=proposal['predicted_graph'],
                            arrows=proposal['arrows'] if strategy=='arrows' else None,
                            net_edits=edits)
