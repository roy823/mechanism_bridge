"""Geometry proposals and optional electron-source/sink coupling, without PES bias."""
import numpy as np
from scipy.optimize import least_squares
from ase.data import covalent_radii
from .event_graph import bond_orders
from .encounters import orient_reactive_encounter


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


def cosine_angle(x, ids):
    a, b, c = x[list(ids)]
    u, v = a-b, c-b
    return float(np.dot(u, v) / max(np.linalg.norm(u)*np.linalg.norm(v), 1e-12))


def angle_targets(numbers, edits, sample):
    """Soft, sampled geometrical priors, available to both symbolic controls."""
    targets = []
    broken = [e for e in edits if e['after'] == 0]
    formed = [e for e in edits if e['before'] == 0]
    for old in broken:
        if old['before']!=1:continue
        for new in formed:
            shared=set(old['atoms']) & set(new['atoms'])
            if len(shared)!=1:continue
            center=next(iter(shared))
            leave=next(i for i in old['atoms'] if i!=center)
            donor=next(i for i in new['atoms'] if i!=center)
            if numbers[center]==6 and numbers[leave] in (7,8) and numbers[donor] in (7,8):
                targets.append(dict(atoms=[donor,center,leave],degrees=(150.,165.,180.)[sample%3],
                                    kind='nucleophilic_substitution_approach'))
    for old in broken:
        for new in formed:
            shared = set(old['atoms']) & set(new['atoms'])
            if len(shared) == 1:
                h = next(iter(shared))
                if numbers[h] == 1:
                    donor = next(i for i in old['atoms'] if i != h)
                    acceptor = next(i for i in new['atoms'] if i != h)
                    angle = (110., 140., 170.)[sample % 3]
                    targets.append(dict(atoms=[donor,h,acceptor], degrees=angle,
                                        kind='hydrogen_transfer'))
    for edit in edits:
        if edit['before'] != 2 or edit['after'] != 1:
            continue
        for c in edit['atoms']:
            hetero = next(i for i in edit['atoms'] if i != c)
            if numbers[c] != 6 or numbers[hetero] not in (7,8):
                continue
            for new in formed:
                if c in new['atoms']:
                    nu = next(i for i in new['atoms'] if i != c)
                    targets.append(dict(atoms=[nu,c,hetero],degrees=(100.,110.,120.)[sample % 3],
                                        kind='carbonyl_addition'))
    return targets


def reaction_direction(x, ij, delta, unchanged):
    """Damped internal-coordinate tangent at the actual seed geometry."""
    rows, rhs = [], []
    for (i,j), value, weight in [(p, v, 1.) for p,v in zip(ij,delta)] + [
                                (p, 0., .5) for p in unchanged]:
        row = np.zeros_like(x)
        unit = (x[i]-x[j])/max(np.linalg.norm(x[i]-x[j]),1e-10)
        row[i], row[j] = weight*unit, -weight*unit
        rows.append(row.ravel()); rhs.append(weight*value)
    jac = np.asarray(rows)
    tangent = np.linalg.lstsq(np.vstack([jac, .1*np.eye(x.size)]),
                             np.r_[rhs, np.zeros(x.size)],rcond=1e-10)[0].reshape(x.shape)
    return internal_direction(x, tangent)


def make_seed(atoms, mol, strategy, proposal, sample, random_seed, symbolic_seed_scale=1.0):
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
    if symbolic_seed_scale <= 0:
        raise ValueError('Symbolic seed scale must be positive')
    amplitude *= symbolic_seed_scale
    edits = proposal['edits']
    x, encounter = orient_reactive_encounter(atoms.numbers,x,mol,edits,random_seed)
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
    # Matched across net-edit and arrow controls: synchronous, decrease-leading,
    # and increase-leading trajectories. These are proposals, not assumed paths.
    phase = (0., -.15, .15)[sample % 3]
    progress_target = np.clip(fraction + phase*np.array([
        np.sign(e['after']-e['before']) for e in edits]), .15, .9)
    angles = angle_targets(atoms.numbers, edits, sample)
    unchanged = [(i, j) for i, j in old_bonds if (i, j) not in changed]
    links = []
    index = {tuple(e['atoms']): k for k, e in enumerate(edits)}
    if strategy == 'arrows':
        for arrow in proposal['arrows']:
            a, b = tuple(sorted(arrow['source'])), tuple(sorted(arrow['sink']))
            if a in index and b in index:
                links.append((index[a], index[b]))
    # Source/sink progress differences follow the sampled asynchronous targets.
    def residual(flat):
        y = flat.reshape(x.shape)
        distances = np.linalg.norm(y[ij[:, 0]] - y[ij[:, 1]], axis=1)
        progress = (distances - r0) / delta
        out = list((distances - r0 - progress_target * delta) / .5)
        out.extend((np.linalg.norm(y[i]-y[j]) - np.linalg.norm(x[i]-x[j])) / 1. for i,j in unchanged)
        out.extend((.12 * (y-x)).ravel())
        out.extend(0.7 * ((progress[a] - progress[b]) -
                   (progress_target[a]-progress_target[b])) for a,b in links)
        out.extend((cosine_angle(y,a['atoms'])-np.cos(np.deg2rad(a['degrees']))) / .5 for a in angles)
        for i in range(len(y)):
            for j in range(i):
                floor = .6 * (covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]])
                out.append(3 * max(0., floor - np.linalg.norm(y[i]-y[j])))
        return out
    initial = x + rng.normal(0., .015, x.shape)
    fit = least_squares(residual, initial.ravel(), max_nfev=200,
                        ftol=1e-5, xtol=1e-5, gtol=1e-5)
    if not fit.success or not np.isfinite(fit.x).all():
        raise ValueError(f'Seed geometry fit did not converge: {fit.message}')
    y = fit.x.reshape(x.shape)
    displacement = internal_direction(x, y-x)
    seed = x + amplitude * displacement
    direction = reaction_direction(seed, ij, delta * progress_target, unchanged)
    return seed, direction, dict(fraction=fraction, template_id=proposal['template_id'],
                            displacement_norm_A=amplitude,
                            unnormalized_fit_displacement_A=float(np.linalg.norm(y-x)),
                            progress_targets=progress_target.tolist(), angle_targets=angles,
                            angle_cosines_at_seed=[cosine_angle(seed,a['atoms']) for a in angles],
                            mode_policy='damped_internal_coordinate_tangent_at_seed',
                            template_source_graph=proposal.get('template_source_graph'),
                            transferred=proposal.get('transferred'),
                            origin=proposal.get('origin'), pattern_smarts=proposal.get('pattern_smarts'),
                            encounter_orientation=encounter,
                            intermolecular=proposal.get('intermolecular',False),
                            source_sink_couplings=links, geometric_fit_cost=float(fit.cost),
                            geometric_fit_nfev=fit.nfev, geometric_fit_optimality=float(fit.optimality),
                            predicted_graph=proposal['predicted_graph'],
                            arrows=proposal['arrows'] if strategy=='arrows' else None,
                            net_edits=edits)
