"""Arrow-only seed features (seed_features='arrow_features_v1').

Everything here needs proposal['arrows']; priors derivable from the net bond
edits stay in search_seeds and are shared with the bond_edits control. The
features are: (a) lone-pair direction of an electron-pair donor, (c) the
directed push-pull order of the arrow chain, used both as an asynchronous
progress shift and as generalized source/sink coupling, and (f) donor
lone-pair alignment in the rigid encounter score. The bond_edits control gets
the same progress-shift axis with a random bond order (no chemical order).
"""
import numpy as np

VALENCE_ELECTRONS = {1: 1, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 15: 5, 16: 6, 17: 7, 35: 7}
CHAIN_SHIFTS = (0., -.3, .3)        # synchronous, push-end first, pull-end first


def lone_pair_count(mol, i):
    """Lone pairs on atom i from valence electrons, explicit bond orders and formal charge."""
    atom = mol.GetAtomWithIdx(int(i))
    bonded = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
    free = VALENCE_ELECTRONS.get(atom.GetAtomicNum(), 0) - atom.GetFormalCharge() - bonded
    return max(0, int(round(free)) // 2)


def lone_pair_directions(y, mol, i):
    """VSEPR lone-pair unit vectors on atom i at geometry y.

    [] means no directional prior (no lone pair or no neighbours); None means a
    single neighbour, where only a cone angle is defined.
    """
    neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(int(i)).GetNeighbors()]
    pairs = lone_pair_count(mol, i)
    if pairs == 0 or not neighbours:
        return []
    if len(neighbours) == 1:
        return None
    units = [(y[j]-y[i])/max(np.linalg.norm(y[j]-y[i]), 1e-10) for j in neighbours]
    s = -np.sum(units, axis=0)
    bisector = s/max(np.linalg.norm(s), 1e-10)
    if pairs == 1 or len(neighbours) >= 3:
        return [bisector]
    normal = np.cross(units[0], units[1])
    normal /= max(np.linalg.norm(normal), 1e-10)
    half = np.deg2rad(54.75)                        # two sp3 lone pairs (water, ether)
    return [np.cos(half)*bisector + np.sin(half)*normal, np.cos(half)*bisector - np.sin(half)*normal]


def lone_pair_terms(arrows, mol, x):
    """(a) Donor d gives a lone pair into the new bond d-a; choose the closest lone pair at x."""
    terms = []
    for arrow in arrows:
        if len(arrow['source']) != 1 or len(arrow['sink']) != 2:
            continue
        donor = int(arrow['source'][0])
        if donor not in arrow['sink']:
            continue
        acceptor = int(next(k for k in arrow['sink'] if k != donor))
        directions = lone_pair_directions(x, mol, donor)
        if directions == []:
            continue
        if directions is None:
            neighbour = next(iter(mol.GetAtomWithIdx(donor).GetNeighbors())).GetIdx()
            steric = 1 + lone_pair_count(mol, donor)
            terms.append(dict(kind='lp_cone', donor=donor, acceptor=acceptor, neighbor=int(neighbour),
                              degrees={4: 109.5, 3: 120., 2: 180.}.get(steric, 109.5)))
        else:
            v = (x[acceptor]-x[donor])/max(np.linalg.norm(x[acceptor]-x[donor]), 1e-10)
            pick = int(np.argmax([float(np.dot(l, v)) for l in directions]))
            terms.append(dict(kind='lp_axis', donor=donor, acceptor=acceptor, pick=pick))
    return terms


def _cosine(y, i, j, k):
    u, v = y[i]-y[j], y[k]-y[j]
    return float(np.dot(u, v)/max(np.linalg.norm(u)*np.linalg.norm(v), 1e-12))


def lone_pair_residuals(y, mol, terms, sample):
    """Hinge residuals: penalize only misalignment beyond the sampled tolerance."""
    tolerance = np.cos(np.deg2rad((10., 20., 30.)[sample % 3]))
    out = []
    for t in terms:
        d, a = t['donor'], t['acceptor']
        if t['kind'] == 'lp_axis':
            directions = lone_pair_directions(y, mol, d)
            v = (y[a]-y[d])/max(np.linalg.norm(y[a]-y[d]), 1e-10)
            out.append(max(0., tolerance - float(np.dot(directions[t['pick']], v)))/.3)
        else:
            out.append((_cosine(y, t['neighbor'], d, a) - np.cos(np.deg2rad(t['degrees'])))/.3)
    return out


def arrow_chain(arrows):
    """(c) Directed push-pull order: per-arrow BFS rank, successors and a cyclic flag.

    The head of an arrow is the sink atom not in its source (or the sink atom
    itself); arrow k precedes m when its head starts m. A cycle without a
    first arrow starts from list position 0 (weak information, flagged).
    """
    def head(arrow):
        atoms = [x for x in arrow['sink'] if x not in arrow['source']]
        return set(atoms or arrow['sink'])
    n = len(arrows)
    succ = {k: [m for m in range(n) if m != k and head(arrows[k]) & set(arrows[m]['source'])]
            for k in range(n)}
    indegree = {m: sum(m in v for v in succ.values()) for m in range(n)}
    starts = [k for k in range(n) if indegree[k] == 0]
    cyclic = not starts
    rank, queue = {k: 0 for k in (starts or [0])}, list(starts or [0])
    while queue:
        k = queue.pop(0)
        for m in succ[k]:
            if m not in rank:
                rank[m] = rank[k] + 1
                queue.append(m)
    top = max(rank.values()) if rank else 0
    for k in range(n):
        rank.setdefault(k, top + 1)
    return rank, succ, cyclic


def touched_bonds(arrow):
    return [tuple(sorted(int(a) for a in side)) for side in (arrow['source'], arrow['sink']) if len(side) == 2]


def chain_progress_shift(edits, arrows, sample, random_seed):
    """Per-edit progress shift lambda*rho; arrows=None gives the matched random-order control."""
    lam = CHAIN_SHIFTS[(sample + 1) % 3]           # Latin pairing with the sample%3 axes
    keys = [tuple(sorted(int(a) for a in e['atoms'])) for e in edits]
    # Separate stream, drawn by both arms, so the main seed RNG stays matched.
    permutation = np.random.default_rng([int(random_seed), 7]).permutation(len(keys)).astype(float)
    if arrows is None:
        order = permutation
    else:
        rank, _, _ = arrow_chain(arrows)
        acc = {k: [] for k in keys}
        for i, arrow in enumerate(arrows):
            for bond in touched_bonds(arrow):
                if bond in acc:
                    acc[bond].append(rank[i])
        fallback = float(np.mean(list(rank.values()))) if rank else 0.
        order = np.array([np.mean(acc[k]) if acc[k] else fallback for k in keys], dtype=float)
    span = order.max() - order.min() if len(order) else 0.
    if span < 1e-9:
        return np.zeros(len(keys)), lam
    return lam*((order - order.min())/span - .5), lam


def generalized_links(arrows, index):
    """Couple every bond touched by arrow k with every bond touched by its successor m."""
    _, succ, _ = arrow_chain(arrows)
    links = set()
    for k, followers in succ.items():
        for m in followers:
            for bk in touched_bonds(arrows[k]):
                for bm in touched_bonds(arrows[m]):
                    if bk != bm and bk in index and bm in index:
                        links.add((index[bk], index[bm]))
    return sorted(links)


def lone_pair_alignment_score(y, mol, arrows, weight=2.):
    """(f) Rigid-encounter addend: donor lone pair pointing at its acceptor (best lone pair)."""
    score = 0.
    for arrow in arrows:
        if len(arrow['source']) != 1 or len(arrow['sink']) != 2:
            continue
        donor = int(arrow['source'][0])
        if donor not in arrow['sink']:
            continue
        acceptor = int(next(k for k in arrow['sink'] if k != donor))
        directions = lone_pair_directions(y, mol, donor)
        if not directions:
            continue
        v = (y[acceptor]-y[donor])/max(np.linalg.norm(y[acceptor]-y[donor]), 1e-10)
        best = max(float(np.dot(l, v)) for l in directions)
        score += weight*max(0., .8 - best)**2
    return score
