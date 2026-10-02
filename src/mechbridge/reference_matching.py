"""Compare observed endpoints with held-out reference reactions (benchmark scoring only).

Atom order is shared between a benchmark start and its reference, so
connectivities can be compared index by index after an automorphism of the
reactant graph (equivalent hydrogens, symmetric fragments). Nothing here is
used during search.
"""
import math

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem

MAX_AUTOMORPHISMS = 5000


def stereo_free(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Cannot parse graph SMILES: {smiles}')
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def bond_set(bonds):
    """Unordered atom pairs from [[i, j], ...] or [[i, j, order], ...]."""
    return frozenset(tuple(sorted(map(int, b[:2]))) for b in bonds)


def mol_bonds(mol):
    return frozenset(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds())


def automorphisms(numbers, bonds, limit=MAX_AUTOMORPHISMS):
    """Element-preserving automorphisms of a connectivity graph (truncated at limit)."""
    graph = nx.Graph()
    graph.add_nodes_from((i, dict(z=int(z))) for i, z in enumerate(numbers))
    graph.add_edges_from(tuple(map(int, b[:2])) for b in bonds)
    matcher = GraphMatcher(graph, graph, node_match=lambda a, b: a['z'] == b['z'])
    maps = []
    for mapping in matcher.isomorphisms_iter():
        maps.append(mapping)
        if len(maps) >= limit:
            break
    return maps


def permute(bonds, mapping):
    return frozenset(tuple(sorted((mapping[i], mapping[j]))) for i, j in bonds)


def edge_matches(bonds_u, bonds_v, react, prod, maps):
    """One reactant automorphism must map the endpoints onto (reactant, product)."""
    return any((permute(bonds_u, m) == react and permute(bonds_v, m) == prod) or
               (permute(bonds_v, m) == react and permute(bonds_u, m) == prod) for m in maps)


def wilson(successes, n, z=1.96):
    """Wilson score interval for a binomial proportion; None when n == 0."""
    if not n:
        return None
    p = successes/n
    center = (p + z*z/(2*n))/(1 + z*z/n)
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1 + z*z/n)
    return [center - half, center + half]
