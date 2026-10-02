"""Attach author TS geometries to the Coley [3+2] references for the representability screen.

The autodE profiles (full_data_profiles.tar.gz) hold each optimized TS in the
author's atom order, which differs from our starts (heavy atoms of both
reactants first, then hydrogens). In a [3+2] cycloaddition every reactant bond
survives at the TS and only the two forming bonds (about 2.2 A) are new, so the
TS connectivity from a distance criterion (RDKit DetermineConnectivity, no bond
orders) equals the reactant connectivity. The covalent-radius factor is reduced
from 1.3 until that graph is isomorphic (element-labelled) to our reactant bonds;
the factor is recorded. Any isomorphism is acceptable because screening and
scoring allow reactant automorphisms. The TS is used only by the post-search screen, never as
a search input. References that cannot be matched are kept with the reason and
admission='no_reference_ts', so the screen skips them.
Usage: prepare_coley_ts_references.py --profiles full_data_profiles.tar.gz --out FILE
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np
from ase.io import read
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.reference_matching import mol_bonds  # noqa: E402

FACTORS = (1.3, 1.25, 1.2, 1.15, 1.1)


def labelled_graph(numbers, bonds):
    graph = nx.Graph()
    graph.add_nodes_from((i, dict(z=int(z))) for i, z in enumerate(numbers))
    graph.add_edges_from(tuple(map(int, b[:2])) for b in bonds)
    return graph


def author_ts(members, rxn_id):
    """The optimized TS (ASE Atoms) of one reaction, in the author's atom order."""
    names = [n for n in members if n.startswith(f'full_dataset_profiles/{rxn_id}/TS_')
             and n.endswith('.xyz') and 'imag_mode' not in n]
    if len(names) != 1:
        raise ValueError(f'{len(names)} TS files')
    return read(io.StringIO(members[names[0]]), format='xyz')


def connectivity(numbers, positions, factor):
    table = Chem.GetPeriodicTable()
    block = f'{len(numbers)}\n\n' + '\n'.join(f'{table.GetElementSymbol(int(z))} {x:.10f} {y:.10f} {w:.10f}'
                                             for z, (x, y, w) in zip(numbers, positions))
    mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(mol, useHueckel=False, covFactor=factor)
    return mol_bonds(mol)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', type=Path,
                        default=ROOT/'data/raw/coley_dipolar/full_data_profiles.tar.gz')
    parser.add_argument('--references', type=Path, default=ROOT/'data/processed/coley_benchmark_references.jsonl')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    RDLogger.DisableLog('rdApp.*')
    references = [json.loads(l) for l in args.references.read_text(encoding='utf-8').splitlines()]
    wanted = {r['rxn_id'] for r in references}
    members = {}
    with tarfile.open(args.profiles, 'r:gz') as archive:
        for info in archive:
            parts = info.name.split('/')
            if info.isfile() and len(parts) == 3 and parts[1] in wanted and info.name.endswith('.xyz'):
                members[info.name] = archive.extractfile(info).read().decode()
    rows, failures = [], {}
    for reference in references:
        out = dict(reference)
        try:
            ts = author_ts(members, reference['rxn_id'])
            ours = labelled_graph(reference['atomic_numbers'], reference['reactant_bonds'])
            mapping = factor = None
            for factor in FACTORS:
                theirs = labelled_graph(ts.numbers, connectivity(ts.numbers, ts.positions, factor))
                matcher = isomorphism.GraphMatcher(theirs, ours, node_match=lambda a, b: a['z'] == b['z'])
                mapping = next(matcher.isomorphisms_iter(), None)
                if mapping is not None:
                    break
            if mapping is None:
                raise ValueError('TS connectivity never matches the reference reactant')
            order = np.empty(len(ts), dtype=int)
            for author, mine in mapping.items():
                order[mine] = author
            out['ts_positions_A'] = ts.positions[order].tolist()
            out['ts_source'] = ('Coley/Stuyver autodE profile TS (optimized); author order mapped by an '
                                f'isomorphism of the TS connectivity at covalent factor {factor}')
        except (ValueError, RuntimeError, KeyError) as exc:
            # This file only feeds the screen, which reads admitted references; scoring
            # keeps the original references file, where the reaction stays admitted.
            out.update(ts_positions_A=None, ts_source=f'unavailable: {exc}', admission='no_reference_ts')
            failures[reference['id']] = str(exc)
        rows.append(out)
    args.out.write_text(''.join(json.dumps(r) + '\n' for r in rows), encoding='utf-8')
    print(json.dumps(dict(references=len(rows), with_ts=sum(r['ts_positions_A'] is not None for r in rows),
                          failures=failures,
                          profiles_sha256=hashlib.sha256(args.profiles.read_bytes()).hexdigest(),
                          out_sha256=hashlib.sha256(args.out.read_bytes()).hexdigest()), indent=2))


if __name__ == '__main__':
    main()
