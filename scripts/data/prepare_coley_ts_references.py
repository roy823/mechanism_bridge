"""Attach author TS geometries to the Coley [3+2] references for the representability screen.

The autodE profiles (full_data_profiles.tar.gz) store each TS with the atoms of
reactant r0 followed by r1. Our starts follow the mapped reaction SMILES, so the
author order is matched to ours by an element-labelled isomorphism between the
connectivity perceived from the author's r0 and r1 geometries and our reactant
bonds; any isomorphism is acceptable because screening and scoring allow
reactant automorphisms. The TS is used only by the post-search screen, never as
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
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol  # noqa: E402
from mechbridge.reference_matching import mol_bonds  # noqa: E402


def labelled_graph(numbers, bonds):
    graph = nx.Graph()
    graph.add_nodes_from((i, dict(z=int(z))) for i, z in enumerate(numbers))
    graph.add_edges_from(tuple(map(int, b[:2])) for b in bonds)
    return graph


def author_reaction(members, rxn_id):
    """(numbers, reactant bonds, TS positions) in the author's r0+r1 order."""
    def xyz(prefix):
        names = [n for n in members if n.startswith(f'full_dataset_profiles/{rxn_id}/{prefix}')
                 and n.endswith('.xyz') and 'imag_mode' not in n]
        if len(names) != 1:
            raise ValueError(f'{len(names)} files for {prefix}')
        return read(io.StringIO(members[names[0]]), format='xyz')
    r0, r1, ts = xyz('r0_'), xyz('r1_'), xyz('TS_')
    numbers = list(r0.numbers) + list(r1.numbers)
    if list(ts.numbers) != numbers:
        raise ValueError('TS atoms are not r0 followed by r1')
    bonds = set(mol_bonds(geometry_mol(r0.numbers, r0.positions, 0)))
    bonds |= {(i + len(r0), j + len(r0)) for i, j in mol_bonds(geometry_mol(r1.numbers, r1.positions, 0))}
    return numbers, bonds, ts.positions


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
            numbers, bonds, ts = author_reaction(members, reference['rxn_id'])
            ours = labelled_graph(reference['atomic_numbers'], reference['reactant_bonds'])
            theirs = labelled_graph(numbers, bonds)
            matcher = isomorphism.GraphMatcher(theirs, ours, node_match=lambda a, b: a['z'] == b['z'])
            mapping = next(matcher.isomorphisms_iter(), None)
            if mapping is None:
                raise ValueError('author reactant connectivity differs from the reference reactant')
            order = np.empty(len(numbers), dtype=int)
            for author, mine in mapping.items():
                order[mine] = author
            out['ts_positions_A'] = ts[order].tolist()
            out['ts_source'] = 'Coley/Stuyver autodE profile TS, author order mapped by isomorphism'
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
