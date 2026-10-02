"""Symbolic-library coverage of benchmark reactions (an upper bound for symbolic arms).

For each admitted start: the number of symbolic proposals the library makes for
the reactant graph, the number of distinct predicted products, and whether any
predicted product equals the held-out reference product (stereo-free key). The
reference is read only to score coverage, never to make proposals. A reaction
the library cannot propose can still be found by the geometry arm or as an
unexpected endpoint, so this bounds only "proposed and then found".
Usage: symbolic_coverage.py --starts F --references F [--ids F] --out FILE
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol  # noqa: E402
from mechbridge.reference_matching import stereo_free, wilson  # noqa: E402
from mechbridge.symbolic_library import ArrowLibrary  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--ids', type=Path, help='JSON with representable_ids (restricts the set)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    library = ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    references = {r['id']: r for r in map(json.loads, args.references.read_text(encoding='utf-8').splitlines())}
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data['representable_ids'] if isinstance(data, dict) else data)
    rows = []
    for start in map(json.loads, args.starts.read_text(encoding='utf-8').splitlines()):
        if keep is not None and start['id'] not in keep:
            continue
        reference = references[start['id']]
        mol = geometry_mol(start['atomic_numbers'], np.asarray(start['positions_A']), start['charge'])
        try:
            proposals = library.propose(mol)
            error = None
        except ValueError as exc:
            proposals, error = [], f'{type(exc).__name__}: {exc}'
        products = sorted({stereo_free(p['predicted_graph']) for p in proposals})
        origins = sorted({str(p.get('origin')) for p in proposals})
        rows.append(dict(id=start['id'], proposals=len(proposals), predicted_products=len(products),
                         reference_product_proposed=reference['product_key'] in products,
                         origins=origins, error=error))
    n = len(rows)
    covered = sum(r['reference_product_proposed'] for r in rows)
    with_any = sum(r['proposals'] > 0 for r in rows)
    summary = dict(starts=str(args.starts), reactions=n, any_proposal=with_any,
                   any_proposal_wilson95=wilson(with_any, n), reference_product_proposed=covered,
                   reference_product_proposed_wilson95=wilson(covered, n),
                   median_proposals=float(np.median([r['proposals'] for r in rows])) if rows else None,
                   rows=rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
