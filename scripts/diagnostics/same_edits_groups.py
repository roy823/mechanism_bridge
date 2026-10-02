"""Fig. 5b candidates: proposals that share net bond edits but differ in their arrows.

For each start, the FP-JCTC-1 symbolic library (published arrows and analyst
grammar, up to eight resonance forms) proposes actions at the start graph.
Proposals are grouped by their net bond edits; groups with at least two
distinct arrow sets are the "same edits, different arrows" cases of Fig. 5b.
No PES evaluation is made.
Usage: same_edits_groups.py --starts FILE [FILE ...] [--ids F.json] --out FILE
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.symbolic_library import ArrowLibrary, ResonanceAwareLibrary  # noqa: E402


def signature(items):
    return json.dumps(items, sort_keys=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, nargs='+', required=True)
    parser.add_argument('--ids', type=Path, help='JSON list or {"ids": [...]}')
    parser.add_argument('--resonance-forms', type=int, default=8)
    parser.add_argument('--limit', type=int, default=64, help='Proposals requested per start')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data if isinstance(data, list) else data['ids'])
    library = ResonanceAwareLibrary(ArrowLibrary(ROOT/'data/raw/synepd/polar.json'), args.resonance_forms)
    rows = []
    for path in args.starts:
        for line in path.read_text(encoding='utf-8').splitlines():
            start = json.loads(line)
            if keep is not None and start['id'] not in keep:
                continue
            try:
                mol = geometry_mol(start['atomic_numbers'], np.asarray(start['positions_A']), start.get('charge', 0))
                proposals = library.propose(mol, limit=args.limit)
            except (ValueError, RuntimeError) as exc:
                rows.append(dict(start=start['id'], error=str(exc)))
                continue
            groups = {}
            for p in proposals:
                groups.setdefault(signature(p['edits']), {}).setdefault(signature(p['arrows']), p)
            multi = [dict(edits=json.loads(k), predicted_graph=next(iter(v.values()))['predicted_graph'],
                          arrow_sets=[dict(template_id=p['template_id'], origin=p.get('origin'),
                                           arrows=p['arrows']) for p in v.values()])
                     for k, v in groups.items() if len(v) >= 2]
            rows.append(dict(start=start['id'], graph=graph_smiles(mol), proposals=len(proposals),
                             edit_groups=len(groups), multi_arrow_groups=len(multi), groups=multi))
            print(json.dumps(dict(start=start['id'], proposals=len(proposals), edit_groups=len(groups),
                                  multi_arrow_groups=len(multi))), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
