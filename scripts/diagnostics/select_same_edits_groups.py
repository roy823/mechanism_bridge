"""Pre-registered Fig. 5b selection of Coley "drawing" groups (same edits, different arrows).

From a same_edits_groups.py output, Coley groups are split into those with at
least one arrow set from a resonance form and those without; --per-kind groups
of each kind are drawn with default_rng(seed) and written as
{start_id: [edits, ...]} for same_edits_arrows.py --edit-keys.
Usage: select_same_edits_groups.py GROUPS.json --out FILE [--per-kind 10 --seed 20261002 --prefix coley_]
"""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('groups', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--per-kind', type=int, default=10)
    parser.add_argument('--seed', type=int, default=20261002)
    parser.add_argument('--prefix', default='coley_')
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    rows = json.loads(args.groups.read_text(encoding='utf-8'))
    pools = dict(resonance=[], direct=[])
    for row in rows:
        if not row['start'].startswith(args.prefix):
            continue
        for g in row.get('groups', []):
            kind = 'resonance' if any('resonance' in (a.get('origin') or '') for a in g['arrow_sets']) else 'direct'
            pools[kind].append((row['start'], g['edits']))
    rng = np.random.default_rng(args.seed)
    selection, chosen = {}, {}
    for kind in ('resonance', 'direct'):          # fixed order of draws
        pool = pools[kind]
        picks = rng.choice(len(pool), size=min(args.per_kind, len(pool)), replace=False)
        chosen[kind] = [pool[k][0] for k in sorted(picks)]
        for k in sorted(picks):
            selection.setdefault(pool[k][0], []).append(pool[k][1])
    args.out.write_text(json.dumps(selection, indent=2), encoding='utf-8')
    print(json.dumps(dict(pool_sizes={k: len(v) for k, v in pools.items()}, chosen=chosen,
                          starts=len(selection)), indent=2))


if __name__ == '__main__':
    main()
