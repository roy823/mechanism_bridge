"""Group Transition1x reference reactions that share one reactant geometry.

Transition1x enumerates several reactions from the same optimized reactant, so a
blind search from that reactant is the same run for all of them. Reactions in
the given ID list are grouped by identical atomic numbers and positions; each
group is represented by its first ID in sorted order. A campaign then runs the
representatives only ('ids'), and score_reference_recovery.py --groups scores
each representative run against every reaction of its group.
Usage: group_t1x_reactants.py --ids REPRESENTABLE.json --out GROUPS.json [--starts FILE]
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, default=ROOT/'data/processed/t1x_test_starts.jsonl')
    parser.add_argument('--ids', type=Path, required=True, help='JSON with "representable_ids" or "ids"')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    data = json.loads(args.ids.read_text(encoding='utf-8'))
    keep = set(data.get('representable_ids') or data['ids'])
    groups = {}
    for line in args.starts.read_text(encoding='utf-8').splitlines():
        start = json.loads(line)
        if start['id'] not in keep:
            continue
        key = hashlib.sha256(json.dumps([start['atomic_numbers'], start['positions_A']]).encode()).hexdigest()
        groups.setdefault(key, []).append(start['id'])
    if sum(map(len, groups.values())) != len(keep):
        raise ValueError('some IDs have no start record')
    members = {min(ids): sorted(ids) for ids in groups.values()}
    representatives = sorted(members)
    sizes = {}
    for ids in members.values():
        sizes[len(ids)] = sizes.get(len(ids), 0) + 1
    result = dict(purpose=__doc__.splitlines()[0], starts=str(args.starts.name),
                  starts_sha256=hashlib.sha256(args.starts.read_bytes()).hexdigest(),
                  ids_file=str(args.ids.name), ids_sha256=hashlib.sha256(args.ids.read_bytes()).hexdigest(),
                  reactions=len(keep), reactants=len(representatives), group_sizes=dict(sorted(sizes.items())),
                  ids=representatives, groups=members)
    args.out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('ids', 'groups')}, indent=2))


if __name__ == '__main__':
    main()
