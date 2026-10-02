"""Write a campaign manifest: one row per (start, strategy, seed) MLIP exploration.

The manifest is the pre-registered run list for a SLURM array: rows are in a
fixed order (start file, start id, strategy, seed), task ids start at 1, and
the manifest records the protocol JSON hash so that every array task can
refuse a different protocol. Starts may be restricted to an ID list (e.g. the
frozen representable set).

Usage: build_campaign_manifest.py --name NAME --starts FILE [--ids FILE] \
         --strategies geometry arrows ... --seeds 17 29 --protocol P.json --potential aimnet2-rxn --out DIR
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STRATEGIES = ('geometry', 'center_random', 'bond_edits', 'arrows', 'hybrid',
              'pygsm_edits', 'pygsm_b2f2', 'pygsm_oracle')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True)
    parser.add_argument('--starts', type=Path, nargs='+', required=True)
    parser.add_argument('--ids', type=Path, help='JSON with "representable_ids" or a list of IDs')
    parser.add_argument('--limit', type=int, help='Keep the first N starts per file (pilots)')
    parser.add_argument('--strategies', nargs='+', choices=STRATEGIES, required=True)
    parser.add_argument('--seeds', type=int, nargs='+', required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--potential', default='aimnet2-rxn')
    parser.add_argument('--extra-args', default='',
                        help='Extra runner arguments recorded in campaign.json, e.g. "--oracle-references F"')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    manifest = args.out/'manifest.tsv'
    if manifest.exists():
        raise FileExistsError(manifest)
    keep = None
    if args.ids:
        data = json.loads(args.ids.read_text(encoding='utf-8'))
        keep = set(data if isinstance(data, list) else data.get('representable_ids') or data['ids'])
    rows = []
    for starts in args.starts:
        ids = [json.loads(line)['id'] for line in starts.read_text(encoding='utf-8').splitlines()]
        ids = [i for i in ids if keep is None or i in keep][:args.limit]
        for start in ids:
            for strategy in args.strategies:
                for seed in args.seeds:
                    rows.append(dict(start_file=starts.resolve().relative_to(ROOT).as_posix(),
                                     start_id=start, strategy=strategy, seed=seed,
                                     potential=args.potential,
                                     outdir=f'{start}/{strategy}/s{seed}'))
    args.out.mkdir(parents=True, exist_ok=True)
    with manifest.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['task_id', *rows[0]], delimiter='\t')
        writer.writeheader()
        for task_id, row in enumerate(rows, start=1):
            writer.writerow(dict(task_id=task_id, **row))
    protocol_sha = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    meta = dict(name=args.name, tasks=len(rows), strategies=args.strategies, seeds=args.seeds,
                starts=[s.resolve().relative_to(ROOT).as_posix() for s in args.starts],
                ids_filter=None if args.ids is None else str(args.ids), limit=args.limit,
                protocol=str(args.protocol), protocol_sha256=protocol_sha, potential=args.potential,
                extra_args=args.extra_args,
                manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest())
    (args.out/'campaign.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print(json.dumps(meta, indent=2))


if __name__ == '__main__':
    main()
