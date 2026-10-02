"""Replay saved attempt seeds under alternative TS / connection protocols (P1, P2).

Each saved attempt directory holds the exact seed geometry and initial
direction. A replay rebuilds that run's protocol and potential, then runs
search_connection from the identical seed under each requested variant.
The 'legacy' variant repeats the stored protocol and is compared with the
historical result.json (a reproducibility check). Nothing is registered into
a network; each replay is one billed attempt.

Sample once (deterministic), then run array tasks over the sample:
  replay_attempts.py sample --out sample.json [--p1 100 --p2 60 --seed 20261002]
  replay_attempts.py run --sample sample.json --variants legacy sella --task 0 --tasks 16 --out DIR
P1c (TS stage on large systems) samples saved seeds of campaign runs instead:
  replay_attempts.py sample-campaign CAMPAIGN_DIR [...] --out sample.json [--large 120 --small 120]
Seeds are deduplicated by content (T1x reactions that share a reactant repeat
the same seeds) and split by size (--large-atoms); within each stratum the
draw is round-robin over (start, strategy) after a seeded shuffle.
"""
import argparse
import dataclasses
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.potentials import load_potential  # noqa: E402
from mechbridge.reaction_network import CountedCalculator, SearchProtocol, aligned_rmsd, search_connection  # noqa: E402

VARIANTS = {
    'legacy': {},
    'sella': dict(ts_optimizer='dimer+sella'),
    'irc': dict(connection_protocol='irc'),
    'sella_irc': dict(ts_optimizer='dimer+sella', connection_protocol='irc'),
    # P3 fmax sensitivity around the candidate protocol (TS and endpoint gates move together).
    'sella_irc_f003': dict(ts_optimizer='dimer+sella', connection_protocol='irc', fmax=.003),
    'sella_irc_f010': dict(ts_optimizer='dimer+sella', connection_protocol='irc', fmax=.010),
    # P1c, applied on top of a run's stored (draft) protocol: earlier Dimer->Sella
    # hand-off with more Sella steps, near Sella-only, and a longer Dimer.
    'p1c_handoff03': dict(dimer_fmax=.3, sella_steps=300),
    'p1c_handoff10': dict(dimer_fmax=1., sella_steps=400),
    'p1c_long_dimer': dict(ts_steps=400),
}
P1C_VARIANTS = ('legacy', 'p1c_handoff03', 'p1c_handoff10', 'p1c_long_dimer')
VALIDATED = ('validated_descents', 'validated_core_descents')


def saved_attempts():
    rows = []
    for direction in sorted(ROOT.glob('reports/**/attempt_*/seed_direction.npy')):
        attempt = direction.parent
        result_path, network_path = attempt/'result.json', attempt.parent/'network.json'
        manifest_path = attempt.parent.parent.parent/'manifest.json'
        if not (result_path.exists() and network_path.exists() and manifest_path.exists()
                and (attempt/'seed.xyz').exists()):
            continue
        status = json.loads(result_path.read_text(encoding='utf-8')).get('status')
        if status in (None, 'seed_generation_failed'):
            continue
        rows.append(dict(attempt=attempt.relative_to(ROOT).as_posix(), family=attempt.parts[len(ROOT.parts)+1],
                         historical_status=status))
    return rows


def round_robin(rows, count, rng):
    """Spread the sample over run families instead of the largest campaign."""
    groups = {}
    for row in rows:
        groups.setdefault(row['family'], []).append(row)
    for members in groups.values():
        rng.shuffle(members)
    chosen, families = [], sorted(groups)
    while len(chosen) < count and any(groups.values()):
        for family in families:
            if groups[family] and len(chosen) < count:
                chosen.append(groups[family].pop())
    return chosen


def sample(args):
    rng = np.random.default_rng(args.seed)
    rows = saved_attempts()
    p1 = round_robin([dict(r) for r in rows], args.p1, rng)
    taken = {r['attempt'] for r in p1}
    p2 = round_robin([dict(r) for r in rows if r['historical_status'] in VALIDATED
                      and r['attempt'] not in taken], args.p2, rng)
    plan = dict(seed=args.seed, pool=len(rows), P1=p1, P2=p2,
                rule='round-robin over run families after a seeded shuffle; '
                     'P2 draws only historically validated attempts not already in P1')
    args.out.write_text(json.dumps(plan, indent=2), encoding='utf-8')
    print(json.dumps(dict(pool=len(rows), P1=len(p1), P2=len(p2)), indent=2))


def campaign_attempts(campaigns, large_atoms):
    """Saved seeds of campaign runs, deduplicated by seed content, with size strata."""
    rows, seen = [], set()
    for campaign in campaigns:
        for direction in sorted(campaign.resolve().glob('runs/*/*/*/*/*/attempt_*/seed_direction.npy')):
            attempt = direction.parent
            if not ((attempt/'seed.xyz').exists() and (attempt/'result.json').exists()):
                continue
            status = json.loads((attempt/'result.json').read_text(encoding='utf-8')).get('status')
            if status in (None, 'seed_generation_failed'):
                continue
            digest = hashlib.sha256((attempt/'seed.xyz').read_bytes() + direction.read_bytes()).hexdigest()
            if digest in seen:
                continue
            seen.add(digest)
            atoms = len(read(attempt/'seed.xyz'))
            start, strategy = attempt.parts[-3], attempt.parts[-2]
            rows.append(dict(attempt=str(attempt), family=f'{start}/{strategy}', historical_status=status,
                             atoms=atoms, stratum='large' if atoms >= large_atoms else 'small'))
    return rows


def sample_campaign(args):
    rng = np.random.default_rng(args.seed)
    rows = campaign_attempts(args.campaigns, args.large_atoms)
    plan = dict(seed=args.seed, campaigns=[str(c) for c in args.campaigns], large_atoms=args.large_atoms,
                pool={s: sum(r['stratum'] == s for r in rows) for s in ('large', 'small')},
                P1c_large=round_robin([r for r in rows if r['stratum'] == 'large'], args.large, rng),
                P1c_small=round_robin([r for r in rows if r['stratum'] == 'small'], args.small, rng),
                rule='content-deduplicated saved seeds; round-robin over start/strategy after a seeded shuffle')
    args.out.write_text(json.dumps(plan, indent=2), encoding='utf-8')
    print(json.dumps(dict(pool=plan['pool'], large=len(plan['P1c_large']), small=len(plan['P1c_small'])),
                     indent=2))


def replay(entry, variant, backends, outdir, threads):
    attempt = ROOT/entry['attempt']
    network = json.loads((attempt.parent/'network.json').read_text(encoding='utf-8'))
    manifest = json.loads((attempt.parent.parent.parent/'manifest.json').read_text(encoding='utf-8'))
    protocol = dataclasses.replace(SearchProtocol(**network['protocol']), **VARIANTS[variant])
    potential = manifest['model']
    if potential not in backends:
        backends[potential] = load_potential(potential, ROOT, threads=threads)[0]
    seed = read(attempt/'seed.xyz')
    if (attempt/'seed_positions_A.npy').exists():     # seed.xyz keeps 8 decimals only
        seed.positions = np.load(attempt/'seed_positions_A.npy')
    direction = np.load(attempt/'seed_direction.npy')
    start = network['start']
    # As in explore(): the validator also sets the backend's charge/multiplicity.
    validator = getattr(backends[potential], 'validate_system', None)
    if validator is not None:
        validator(seed.numbers, start['charge'], start['multiplicity'])
    calculator = CountedCalculator(backends[potential], protocol.evaluations_per_attempt)
    result = search_connection(seed, direction, calculator, outdir, protocol, start['charge'])
    historical = json.loads((attempt/'result.json').read_text(encoding='utf-8'))
    row = dict(attempt=entry['attempt'], variant=variant, potential=potential, status=result['status'],
               evaluations=result['evaluations'], seconds=result['seconds'],
               historical_status=historical['status'], historical_evaluations=historical.get('evaluations'),
               ts_energy_eV=(result.get('ts') or {}).get('energy_eV'),
               imaginary_count=(result.get('ts') or {}).get('imaginary_count'),
               endpoint_graphs=[e.get('graph_smiles') for e in result['endpoints']],
               is_IRC=result['is_IRC'], ts_optimization=result.get('ts_optimization'))
    if variant == 'legacy':
        row['historical_match'] = (result['status'] == historical['status'] and
                                   result['evaluations'] == historical.get('evaluations'))
    if 'ts_positions_A' in result and 'ts_positions_A' in historical:
        row['ts_rmsd_to_historical_A'] = aligned_rmsd(np.asarray(result['ts_positions_A']),
                                                      np.asarray(historical['ts_positions_A']))
    return row


def run(args):
    plan = json.loads(args.sample.read_text(encoding='utf-8'))
    allowed = {('P1', 'legacy'), ('P1', 'sella'), ('P1', 'sella_irc'),
               ('P2', 'legacy'), ('P2', 'irc'), ('P2', 'sella_irc'),
               ('P2', 'sella_irc_f003'), ('P2', 'sella_irc_f010'),
               *((g, v) for g in ('P1c_large', 'P1c_small') for v in P1C_VARIANTS)}
    jobs = [(group, entry, variant) for group in ('P1', 'P2', 'P1c_large', 'P1c_small')
            for entry in plan.get(group, []) for variant in args.variants if (group, variant) in allowed]
    mine = jobs[args.task::args.tasks]
    args.out.mkdir(parents=True, exist_ok=True)
    backends, rows = {}, []
    for group, entry, variant in mine:
        key = entry['attempt'].replace('/', '__')
        outdir = args.out/'attempts'/group/key/variant
        row = replay(entry, variant, backends, outdir, args.threads)
        row['group'] = group
        rows.append(row)
        print(json.dumps(row), flush=True)
    (args.out/f'task_{args.task:03d}.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    s = sub.add_parser('sample')
    s.add_argument('--out', type=Path, required=True)
    s.add_argument('--p1', type=int, default=100)
    s.add_argument('--p2', type=int, default=60)
    s.add_argument('--seed', type=int, default=20261002)
    c = sub.add_parser('sample-campaign')
    c.add_argument('campaigns', type=Path, nargs='+')
    c.add_argument('--out', type=Path, required=True)
    c.add_argument('--large', type=int, default=120)
    c.add_argument('--small', type=int, default=120)
    c.add_argument('--large-atoms', type=int, default=25)
    c.add_argument('--seed', type=int, default=20261002)
    r = sub.add_parser('run')
    r.add_argument('--sample', type=Path, required=True)
    r.add_argument('--variants', nargs='+', choices=sorted(VARIANTS), required=True)
    r.add_argument('--task', type=int, default=0)
    r.add_argument('--tasks', type=int, default=1)
    r.add_argument('--threads', type=int, default=2)
    r.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command in ('sample', 'sample-campaign'):
        if args.out.exists():
            raise FileExistsError(args.out)
        (sample if args.command == 'sample' else sample_campaign)(args)
    else:
        if not 0 <= args.task < args.tasks:
            parser.error('--task must be in [0, --tasks)')
        run(args)


if __name__ == '__main__':
    main()
