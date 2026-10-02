"""Fig. 5b: same net bond edits, different arrows (G doc 4.2-4.3).

For one start: relax the root on the MLIP once, propose at the root with the
FP-JCTC-1 library, and keep the groups of proposals that share net bond edits
but differ in their arrows. Every arrow set of a group gets samples x seeds
seeds (strategy 'arrows'); bond_edits gets samples x seeds x (number of arrow
sets) seeds with further random seeds, so both arms spend the same number of
attempts. The attempt RNG is derived exactly as in explore(), so the arrow
sets of a group share their random numbers and differ only in the arrows.
Each seed runs search_connection under the protocol (one billed attempt).
TS identity (G doc 4.3) is evaluated by summarize_same_edits.py.
Usage: same_edits_arrows.py --starts FILE --start-id ID --protocol-json P --outdir DIR
         [--seeds 17 29 43 59] [--samples 0 1 2] [--max-groups N]
"""
import argparse
import dataclasses
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts/exploration'))
from mechbridge.event_graph import geometry_mol, graph_smiles  # noqa: E402
from mechbridge.potentials import load_potential  # noqa: E402
from mechbridge.provenance import runtime_provenance  # noqa: E402
from mechbridge.reaction_network import (CountedCalculator, initialize_root, is_recoverable_failure,  # noqa: E402
                                         search_connection)
from mechbridge.search_seeds import make_seed  # noqa: E402
from mechbridge.symbolic_library import ArrowLibrary, ResonanceAwareLibrary  # noqa: E402
from run_network_exploration import protocol_from_json  # noqa: E402


def key(items):
    return json.dumps(items, sort_keys=True)


def seed_rng(graph, positions, edits, seed, variant):
    """Attempt RNG seed exactly as reaction_network.explore derives it."""
    identity = json.dumps(dict(graph=graph, positions=np.round(positions, 6).tolist(), edits=edits,
                               seed=seed, variant=variant), sort_keys=True)
    return int.from_bytes(hashlib.sha256(identity.encode()).digest()[:4], 'little')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--starts', type=Path, required=True)
    parser.add_argument('--start-id', required=True)
    parser.add_argument('--protocol-json', type=Path, required=True)
    parser.add_argument('--potential', default='aimnet2-rxn')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--seeds', type=int, nargs='+', default=[17, 29, 43, 59])
    parser.add_argument('--samples', type=int, nargs='+', default=[0, 1, 2])
    parser.add_argument('--max-groups', type=int, help='Keep the first N groups (proposal order)')
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(args.outdir)
    args.outdir.mkdir(parents=True)
    RDLogger.DisableLog('rdApp.*')
    start = next(s for s in map(json.loads, args.starts.read_text(encoding='utf-8').splitlines())
                 if s['id'] == args.start_id)
    protocol = protocol_from_json(args.protocol_json, args.seeds[0])
    backend, potential = load_potential(args.potential, ROOT, threads=args.threads)
    validator = getattr(backend, 'validate_system', None)
    if validator is not None:
        validator(start['atomic_numbers'], start['charge'], start['multiplicity'])
    calculator = CountedCalculator(backend, protocol.total_evaluations)
    root, _, init_evaluations = initialize_root(start, calculator, args.outdir, protocol)
    if not root['full_system_minimum']:
        raise SystemExit('root minimum unresolved')
    numbers, positions = start['atomic_numbers'], np.asarray(root['positions_A'])
    mol = geometry_mol(numbers, positions, start['charge'])
    graph = graph_smiles(mol)
    library = ResonanceAwareLibrary(ArrowLibrary(ROOT/'data/raw/synepd/polar.json'),
                                    protocol.proposal_resonance_forms)
    groups = {}
    for proposal in library.propose(mol, limit=64):
        groups.setdefault(key(proposal['edits']), {}).setdefault(key(proposal['arrows']), proposal)
    groups = [list(v.values()) for v in groups.values() if len(v) >= 2][:args.max_groups]
    rows = []
    def attempt(group, arm, proposal, sample, seed, strategy):
        rng = seed_rng(graph, positions, proposal['edits'], seed, sample)
        state = Atoms(numbers=numbers, positions=positions)
        dest = args.outdir/f'g{group:02d}'/f'{arm}_s{sample}_r{seed}'
        row = dict(group=group, arm=arm, strategy=strategy, sample=sample, seed=seed, random_seed=rng,
                   template_id=proposal['template_id'], origin=proposal.get('origin'))
        try:
            x, direction, meta = make_seed(state, mol, strategy, proposal, sample, rng, protocol.symbolic_seed_scale,
                                           protocol.encounter_policy, protocol.seed_features,
                                           protocol.seed_fit_max_nfev, protocol.seed_feature_ablation)
        except Exception as exc:
            if not is_recoverable_failure(exc):
                raise
            rows.append(dict(row, status='seed_generation_failed', evaluations=0, error=str(exc)))
            return
        seed_atoms = Atoms(numbers=numbers, positions=x)
        run_protocol = dataclasses.replace(protocol, random_seed=seed)    # as in a campaign run with this seed
        result = search_connection(seed_atoms, direction, CountedCalculator(backend, protocol.evaluations_per_attempt),
                                   dest, run_protocol, start['charge'])
        rows.append(dict(row, status=result['status'], evaluations=result['evaluations'],
                         ts_energy_eV=(result.get('ts') or {}).get('energy_eV'),
                         ts_positions_A=result.get('ts_positions_A'),
                         endpoint_graphs=sorted(e.get('graph_smiles') for e in result['endpoints']),
                         arrow_terms=meta.get('arrow_terms')))
        print(json.dumps({k: v for k, v in rows[-1].items() if k not in ('ts_positions_A', 'arrow_terms')}),
              flush=True)
    extra_seeds = [args.seeds[-1] + 1000*k for k in range(1, 64)]          # bond_edits budget match
    for g, proposals in enumerate(groups):
        for a, proposal in enumerate(proposals):
            for sample in args.samples:
                for seed in args.seeds:
                    attempt(g, f'arrows{a}', proposal, sample, seed, 'arrows')
        control_seeds = (args.seeds + extra_seeds)[:len(args.seeds)*len(proposals)]
        for sample in args.samples:
            for seed in control_seeds:
                attempt(g, 'bond_edits', proposals[0], sample, seed, 'bond_edits')
    report = dict(start=start['id'], graph=graph, atomic_numbers=numbers, charge=start['charge'],
                  seeds=args.seeds, samples=args.samples, protocol_json=str(args.protocol_json),
                  protocol_sha256=hashlib.sha256(args.protocol_json.read_bytes()).hexdigest(),
                  potential=potential, provenance=runtime_provenance(ROOT), initialization_evaluations=init_evaluations,
                  root_positions_A=positions.tolist(),
                  groups=[dict(group=g, edits=p[0]['edits'], predicted_graph=p[0]['predicted_graph'],
                               arrow_sets=[dict(template_id=q['template_id'], origin=q.get('origin'),
                                                arrows=q['arrows']) for q in p]) for g, p in enumerate(groups)],
                  rows=rows)
    (args.outdir/'same_edits.json').write_text(json.dumps(report, indent=2, default=list), encoding='utf-8')


if __name__ == '__main__':
    main()
