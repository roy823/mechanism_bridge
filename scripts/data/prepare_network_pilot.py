"""Join only starting reactants in existing RGD1 pool to replayable SynEPD arrows."""
import json
from pathlib import Path
import sys
from rdkit import RDLogger
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from mechbridge.symbolic_library import ArrowLibrary
from mechbridge.event_graph import geometry_mol


def main():
    RDLogger.DisableLog('rdApp.*')
    library = ArrowLibrary(ROOT/'data/raw/synepd/polar.json')
    pool = [json.loads(l) for l in (ROOT/'data/processed/feasibility_pool.jsonl').read_text().splitlines()]
    starts, matches, seen = [], [], set()
    for record in sorted(pool, key=lambda r: (len(r['atomic_numbers']), r['event_id'])):
        mol = geometry_mol(record['atomic_numbers'], record['positions_A']['reactant'], 0)
        proposals = library.propose(mol)
        if not proposals:
            continue
        matches.append(dict(rgd1_id=record['event_id'], reactant=record['reactant_smiles'],
                            templates=sorted({p['template_id'] for p in proposals})))
        if record['reactant_smiles'] in seen:
            continue
        seen.add(record['reactant_smiles'])
        starts.append(dict(id=record['event_id'], atomic_numbers=record['atomic_numbers'],
            positions_A=record['positions_A']['reactant'], charge=0, multiplicity=1,
            provenance=dict(dataset='RGD1', source_group=record['source_group'],
                            sha256=record['source_sha256'])))
    dest = ROOT/'data/processed/network_starts.jsonl'
    dest.write_text(''.join(json.dumps(s)+'\n' for s in starts), encoding='utf-8')
    report = dict(symbolic_audit=dict(library.audit), pool_size=len(pool),
                  matched_starting_records=matches, selected_starts=[s['id'] for s in starts],
                  matching='reactant_only; no requirement to match symbolic product or TS',
                  selection='all unique matching reactants in pre-existing <=8-atom pilot pool',
                  held_out_fields=['RGD1 product geometry', 'RGD1 TS geometry', 'RGD1 barriers'],
                  synepd_release=json.loads((ROOT/'data/raw/synepd/release-manifest.json').read_text())['dataset_release'])
    (ROOT/'reports/network_data_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
