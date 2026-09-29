"""Recheck the prepared pool with the exact RDKit used by the QC environment."""
import json
from pathlib import Path
import sys
from rdkit import rdBase

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.event_graph import geometry_mol,graph_smiles
from mechbridge.chemistry import canonical

records=[json.loads(x) for x in (ROOT/'data/processed/feasibility_pool.jsonl').read_text().splitlines()]
failures=[]
for record in records:
    for side in ['reactant','product']:
        mol=geometry_mol(record['atomic_numbers'],record['positions_A'][side],record['charge'])
        if graph_smiles(mol)!=canonical(record[side+'_smiles']):
            failures.append({'event_id':record['event_id'],'side':side})
result={'rdkit_version':rdBase.rdkitVersion,'records':len(records),'failures':failures,'passed':not failures}
(ROOT/'reports/cohort_graph_audit_qc.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result),flush=True)
if failures: raise SystemExit(1)
