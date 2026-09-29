#!/usr/bin/env python3
"""Export small physical events to a reviewable QC queue. Never guess charge/spin."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from ase import Atoms
from ase.io import write
from mechbridge.io import read_jsonl,write_jsonl

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('events');p.add_argument('outdir',type=Path)
p.add_argument('--max-atoms',type=int,default=12);p.add_argument('--limit',type=int,default=50)
p.add_argument('--charge',type=int,required=True);p.add_argument('--multiplicity',type=int,required=True)
p.add_argument('--method',default='wb97x');p.add_argument('--basis',default='6-31g(d)')
a=p.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
jobs=[]
for event in read_jsonl(a.events):
    physical=event['physical'];numbers=physical.get('atomic_numbers',[])
    if not numbers or len(numbers)>a.max_atoms or 'ts' not in physical.get('positions_A',{}):continue
    known_charge=event['system'].get('charge')
    if known_charge is not None and known_charge!=a.charge:continue
    if (sum(numbers)-a.charge-(a.multiplicity-1))%2:continue
    folder=a.outdir/f'job_{len(jobs):05d}';folder.mkdir(exist_ok=True)
    write(folder/'ts.xyz',Atoms(numbers=numbers,positions=physical['positions_A']['ts']))
    config={'event_id':event['event_id'],'charge':a.charge,'multiplicity':a.multiplicity,
            'method':a.method,'basis':a.basis,'environment':'gas_phase',
            'electronic_state_assignment':'explicit_CLI_assertion_requires_scientific_review',
            'source_provenance':event['provenance'],'status':'pending_calculation',
            'command':['mechbridge','verify',str(folder/'ts.xyz'),str(folder/'results'),
                       '--dimer','--charge',str(a.charge),'--multiplicity',str(a.multiplicity),
                       '--method',a.method,'--basis',a.basis]}
    (folder/'job.json').write_text(json.dumps(config,indent=2));jobs.append(config)
    if len(jobs)>=a.limit:break
write_jsonl(a.outdir/'jobs.jsonl',jobs)
print(json.dumps({'jobs_prepared':len(jobs),'jobs_executed':0}))
