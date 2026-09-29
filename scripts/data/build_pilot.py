#!/usr/bin/env python3
"""Build a deterministic real-data pilot from downloaded archives; no fabricated pairs."""
import argparse
from collections import Counter
import io
import json
from pathlib import Path
import zipfile
from mechbridge.adapters import flow_lines, provenance, mech_csv, rgd1, transition1x
from mechbridge.io import write_jsonl, read_jsonl
from mechbridge.pairing import pair_candidates, split_records

ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--limit',type=int,default=1000)
a=p.parse_args()
out=ROOT/'data/processed';out.mkdir(exist_ok=True)
summary={'pilot_limit_per_source':a.limit,'sampling':'first records, engineering smoke only, not a representative benchmark'}
symbolic=[];physical=[]
archive=ROOT/'data/raw/flower/data.zip'
if archive.exists():
    with zipfile.ZipFile(archive) as z:
        member='data/flower_dataset/train.txt'
        prov={**provenance(archive,'10.6084/m9.figshare.28359407.v3'),'zip_member':member,
              'source_split':'train','variant':'flower_dataset'}
        with z.open(member) as f:
            symbolic=list(flow_lines(io.TextIOWrapper(f,encoding='utf8'),'flower',prov,a.limit,'flower_dataset:'))
    summary['flower']=write_jsonl(out/'flower.jsonl',symbolic)
mech=ROOT/'data/raw/mech_uspto/mech-USPTO-31k.csv'
if mech.exists():
    records=list(mech_csv(mech,'updated_reaction','mechanistic_label',limit=a.limit))
    summary['mech_uspto']=write_jsonl(out/'mech_uspto.jsonl',records)
    summary['mech_arrow_parse_failures']=sum('arrow_parse_error' in r['symbolic'] for r in records)
    # Excluded from elementary-event matching until a supported grouping is available.
rgd=ROOT/'data/raw/rgd1_zenodo/RGD1_allrxns.h5'
if not rgd.exists():rgd=ROOT/'data/raw/rgd1/RGD1_CHNO.h5'
if rgd.exists():
    mapped_csv=ROOT/'data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv'
    physical=list(rgd1(rgd,a.limit,mapped_csv if mapped_csv.exists() else None))
    if 'zenodo' in str(rgd):
        for r in physical:r['provenance']['retrieval_doi']='10.5281/zenodo.7860446'
    summary['rgd1']=write_jsonl(out/'rgd1.jsonl',physical)
    summary['rgd1_mapped_csv_graph_agreement']=dict(Counter(str(r['provenance'].get('mapped_csv_endpoint_graphs_agree')) for r in physical))
t1x=ROOT/'data/raw/transition1x/Transition1x.h5'
if t1x.exists():
    records=list(transition1x(t1x,a.limit))
    summary['transition1x']=write_jsonl(out/'transition1x.jsonl',records)
summary['exact_graph_candidates']=write_jsonl(out/'pair_candidates.jsonl',pair_candidates(symbolic,physical))
summary['verified_paired_events']=0
summary['flower_audit']={key:sum(bool(r.get('graph_audit',{}).get(key)) for r in symbolic)
                          for key in ['composition_conserved','charge_conserved','complete_unique_atom_maps','mapped_atom_identity_conserved']}
write_jsonl(out/'pilot_split.jsonl',split_records(symbolic+physical))
(ROOT/'reports/pilot_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
