#!/usr/bin/env python3
"""Author's Zenodo deposit. HDF5 matches Figshare v6 byte-for-byte (MD5 checked)."""
import concurrent.futures
import json
import urllib.request
from pathlib import Path
from download_data import ROOT, download

root=ROOT/'data'
path=root/'raw/rgd1_zenodo/record.json'
path.parent.mkdir(parents=True,exist_ok=True)
with urllib.request.urlopen('https://zenodo.org/api/records/7860446',timeout=120) as r:meta=json.load(r)
path.write_text(json.dumps(meta,indent=2))
info={'doi':meta['doi'],'license':meta['metadata']['license']}
files=[]
for file in meta['files']:
    if file['key'] in ['RGD1_allrxns.h5','RGD1CHNO_AMsmiles.csv']:
        files.append({'name':file['key'],'id':file['id'],'size':file['size'],
                      'computed_md5':file['checksum'].split(':')[1], 'download_url':file['links']['self']})
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results=list(pool.map(lambda f: download('rgd1_zenodo',f,info,root),files))
if any(r['status']!='downloaded_verified' for r in results):raise SystemExit(1)
