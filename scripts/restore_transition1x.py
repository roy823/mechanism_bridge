#!/usr/bin/env python3
"""Reassemble independently gzipped Transition1x parts; verify each part and whole SHA256."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

def restore(manifest_path, output):
    manifest_path=Path(manifest_path)
    manifest=json.loads(manifest_path.read_text())
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(f"Refusing to overwrite {output}")
    partial=output.with_name(output.name+'.assembling')
    whole=hashlib.sha256();total=0
    try:
        with partial.open('wb') as out:
            for expected_index,part in enumerate(manifest['parts'],1):
                if part['index']!=expected_index or Path(part['file']).name!=part['file']:
                    raise ValueError('Invalid part ordering or path')
                h=hashlib.sha256();n=0
                with gzip.open(manifest_path.parent/part['file'],'rb') as f:
                    while block:=f.read(4*1024*1024):
                        out.write(block);h.update(block);whole.update(block);n+=len(block)
                if n!=part['raw_bytes'] or h.hexdigest()!=part['raw_sha256']:
                    raise ValueError(f"Part failed verification: {part['file']}")
                total+=n
        if total!=manifest['raw_bytes'] or whole.hexdigest()!=manifest['raw_sha256']:
            raise ValueError('Reassembled file failed whole-file verification')
        partial.replace(output)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    print(json.dumps({'output':str(output),'bytes':total,'sha256':whole.hexdigest(),'verified':True}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path);p.add_argument('--output',type=Path,default=Path('Transition1x.h5'))
    a=p.parse_args();restore(a.manifest,a.output)
