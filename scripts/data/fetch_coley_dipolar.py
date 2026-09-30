"""Download the pinned author table; no reference geometries enter exploration."""
import hashlib
import argparse
import json
from pathlib import Path
from pathlib import PurePosixPath
import tarfile
import requests
ROOT=Path(__file__).resolve().parents[2]
REVISION='1608a64499779bf8f89a886308d97c6a830d084c'
FILES={
    'LICENSE':'4ef9682f72ede67442b9fd976812c34884ce9bbc35b4eb56cb1642f6bdca299d',
    'README.md':'ca5df496744c9d53775885d04b0d79e464a90c5217f554b5b1ad0c282bae69d6',
    'final_data_files/full_dataset.csv':'5f08c6b2dc04dd3ae11c74b2e0251a595b230915c7196d7a67e96bfac7f1e12d'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-profiles',action='store_true',
        help='Also fetch the 94 MB author geometry archive for separate oracle diagnostics')
    args=parser.parse_args()
    directory=ROOT/'data/raw/coley_dipolar';directory.mkdir(parents=True,exist_ok=True)
    receipts=[]
    files=dict(FILES)
    if args.reference_profiles:
        files['final_data_files/full_data_profiles.tar.gz']='c42f6cbe0031521d159844f845584c27435f778d4de6516faa117717dd6766ac'
    for source,digest in files.items():
        url=f'https://raw.githubusercontent.com/coleygroup/dipolar_cycloaddition_dataset/{REVISION}/{source}'
        path=directory/Path(source).name
        if path.exists():data=path.read_bytes()
        else:
            response=requests.get(url,timeout=120);response.raise_for_status();data=response.content
        if hashlib.sha256(data).hexdigest()!=digest:raise ValueError('Pinned checksum mismatch: '+source)
        if not path.exists():path.write_bytes(data)
        receipts.append(dict(url=url,file=path.relative_to(ROOT).as_posix(),bytes=len(data),sha256=digest))
    (directory/'receipt.json').write_text(json.dumps(dict(revision=REVISION,files=receipts),indent=2),encoding='utf-8')
    if args.reference_profiles:
        selected=[]
        with tarfile.open(directory/'full_data_profiles.tar.gz','r:gz') as archive:
            for member in archive:
                parts=PurePosixPath(member.name).parts
                if (len(parts)!=3 or parts[0]!='full_dataset_profiles' or parts[1] not in ('3216','3217')
                    or not member.isfile() or parts[2]=='.DS_Store'):
                    continue
                if '\\' in parts[2] or parts[2] in ('.','..'):
                    raise ValueError('Unsafe reference profile filename')
                target=directory/'selected_profiles'/parts[1]/parts[2]
                data=archive.extractfile(member).read()
                if target.exists() and target.read_bytes()!=data:raise ValueError('Existing reference profile differs')
                target.parent.mkdir(parents=True,exist_ok=True)
                if not target.exists():target.write_bytes(data)
                selected.append(dict(member=member.name,file=target.relative_to(ROOT).as_posix(),
                    sha256=hashlib.sha256(data).hexdigest()))
        if len(selected)!=12:raise ValueError('Unexpected pinned profile contents')
        (directory/'selected_profiles/receipt.json').write_text(json.dumps(selected,indent=2),encoding='utf-8')
    print(json.dumps(receipts,indent=2))


if __name__=='__main__':main()
