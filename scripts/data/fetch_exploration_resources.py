"""Fetch pinned author data and AIMNet2-rxn weights with SHA256 receipts."""
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[2]


def get(url):
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    return r


def main():
    gh = '4fefc016fd4d4e4305dec92586e593bafe3b3e5b'
    hf = '13e9de20a1e878a18a371e94dfb2176629281cb6'
    files = []
    for name in ['data/polar.json', 'data/release-manifest.json', 'LICENSE']:
        files.append((f'https://raw.githubusercontent.com/TieuLongPhan/SynEPD/{gh}/{name}',
                      ROOT / 'data/raw/synepd' / Path(name).name, gh))
    for name in ['config.json', 'ensemble_0.safetensors', 'README.md']:
        files.append((f'https://huggingface.co/isayevlab/aimnet2-rxn/resolve/{hf}/{name}',
                      ROOT / 'models/aimnet2-rxn' / name, hf))
    receipt = []
    for url, path, revision in files:
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            data = get(url).content
            temp = path.with_suffix(path.suffix + '.part')
            temp.write_bytes(data)
            temp.replace(path)
        data = path.read_bytes()
        receipt.append(dict(path=path.relative_to(ROOT).as_posix(), url=url,
                            revision=revision, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
        print(json.dumps(receipt[-1]), flush=True)
    manifest = json.loads((ROOT/'data/raw/synepd/release-manifest.json').read_text())
    for spec in manifest['sources']:
        path = ROOT/'data/raw/synepd'/spec['name']
        if hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
            raise ValueError('SynEPD data differs from author release manifest')
    (ROOT / 'reports/exploration_resource_receipt.json').write_text(
        json.dumps(receipt, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
