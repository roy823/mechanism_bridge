#!/usr/bin/env python3
"""Official Figshare article-archive fallback for expiring single-file redirects.

Archive may include weights; only requested dataset files are extracted.
"""
import argparse
import concurrent.futures
import json
from pathlib import Path
import shutil
import urllib.request
import zipfile
from download_data import ROOT, SOURCES, metadata, digest

def fetch(name, root):
    meta = metadata(name, root)
    aid, version = SOURCES[name]
    url = f"https://ndownloader.figshare.com/articles/{aid}/versions/{version}"
    folder = root / "raw" / name; folder.mkdir(parents=True, exist_ok=True)
    archive = folder / f"article_v{version}.zip"
    if archive.exists() and not zipfile.is_zipfile(archive):
        archive.unlink()  # Interrupted or malformed archive; never treat as complete data.
    if not archive.exists():
        part = archive.with_suffix(".zip.part")
        with urllib.request.urlopen(url, timeout=180) as response, part.open("wb") as out:
            shutil.copyfileobj(response, out, 4*1024*1024)
        part.replace(archive)
    if not zipfile.is_zipfile(archive):
        raise ValueError("Article archive incomplete; rerun to retry")
    results = []
    with zipfile.ZipFile(archive) as z:
        for f in meta["files"]:
            if "checkpoint" in f["name"]:
                continue
            target = folder / f["name"]
            if target.parent != folder:
                raise ValueError("Unsafe archive filename")
            # Only exact metadata-listed filenames, no extractall/path traversal.
            info = z.getinfo(f["name"])
            if info.file_size != f["size"]:
                raise ValueError("Archive entry differs from authoritative file size")
            temp = target.with_suffix(target.suffix + ".part")
            with z.open(info) as inp, temp.open("wb") as out:
                shutil.copyfileobj(inp, out, 4*1024*1024)
            if digest(temp, "md5") != f["computed_md5"]:
                temp.unlink(); raise ValueError("Extracted data failed MD5")
            temp.replace(target)
            receipt = {"dataset": name, "doi": meta["doi"], "license": meta["license"],
                       "url": url, "filename": f["name"], "file_id": f["id"],
                       "bytes": f["size"], "md5": f["computed_md5"], "sha256": digest(target),
                       "status": "downloaded_verified", "via": "official_article_archive"}
            target.with_name(target.name + ".receipt.json").write_text(json.dumps(receipt, indent=2))
            print(json.dumps(receipt), flush=True); results.append(receipt)
    archive.unlink()
    return results

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--datasets", nargs="+", choices=SOURCES, required=True)
    p.add_argument("--root", type=Path, default=ROOT / "data")
    a = p.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda n: fetch(n, a.root), a.datasets))
