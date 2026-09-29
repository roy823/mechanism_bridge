#!/usr/bin/env python3
"""Pinned Figshare downloads. No credentials; MD5+size verification; resumable bytes."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {"mech_uspto": (24797220, 2), "rgd1": (21066901, 6),
           "transition1x": (19614657, 4), "flower": (28359407, 3)}

def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def metadata(name, root):
    aid, version = SOURCES[name]
    path = root / "metadata" / f"{name}.json"
    if path.exists():
        data = json.loads(path.read_text())
    else:
        url = f"https://api.figshare.com/v2/articles/{aid}/versions/{version}"
        with urllib.request.urlopen(url, timeout=120) as r:
            data = json.load(r)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2))
    if data["id"] != aid or data["version"] != version:
        raise ValueError("Metadata does not match pinned dataset version")
    return data

def download(name, f, meta, root, retries=3):
    target = root / "raw" / name / f["name"]
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    result = {"dataset": name, "filename": f["name"], "file_id": f["id"],
              "doi": meta["doi"], "license": meta["license"],
              "url": f["download_url"], "expected_bytes": f["size"],
              "expected_md5": f["computed_md5"], "path": str(target.relative_to(root))}
    for attempt in range(retries):
        try:
            if not target.exists():
                offset = part.stat().st_size if part.exists() else 0
                if offset > f["size"]:
                    part.unlink(); offset = 0
                if offset < f["size"]:
                    # Figshare emits short-lived redirects. Request a fresh public link
                    # on each retry so intermediary caches cannot replay expired signatures.
                    fresh_url = f["download_url"] + "?download=1&request_time=" + str(time.time_ns())
                    req = urllib.request.Request(fresh_url, headers={"Range": f"bytes={offset}-"} if offset else {})
                    with urllib.request.urlopen(req, timeout=120) as response:
                        if offset and response.status == 206:
                            if not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                                raise ValueError("Invalid Content-Range")
                            mode = "ab"
                        elif response.status == 200:
                            mode = "wb"
                        else:
                            raise ValueError(f"Unexpected HTTP status {response.status}")
                        with part.open(mode) as out:
                            while chunk := response.read(4 * 1024 * 1024):
                                out.write(chunk)
                if part.stat().st_size != f["size"]:
                    raise ValueError("Size mismatch")
                if digest(part, "md5") != f["computed_md5"]:
                    part.unlink()
                    raise ValueError("MD5 mismatch; deleted corrupt partial")
                part.replace(target)
            if target.stat().st_size != f["size"] or digest(target, "md5") != f["computed_md5"]:
                raise ValueError("Existing file fails checksum; remove or rename it before retry")
            result.update(status="downloaded_verified", bytes=target.stat().st_size,
                          sha256=digest(target))
            break
        except Exception as exc:
            result.update(status="download_failed", error=f"{type(exc).__name__}: {exc}")
            print(json.dumps({"attempt": attempt + 1, **result}), flush=True)
            if attempt + 1 < retries:
                time.sleep(2)
    (target.parent / (target.name + ".receipt.json")).write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
    return result

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--datasets", nargs="+", choices=SOURCES, default=list(SOURCES))
    p.add_argument("--root", type=Path, default=ROOT / "data")
    p.add_argument("--max-file-mb", type=float, default=1500, help="Explicit download size cap; larger files recorded as skipped")
    p.add_argument("--checkpoints", action="store_true")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--retries", type=int, default=3)
    a = p.parse_args()
    tasks, results = [], []
    for name in a.datasets:
        m = metadata(name, a.root)
        for f in m["files"]:
            if "checkpoint" in f["name"] and not a.checkpoints:
                results.append({"dataset": name, "filename": f["name"], "status": "not_requested_model_weights"})
            elif f["size"] > a.max_file_mb * 1e6:
                results.append({"dataset": name, "filename": f["name"], "status": "skipped_size_limit", "bytes": f["size"]})
            else:
                tasks.append((name, f, m))
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures = [pool.submit(download, n, f, m, a.root, a.retries) for n, f, m in tasks]
        results.extend(f.result() for f in futures)
    (a.root / "download_manifest.json").write_text(json.dumps(results, indent=2))
    if any(r["status"] == "download_failed" for r in results):
        raise SystemExit(1)

if __name__ == "__main__":
    main()
