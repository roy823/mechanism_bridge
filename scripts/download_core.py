#!/usr/bin/env python3
"""Download the five pinned core files with aria2, then verify size, MD5 and SHA256.

Run again to resume. The original inventory supplies expected hashes, never local
completion status. This run writes reports/local_download_inventory.json.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports/local_download_inventory.json"
LOCK = threading.Lock()


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    return md5.hexdigest(), sha.hexdigest()


def check(path, spec):
    if path.stat().st_size != spec["bytes"]:
        raise ValueError(f"Size mismatch: {path.name}")
    md5, sha = hashes(path)
    if md5 != spec.get("md5", spec.get("expected_md5")) or sha != spec["sha256"]:
        raise ValueError(f"Checksum mismatch: {path.name}; retained for inspection")
    return md5, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aria2", type=Path, default=ROOT / "tools/aria2/aria2c.exe")
    parser.add_argument("--connections", type=int, default=8)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.connections <= 16 or args.workers < 1:
        parser.error("connections must be 1..16 and workers must be positive")
    if not args.verify_only and not args.aria2.is_file():
        parser.error("Install the official aria2 binary in tools/aria2 or pass --aria2")
    specs = json.loads((ROOT / "reports/data_inventory.json").read_text(encoding="utf-8"))["files"]
    states = [{"dataset": s["dataset"], "path": s["relative_path"],
               "expected_bytes": s["bytes"], "expected_sha256": s["sha256"],
               "doi": s["doi"], "license": s["license"], "status": "pending"} for s in specs]
    if not args.verify_only:
        (ROOT / "reports/download_process.json").write_text(json.dumps({
            "pid": os.getpid(), "script": "scripts/download_core.py",
            "started_at_utc": datetime.now(timezone.utc).isoformat()
        }, indent=2), encoding="utf-8")

    def update(index, **fields):
        with LOCK:
            states[index].update(fields)
            report = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
                      "expected_total_bytes": sum(s["bytes"] for s in specs),
                      "all_verified": all(s["status"] == "downloaded_verified" for s in states),
                      "files": states}
            temp = REPORT.with_suffix(".json.tmp")
            temp.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
            temp.replace(REPORT)
            print(json.dumps({"dataset": specs[index]["dataset"], **fields}), flush=True)

    def fetch(index):
        spec = specs[index]
        target = (ROOT / spec["relative_path"]).resolve()
        target.relative_to((ROOT / "data/raw").resolve())
        part = target.with_name(target.name + ".part")
        try:
            if not target.exists():
                if args.verify_only:
                    update(index, status="incomplete" if part.exists() else "missing")
                    return
                target.parent.mkdir(parents=True, exist_ok=True)
                if shutil.disk_usage(target.parent).free < spec["bytes"] + 1024**3:
                    raise OSError("Insufficient free disk space")
                if spec["dataset"] == "rgd1_zenodo":
                    url = spec["url"]
                else:
                    url = f"https://ndownloader.figshare.com/files/{spec['file_id']}"
                update(index, status="downloading", url=url)
                log = ROOT / "reports" / ("download_" + target.name + ".log")
                command = [str(args.aria2.resolve()), "--no-conf=true", "--continue=true",
                           "--disable-ipv6=true",
                           "--auto-file-renaming=false", "--allow-overwrite=false",
                           "--file-allocation=none", "--max-tries=10", "--retry-wait=5",
                           "--connect-timeout=30", "--timeout=60", "--min-split-size=8M",
                           "--lowest-speed-limit=128K",
                           f"--max-connection-per-server={args.connections}",
                           f"--split={args.connections}", "--summary-interval=30",
                           "--console-log-level=warn", "--download-result=full",
                           "--show-console-readout=false",
                           "--enable-color=false", "--dir=" + str(target.parent),
                           "--out=" + part.name]
                with log.open("a", encoding="utf-8") as output:
                    failures_without_progress = 0
                    while True:
                        # Figshare redirects expire in seconds. Refresh the public
                        # URL when a transfer fails; aria2 retains completed pieces.
                        fresh_url = url + "?download=1&request_time=" + str(time.time_ns())
                        before = part.stat().st_mtime_ns if part.exists() else None
                        completed = subprocess.run(command + [fresh_url], stdout=output,
                                                   stderr=subprocess.STDOUT)
                        if completed.returncode == 0:
                            break
                        after = part.stat().st_mtime_ns if part.exists() else None
                        failures_without_progress = (0 if after != before else failures_without_progress + 1)
                        if failures_without_progress >= 5:
                            raise RuntimeError(f"aria2 exited {completed.returncode}; see {log.name}")
                        time.sleep(5)
                update(index, status="verifying")
                md5, sha = check(part, spec)
                part.replace(target)
            else:
                update(index, status="verifying")
                md5, sha = check(target, spec)
            receipt = {**states[index], "status": "downloaded_verified", "bytes": target.stat().st_size,
                       "md5": md5, "sha256": sha, "verified_at_utc": datetime.now(timezone.utc).isoformat()}
            target.with_name(target.name + ".receipt.json").write_text(
                json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8")
            update(index, **{k: receipt[k] for k in ["status", "bytes", "md5", "sha256", "verified_at_utc"]})
        except Exception as exc:
            update(index, status="failed", error=f"{type(exc).__name__}: {exc}")

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(fetch, range(len(specs))))
    if not all(s["status"] == "downloaded_verified" for s in states):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
