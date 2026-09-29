import hashlib
import json
from pathlib import Path

def read_jsonl(path):
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)

def write_jsonl(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    count = 0
    with temp.open("w", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n")
            count += 1
    temp.replace(path)
    return count

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def event(source, source_id, kind, provenance):
    return {"schema_version": "0.1.0", "event_id": f"{source}:{source_id}",
            "source": source, "source_id": str(source_id), "kind": kind,
            "system": {"charge": None, "multiplicity": None, "environment": None},
            "symbolic": {}, "physical": {}, "provenance": provenance,
            "validation": {"arrows": "not_available", "stationarity": "not_checked",
                           "saddle_order": "not_checked", "connectivity": "not_checked",
                           "kinetics": "not_checked"}}
