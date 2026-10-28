"""Input schemas, checksums and readable CSV export for the compact replay."""

import csv
import gzip
import hashlib
import json
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
DATA_SHA256 = "d6e7b0461c4027e7000fd253e7f9ae8866a7fbba877c0bf46d5768248a8ecb60"
COORDINATES_SHA256 = "c7fd70dfdf6fe413dd7f8303b1e628ca8ef71494a276e5e285eab80479388221"


def load_payload(name, expected):
    """Load one compressed input only after verifying its recorded SHA256."""
    p = PACKAGE_ROOT / name
    if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
        raise ValueError("Input hash mismatch: " + name)
    return json.loads(gzip.decompress(p.read_bytes()))


def table_rows(d, name):
    """Expand a named column/row schema into ordinary record dictionaries."""
    t = d["tables"][name]
    return [dict(zip(t["columns"], r)) for r in t["rows"]]


def write_csv(path, rs, columns=None):
    """Write a readable CSV while retaining the recorded field order."""
    if not rs and (not columns):
        raise ValueError("Cannot infer empty CSV schema: " + str(path))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns or list(rs[0]))
        w.writeheader()
        w.writerows(rs)


def parse_boolean(x):
    """Accept only explicit Boolean values from the frozen tables."""
    if x is True or str(x).lower() == "true":
        return True
    if x is False or str(x).lower() == "false":
        return False
    raise ValueError("Invalid boolean " + repr(x))
