#!/usr/bin/env python3
"""Read-only package integrity verification; no numerical calls."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def verify():
    identity = json.loads((ROOT/'provenance/SOURCE_IDENTITY.json').read_text())
    for name, expected in identity['source_files'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
    checks = 0
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        expected, name = line.split('  ', 1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name
        checks += 1
    print(f'PASS: {checks} package hashes; {len(identity["source_files"])} frozen source/build files; zero numerical steps')

if __name__ == '__main__':
    verify()
