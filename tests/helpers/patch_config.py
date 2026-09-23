#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path)
parser.add_argument("assignments", nargs="+")
args = parser.parse_args()
text = args.config.read_text()
for assignment in args.assignments:
    key, value = assignment.split("=", 1)
    pattern = re.compile(rf"(?mi)^(\s*{re.escape(key)}\s*=\s*)[^!\r\n]*(.*)$")
    text, count = pattern.subn(rf"\g<1>{value}\g<2>", text, count=1)
    if count != 1:
        raise SystemExit(f"configuration key not found exactly once: {key}")
args.config.write_text(text)
