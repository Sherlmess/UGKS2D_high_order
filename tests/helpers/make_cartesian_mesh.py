#!/usr/bin/env python3
import argparse
import struct
from pathlib import Path


def record(payload: bytes) -> bytes:
    marker = struct.pack("=i", len(payload))
    return marker + payload + marker


parser = argparse.ArgumentParser()
parser.add_argument("path", type=Path)
parser.add_argument("--nx", type=int, required=True)
parser.add_argument("--ny", type=int, required=True)
parser.add_argument("--surf", type=Path)
args = parser.parse_args()
nxp, nyp = args.nx + 1, args.ny + 1
xs = [i / args.nx for j in range(nyp) for i in range(nxp)]
ys = [j / args.ny for j in range(nyp) for i in range(nxp)]
payload1 = struct.pack("=ii", nxp, nyp)
payload2 = struct.pack("=" + "d" * (2 * nxp * nyp), *(xs + ys))
args.path.write_bytes(record(payload1) + record(payload2))
if args.surf is not None:
    surface_x = [(i + 0.5) / args.nx for i in range(args.nx)]
    surface_y = [0.0 for _ in range(args.nx)]
    surface_flow = [value for _ in range(args.nx) for value in (1.0, 0.1, 0.0, 1.0, 0.5)]
    surface_header = struct.pack("=ii", args.nx, 7)
    surface_payload = struct.pack("=" + "d" * (7 * args.nx), *(surface_x + surface_y + surface_flow))
    args.surf.write_bytes(record(surface_header) + record(surface_payload))
