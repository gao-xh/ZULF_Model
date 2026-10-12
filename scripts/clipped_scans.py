"""Scans of a run that stay inside the ADC range over a time window: a keep list for average_scans.py --keep.

    python scripts/clipped_scans.py RUN_FOLDER OUT.json [--window 0.1,4.1] [--rail 28500]

A magnetometer output with a large, drifting offset can reach the ADC rail; the scan is then flat (no signal)
from that sample on, and averaging it in dilutes the signal and adds a step (diethylamine 2025-10-10: about half
of the 8001 scans clipped within 0.1-4.1 s). A scan is kept when no sample in the window has |value| >= rail.
Writes {"keep": [scan numbers], "first_clip_s": {scan: time or null}, settings}. RUN_FOLDER is read only.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from average_scans import decode_dat, read_ini  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--window", default="0.1,4.1", help="time window lo,hi (s) that must be free of clipping")
    ap.add_argument("--rail", type=float, default=28500.0, help="|value| at or above which a sample is clipped")
    args = ap.parse_args()
    run = Path(args.run)
    files = sorted((int(p.stem), p) for p in run.glob("*.dat") if p.stem.isdecimal())
    fs = float(read_ini(run / f"{files[0][0]}.ini")["sampling_rate_hz"])
    lo, hi = (float(v) for v in args.window.split(","))
    i0, i1 = int(round(lo * fs)), int(round(hi * fs))
    keep, first = [], {}
    for scan, path in files:
        x = decode_dat(path.read_bytes()).astype(np.int64)
        clipped = np.flatnonzero(np.abs(x[i0:]) >= args.rail)
        first[scan] = float((clipped[0] + i0) / fs) if len(clipped) else None
        if not len(clipped) or clipped[0] + i0 >= i1:
            keep.append(scan)
    out = {"run": str(run), "window_s": [lo, hi], "rail": args.rail, "scans_found": len(files),
           "scans_kept": len(keep), "keep": keep, "first_clip_s": first}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"{len(keep)} of {len(files)} scans free of clipping in {lo}-{hi} s -> {args.out}")


if __name__ == "__main__":
    main()
