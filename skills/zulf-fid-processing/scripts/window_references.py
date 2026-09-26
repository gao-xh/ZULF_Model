"""Residual scale for a processing window: no-lines residual and noise floor.

Usage:
    python window_references.py FID.npy [--fs 4000] [--start 200] [--stop 4200] [--sg 201]
        [--zero-fill 4] [--ranges 100-119,121-179,181-293,295-320] [--signal 127-135,142-152]
        [--phase0 RAD --delay S]   (phased real-only scale)

The no-lines residual is what a per-band linear background alone leaves; the
noise floor is the residual expected from noise only, estimated from points
outside the --signal windows. Report fits as (null - r) / (null - floor).
"""
import argparse
import json

import numpy as np

from zulf_core.render.acquisition import Acquisition
from zulf_core.render.phasing import phase_correct
from zulf_core.solver import ObservedSpectrum


def bands(text):
    return [tuple(float(v) for v in item.split("-")) for item in text.split(",") if item]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fid")
    p.add_argument("--fs", type=float, default=4000.0)
    p.add_argument("--start", type=int, default=200)
    p.add_argument("--stop", type=int, default=None)
    p.add_argument("--sg", type=int, default=201)
    p.add_argument("--zero-fill", type=int, default=4)
    p.add_argument("--ranges", default="100-119,121-179,181-293,295-320")
    p.add_argument("--signal", default="127-135,142-152,256-264")
    p.add_argument("--phase0", type=float, default=None)
    p.add_argument("--delay", type=float, default=0.0)
    a = p.parse_args()
    x = np.load(a.fid).astype(float)
    acq = Acquisition(a.fs, len(x), start_sample=a.start, stop_sample=a.stop, sg_window=a.sg, sg_order=2,
                      remove_mean=True)
    obs = ObservedSpectrum.from_fid(x, acq, bands(a.ranges), zero_fill=a.zero_fill)
    f, y, band = obs.frequencies_hz[obs.selected], obs.values[obs.selected], obs.band_index[obs.selected]
    if a.phase0 is not None:
        y = phase_correct(y, f, a.phase0, a.delay, acq).real
    resid = []
    for k in np.unique(band):
        m = band == k
        u = (f[m] - f[m].mean()) / (np.ptp(f[m]) / 2)
        design = np.column_stack([np.ones(m.sum()), u])
        resid.append(y[m] - design @ np.linalg.lstsq(design, y[m], rcond=None)[0])
    resid = np.concatenate(resid)
    signal = np.zeros(len(f), bool)
    for lo, hi in bands(a.signal):
        signal |= (f > lo) & (f < hi)
    floor = np.sqrt(np.mean(np.abs(resid[~signal]) ** 2) * len(y)) / np.linalg.norm(y)
    print(json.dumps({"null": float(np.linalg.norm(resid) / np.linalg.norm(y)), "noise_floor": float(floor),
                      "points": int(len(y)), "mode": "phased real" if a.phase0 is not None else "complex"}))


if __name__ == "__main__":
    main()
