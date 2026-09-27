"""Validate zulf_processing on the confirmed samples: per-dataset plan, switching edge, and phase criteria against
the phase of each sample's known-structure complex fit (regression summary).

    ZULF_DATA_DIR=... python scripts/validate_processing.py --reference runs/regression/final/summary.json \\
        --out runs/processing_validation

Reports, per sample and phase method, the phase error (deg, modulo 180) at 125, 200 and 250 Hz and the delay
minus the edge. The reference is itself a fit (model-dependent); a method is better when it agrees with it more
closely than the leave-one-out median calibration does. Writes OUT/summary.json, OUT/summary.md and a figure per
sample.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_processing import load_fid, process_dataset   # noqa: E402

FREQS = (125.0, 200.0, 250.0)


def error_deg(p, t, p0, t0):
    return [round(float(((np.degrees(p - p0) + 360 * f * (t - t0)) + 90) % 180 - 90), 1) for f in FREQS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", default="runs/regression/final/summary.json")
    ap.add_argument("--out", default="runs/processing_validation")
    ap.add_argument("--methods", default="entropy:edge_prior,entropy:edge_fixed,lines:edge_prior,lines:edge_fixed")
    args = ap.parse_args()
    cfg = json.load(open(ROOT / "configs" / "confirmed_samples.json"))
    ref = {k: (np.radians(v["known"]["phase0_deg"]), v["known"]["delay_ms"] * 1e-3)
           for k, v in json.load(open(args.reference)).items() if "known" in v}
    data = Path(os.environ.get("ZULF_DATA_DIR", ""))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    methods = [m.split(":") for m in args.methods.split(",")]
    summary = {}
    for s in cfg["samples"]:
        sid = s["id"]
        fid = load_fid(data / s["file"])
        rows = {}
        others = [k for k in ref if k != sid]
        if sid in ref and others:
            # robust leave-one-out calibration: median delay; references whose fitted delay is > 1 ms from it
            # (harmonic-band ambiguity, e.g. ethylenediamine +2.2 ms) are left out of the phase mean
            ct = float(np.median([ref[k][1] for k in others]))
            keep = [k for k in others if abs(ref[k][1] - ct) < 1e-3]
            ph = np.array([ref[k][0] + 2 * np.pi * FREQS[1] * (ref[k][1] - ct) for k in keep])
            c0 = float(np.angle(np.exp(2j * ph).mean()) / 2)
            rows["calibration (leave-one-out median)"] = {"errors_deg": error_deg(c0, ct, *ref[sid]),
                                                          "delay_ms": round(ct * 1e3, 3), "used": keep}
        for crit, mode in methods:
            d = process_dataset(fid, cfg["processing"]["sampling_rate_hz"], sid, phase_criterion=crit, delay_mode=mode)
            p = d.phase
            row = {"phase0_deg": round(float(np.degrees(p.phase0_rad)), 1), "delay_ms": round(p.delay_s * 1e3, 3),
                   "edge_ms": round(p.edge_delay_s * 1e3, 3) if p.edge_delay_s else None}
            if sid in ref:
                row["errors_deg"] = error_deg(p.phase0_rad, p.delay_s, *ref[sid])
            rows[f"{crit}:{mode}"] = row
            if crit == methods[0][0] and mode == methods[0][1]:
                d.figure(str(out / f"{sid}.png"))
                d.save_record(str(out / f"{sid}.json"))
        rows["plan"] = {"start_sample": d.plan.start_sample, "stop_sample": d.plan.stop_sample,
                        "reason": d.plan.reasons["start_sample"]}
        summary[sid] = rows
        print(sid, s["compound"], json.dumps(rows), flush=True)
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    lines = ["| sample | method | errors at 125/200/250 Hz (deg) | delay (ms) | edge (ms) |", "|---|---|---|---|---|"]
    for sid, rows in summary.items():
        for m, r in rows.items():
            if m == "plan":
                continue
            lines.append(f"| {sid} | {m} | {r.get('errors_deg', '')} | {r.get('delay_ms', '')} | {r.get('edge_ms', '')} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
