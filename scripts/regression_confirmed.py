"""Regression over the confirmed samples: blind search (skeleton found?) and known-structure fit (chi2, time, 1J,
phase). Run after every change to zulf_hypothesis or the solver and compare with the previous run.

    ZULF_DATA_DIR=/path/to/fids python scripts/regression_confirmed.py --mode known --out runs/regression/NAME
    options: --mode known|blind|both, --samples id1,id2, --workers 4

Writes OUT/summary.json and OUT/summary.md, plus the automatic report of every fit (OUT/<id>_known*, _blind*).
"""
import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zulf_core.render.acquisition import Acquisition                        # noqa: E402
from zulf_core.solver import ObservedSpectrum                               # noqa: E402
from zulf_hypothesis import (blind_settings, fit_settings, fit_structure, propose_hypotheses,  # noqa: E402
                             search_hypotheses, write_report)
from zulf_hypothesis.fit import EXCHANGE_ELEMENTS                            # noqa: E402
from zulf_hypothesis.motifs import MOTIFS, _chain                           # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "configs" / "confirmed_samples.json"
# instrument prior for the fitted delay (configs: processing.phase_delay_bounds_s; env ZULF_DELAY_BOUNDS "lo,hi")
_bounds = os.environ.get("ZULF_DELAY_BOUNDS") or json.load(open(CONFIG))["processing"].get("phase_delay_bounds_s")
DELAY = tuple(float(x) for x in (_bounds.split(",") if isinstance(_bounds, str) else _bounds)) if _bounds else None


def observed_for(sample, proc, data_dir):
    x = np.load(Path(data_dir) / sample["file"]).astype(float)
    acq = Acquisition(proc["sampling_rate_hz"], len(x), start_sample=proc["start_sample"],
                      stop_sample=proc["stop_sample"], sg_window=proc["sg_window"], sg_order=proc["sg_order"],
                      remove_mean=proc["remove_mean"], apodization_rate_per_s=proc["apodization_rate_per_s"])
    return ObservedSpectrum.from_fid(x, acq, [tuple(r) for r in proc["ranges"]], zero_fill=proc["zero_fill"])


def structure_for(sample):
    if "motif" in sample:
        return MOTIFS[sample["motif"]].fragment(sample["one_bond"])
    c = sample["chain"]
    return _chain(sample["compound"], [tuple(g) for g in c["groups"]], [tuple(b) for b in c["bonds"]],
                  tuple(c.get("sym", ())))(sample["one_bond"])


def skeleton(fragment):
    """Multiset of (element, H count) over labelled sites carrying non-exchangeable protons."""
    elements = {s.label: s.element for s in fragment.sites if s.label_isotopes()}
    counts = Counter()
    for p in fragment.protons:
        if p.site in elements and elements[p.site] not in EXCHANGE_ELEMENTS and not p.exchangeable:
            counts[p.site] += p.size
    return sorted(Counter((elements[s], n) for s, n in counts.items()).items())


def one_bond(couplings, threshold=90.0):
    return sorted(round(v, 2) for k, v in couplings.items() if abs(v) >= threshold)


def run_known(sample, obs, workers, out):
    t = time.time()
    fit = fit_structure(structure_for(sample), obs, fit_settings(workers=workers, phase_delay_bounds_s=DELAY), report=str(out / f"{sample['id']}_known"),
                        route="both")
    best = fit.best
    red = best.chi2 / max(best.n - best.k, 1)
    row = {"best": best.key, "reduced_chi2": round(red, 2), "k": best.k, "seconds": round(time.time() - t),
           "one_bond_hz": one_bond(best.couplings), "findings": [f.code for f in best.findings],
           "phase0_deg": round(float(np.degrees(fit.phasing.get("phase0_rad", np.nan))), 1),
           "delay_ms": round(1e3 * fit.phasing.get("delay_s", np.nan), 3)}
    if fit.phased is not None and fit.phased.best is not None:
        pb = fit.phased.best
        row["phased_best"] = pb.key
        row["phased_reduced_chi2"] = round(pb.chi2 / max(pb.n - pb.k, 1), 2)
        gains = [complex(*g) for g in pb.summary["gains"]]
        # the phased model is corrected like the data, so its gain phase should equal the complex fit's phase
        # (modulo 180 deg); report the difference
        diff = np.angle(max(gains, key=abs)) - fit.phasing.get("phase0_rad", 0.0)
        row["phased_residual_phase_deg"] = round(float(np.degrees((diff + np.pi / 2) % np.pi - np.pi / 2)), 1)
    return row


def run_blind(sample, obs, workers, out, proc, rounds=1):
    t = time.time()
    ps = propose_hypotheses(obs, tuple(proc["instrument_lines_hz"]))
    res = search_hypotheses(obs, ps, blind_settings(workers, sample["id"], phase_delay_bounds_s=DELAY, rounds=rounds))
    write_report(res, obs, str(out / f"{sample['id']}_blind"))
    truth = skeleton(structure_for(sample))
    ranked = [e for e in res.ranked() if e.model.fragment is not None]
    hit = next((i + 1 for i, e in enumerate(ranked) if skeleton(e.model.fragment) == truth), None)
    best = res.best
    return {"best": best.key, "best_skeleton": skeleton(best.model.fragment) if best.model.fragment else None,
            "truth_skeleton": truth, "skeleton_correct": skeleton(best.model.fragment) == truth
            if best.model.fragment else False, "first_correct_rank": hit, "seconds": round(time.time() - t),
            "c_hat": round(res.c_hat, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="known", choices=["known", "blind", "both"])
    ap.add_argument("--samples", default="")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--blind-rounds", type=int, default=1, help="extension rounds of the blind search")
    ap.add_argument("--out", default=f"runs/regression/{time.strftime('%Y%m%d_%H%M%S')}")
    args = ap.parse_args()
    cfg = json.load(open(CONFIG))
    data_dir = os.environ.get("ZULF_DATA_DIR", "")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    wanted = set(filter(None, args.samples.split(",")))
    summary = {}
    for sample in cfg["samples"]:
        if wanted and sample["id"] not in wanted:
            continue
        obs = observed_for(sample, cfg["processing"], data_dir)
        row = {"compound": sample["compound"]}
        if args.mode in ("known", "both"):
            row["known"] = run_known(sample, obs, args.workers, out)
            print(sample["id"], "known", row["known"], flush=True)
        if args.mode in ("blind", "both"):
            row["blind"] = run_blind(sample, obs, args.workers, out, cfg["processing"], args.blind_rounds)
            print(sample["id"], "blind", row["blind"], flush=True)
        summary[sample["id"]] = row
        json.dump(summary, open(out / "summary.json", "w"), indent=1, default=str)
    lines = ["| sample | compound | known best | red. chi2 | time s | 1J | blind best | skeleton ok | time s |",
             "|---|---|---|---|---|---|---|---|---|"]
    for sid, r in summary.items():
        k, b = r.get("known", {}), r.get("blind", {})
        lines.append(f"| {sid} | {r['compound']} | {k.get('best', '')} | {k.get('reduced_chi2', '')} | "
                     f"{k.get('seconds', '')} | {k.get('one_bond_hz', '')} | {b.get('best', '')} | "
                     f"{b.get('skeleton_correct', '')} ({b.get('first_correct_rank', '')}) | {b.get('seconds', '')} |")
    (out / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
