"""One command for a new sample: standard processing, overview figure, hypotheses, search, report.

    python scripts/analyze_sample.py FID.npy --id SAMPLE [--out runs/blind/SAMPLE] [--workers 4]
        [--structure '{"motif": "ethyl", "one_bond": {"C1": 131, "C2": 125}}' | '{"chain": {...}, "one_bond": {...}}']

Without --structure: blind search (propose_hypotheses + search_hypotheses) with the settings used for the
confirmed samples; with --structure: fit_structure (every exchange regime and variant, fit-phased route too).
Processing and ranges come from configs/confirmed_samples.json ("processing"). Writes OUT/overview.png,
OUT/blind.{json,md,png} or OUT/structure.{json,md,png} (+ _phased), and prints the ranked table.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from zulf_core.render.acquisition import Acquisition, evaluate_spectrum, process_record   # noqa: E402
from zulf_core.solver import ObservedSpectrum                                              # noqa: E402
from zulf_hypothesis import (blind_settings, fit_settings, fit_structure, propose_hypotheses,  # noqa: E402
                             search_hypotheses, write_report)


def overview(x, acq, path, instrument):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    f = np.arange(5.0, 520.0, acq.sampling_rate_hz / acq.n / 4)
    mag = np.abs(evaluate_spectrum(process_record(x, acq), acq, f))
    mag /= np.median(mag[(f > 300) & (f < 500)])
    fig, ax = plt.subplots(2, 1, figsize=(13, 7))
    ax[0].plot(f, mag, lw=0.7, color="#00798c")
    ax[0].set_yscale("log")
    ax[1].plot(f, mag, lw=0.8, color="#00798c")
    for a in ax:
        for line in instrument:
            a.axvline(line, color="#cccccc", ls=":", lw=0.8, zorder=0)
        a.set_xlim(5, 520)
        a.set_xlabel("frequency (Hz)")
    ax[0].set_title("magnitude / noise (300-500 Hz), log")
    ax[1].set_title("linear (grey: instrument lines)")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fid")
    ap.add_argument("--id", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--structure", default="")
    args = ap.parse_args()
    proc = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]
    out = Path(args.out or f"runs/blind/{args.id}")
    out.mkdir(parents=True, exist_ok=True)
    x = np.load(args.fid).astype(float)
    acq = Acquisition(proc["sampling_rate_hz"], len(x), start_sample=proc["start_sample"],
                      stop_sample=proc["stop_sample"], sg_window=proc["sg_window"], sg_order=proc["sg_order"],
                      remove_mean=proc["remove_mean"], apodization_rate_per_s=proc["apodization_rate_per_s"])
    obs = ObservedSpectrum.from_fid(x, acq, [tuple(r) for r in proc["ranges"]], zero_fill=proc["zero_fill"])
    overview(x, acq, out / "overview.png", proc["instrument_lines_hz"])
    t = time.time()
    if args.structure:
        import regression_confirmed as reg
        spec = json.loads(args.structure)
        spec.setdefault("compound", args.id)
        fit = fit_structure(reg.structure_for(spec), obs, fit_settings(workers=args.workers),
                            report=str(out / "structure"), route="both")
        table, best = fit.table(), fit.best
    else:
        ps = propose_hypotheses(obs, tuple(proc["instrument_lines_hz"]))
        res = search_hypotheses(obs, ps, blind_settings(args.workers, args.id))
        write_report(res, obs, str(out / "blind"))
        table, best = res.table(), res.best
    print(f"{args.id}: {time.time() - t:.0f} s, best {best.key}")
    for r in table[:12]:
        print(f"{r['rank']:3d} {r['hypothesis'][:70]:70s} {r['variant']:6s} {r['status']:9s} d={r['delta']:8.1f} "
              f"k={r['k']:2d} {','.join(r['findings'])}{' DEMOTED' if r.get('demoted') else ''}")


if __name__ == "__main__":
    main()
