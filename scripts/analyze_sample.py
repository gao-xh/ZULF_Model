"""One command for a new sample: standard processing, overview figure, hypotheses, search, report.

    python scripts/analyze_sample.py FID.npy --id SAMPLE [--out runs/blind/SAMPLE] [--workers 4]
        [--structure '{"motif": "ethyl", "one_bond": {"C1": 131, "C2": 125}}' | '{"chain": {...}, "one_bond": {...}}']

Without --structure: blind search (propose_hypotheses + search_hypotheses) with the settings used for the
confirmed samples; with --structure: fit_structure (every exchange regime and variant, fit-phased route too).
Processing: zulf_processing.process_dataset per dataset (defaults from configs/confirmed_samples.json
"processing"; crop after this dataset's ringing; phase per dataset). Writes OUT/processing.{png,json}, OUT/overview.png,
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
from zulf_processing import load_fid, process_dataset                                      # noqa: E402
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


DELAY = None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fid")
    ap.add_argument("--id", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--structure", default="")
    ap.add_argument("--labeling", default="natural", choices=["natural", "15N", "2H-exchange", "unknown"],
                    help="isotope labelling of the sample for --structure (D54); unknown: fit every applicable "
                         "labelling and rank them on one scale")
    ap.add_argument("--phase", default="", help="phase criterion for the processed spectrum: calibration (default "
                    "when the config has phase_calibration), entropy or lines")
    ap.add_argument("--delay-bounds", default="", help="instrument prior for the fitted delay in s, 'lo,hi'")
    args = ap.parse_args()
    proc = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]
    global DELAY
    bounds = args.delay_bounds or proc.get("phase_delay_bounds_s")
    DELAY = tuple(float(x) for x in (bounds.split(",") if isinstance(bounds, str) else bounds)) if bounds else None
    out = Path(args.out or f"runs/blind/{args.id}")
    out.mkdir(parents=True, exist_ok=True)
    x = load_fid(args.fid)
    # per-dataset processing (zulf_processing, D44): plan with reasons, switching edge, per-dataset phase
    criterion = args.phase or ("calibration" if proc.get("phase_calibration") else "entropy")
    data = process_dataset(x, proc["sampling_rate_hz"], args.id, defaults=proc, phase_criterion=criterion)
    data.figure(str(out / "processing.png"))
    data.save_record(str(out / "processing.json"))
    obs, acq = data.observed, data.plan.acquisition()
    overview(x, acq, out / "overview.png", data.plan.instrument_lines_hz)
    print(f"{args.id}: crop {data.plan.start_sample}-{data.plan.stop_sample} ({data.plan.reasons['start_sample']}), "
          f"edge {data.diagnostics.edge_time_s * 1e3 if data.diagnostics.edge_time_s else float('nan'):.3f} ms, "
          f"phase0 {np.degrees(data.phase.phase0_rad):.1f} deg, delay {data.phase.delay_s * 1e3:.3f} ms")
    t = time.time()
    if args.structure:
        import regression_confirmed as reg
        spec = json.loads(args.structure)
        spec.setdefault("compound", args.id)
        settings = fit_settings(workers=args.workers, phase_delay_bounds_s=DELAY)
        if args.labeling == "natural":
            fit = fit_structure(reg.structure_for(spec), obs, settings, report=str(out / "structure"), route="both")
        else:
            from zulf_hypothesis.fit import fit_structure_labelings
            cmp = fit_structure_labelings(reg.structure_for(spec), obs, args.labeling, settings=settings,
                                          report=str(out / "structure"), route="both")
            print(f"{args.id}: labelling hypotheses (common scale, conditional on this structure):")
            for row in cmp.rows:
                print(f"  {row['labeling']:12s} {row['model'][:50]:50s} {row['variant']:6s} d={row['delta']:10.1f} "
                      f"k={row['k']}")
            json.dump(cmp.rows, open(out / "labelings.json", "w"), indent=1)
            fit = cmp.fits[cmp.best]
        table, best = fit.table(), fit.best
    else:
        ps = propose_hypotheses(obs, tuple(proc["instrument_lines_hz"]))
        res = search_hypotheses(obs, ps, blind_settings(args.workers, args.id, phase_delay_bounds_s=DELAY))
        write_report(res, obs, str(out / "blind"))
        table, best = res.table(), res.best
    print(f"{args.id}: {time.time() - t:.0f} s, best {best.key}")
    for r in table[:12]:
        print(f"{r['rank']:3d} {r['hypothesis'][:70]:70s} {r['variant']:6s} {r['status']:9s} d={r['delta']:8.1f} "
              f"k={r['k']:2d} {','.join(r['findings'])}{' DEMOTED' if r.get('demoted') else ''}")


if __name__ == "__main__":
    main()
