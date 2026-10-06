"""Average FID -> processed complex spectrum + one-entry series file for fit_joint_series / j_tuner.

    python scripts/make_series_entry.py --fid DATA/average_fid.npy --id ethylenediamine --out runs/series/eda
        [--crop 0.1] [--record 8] [--apodization 0.3] [--zero-fill 3] [--ranges 85,119.6;120.4,179.6]

Recipe of the isopropylamine / amine analyses (docs/analysis/2026-10-03_isopropylamine_complex-fit.md): crop from
`--crop` s for `--record` s, exponential window, zero fill, phase from the instrument calibration of
configs/confirmed_samples.json (phase0 at the switching edge + delay offset). Without --ranges the fit ranges are the
bands where the smoothed |spectrum| exceeds 5 noise levels (noise from 330-380 Hz), widened by 3 Hz, with the mains
harmonics (n x 60.06 Hz) and the instrument lines of the config +-0.4 Hz left out; check them on OUT/phased.png.
Writes OUT/frequency.npy, OUT/amplitude.npy (complex), OUT/series.json (with "source_fid", the FID it came from, for
figures) and OUT/phased.png.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_core.render.phasing import correction_phasor, reference_delay_s  # noqa: E402
from zulf_processing import plan_for_dataset, process_dataset              # noqa: E402
from zulf_processing.diagnostics import switching_edge                     # noqa: E402


def auto_ranges(f, r, lines, threshold=5.0, noise_band=(330.0, 380.0), lo=20.0, hi=380.0):
    from scipy.ndimage import uniform_filter1d
    mag = uniform_filter1d(np.abs(r), 9)
    nb = (f > noise_band[0]) & (f < noise_band[1])
    noise = 1.4826 * np.median(np.abs(r.real[nb] - np.median(r.real[nb])))
    on = (f > lo) & (f < hi) & (mag > threshold * noise)
    for line in lines:
        on &= ~((f > line - 0.6) & (f < line + 0.6))
    idx = np.flatnonzero(on)
    bands = []
    if len(idx):
        start = prev = f[idx[0]]
        for v in f[idx[1:]]:
            if v - prev > 6.0:
                bands.append([start - 3, prev + 3])
                start = v
            prev = v
        bands.append([start - 3, prev + 3])
    out = []
    for b0, b1 in (b for b in bands if b[1] - b[0] > 7.0):
        a = b0
        for line in sorted(x for x in lines if b0 < x < b1):
            out.append([round(float(a), 1), round(line - 0.4, 1)])
            a = line + 0.4
        out.append([round(float(a), 1), round(float(b1), 1)])
    return [q for q in out if q[1] - q[0] > 1.0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fid", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sampling-rate", type=float, default=4000.0)
    ap.add_argument("--crop", type=float, default=0.1, help="record start (s)")
    ap.add_argument("--record", type=float, default=8.0, help="record length (s)")
    ap.add_argument("--apodization", type=float, default=0.3, help="exponential window (1/s)")
    ap.add_argument("--zero-fill", type=int, default=3)
    ap.add_argument("--ranges", default="", help="lo,hi;lo,hi;... fit ranges (Hz); default: automatic")
    args = ap.parse_args()
    cfg = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]
    cal = cfg["phase_calibration"]
    fs = args.sampling_rate
    y = np.load(args.fid).astype(float)
    edge = switching_edge(y, fs)["edge_time_s"] + cal["delay_offset_s"]
    start = int(round(args.crop * fs))
    plan = plan_for_dataset(len(y), fs, None, {"start_sample": start, "stop_sample": start + int(round(args.record * fs)),
                                               "zero_fill": args.zero_fill,
                                               "apodization_rate_per_s": args.apodization, "ranges": [[20.0, 380.0]]})
    ds = process_dataset(y, fs, plan=plan, phase_criterion=None)
    acq = plan.acquisition()
    f = ds.frequencies_hz
    phi0 = np.radians(cal["phase0_deg"])
    r = ds.spectrum * correction_phasor(f, phi0, -(edge + acq.time_origin_s) + reference_delay_s(acq))
    lines = [60.06 * k for k in range(1, 7)] + [float(v) for v in cfg.get("instrument_lines_hz", []) if v % 60]
    ranges = ([[float(v) for v in q.split(",")] for q in args.ranges.split(";")] if args.ranges
              else auto_ranges(f, r, lines))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "frequency.npy", f)
    np.save(out / "amplitude.npy", r)
    entry = {"id": args.id, "x": 1.0, "freq": str(out / "frequency.npy"), "values": str(out / "amplitude.npy"),
             "record": acq.to_dict(), "phasing": {"phase0_rad": float(phi0), "delay_s": float(-(edge + acq.time_origin_s))},
             "ranges": ranges, "source_fid": str(Path(args.fid).expanduser().resolve())}
    json.dump([entry], open(out / "series.json", "w"), indent=1)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(16, 8))
        for a, (lo, hi) in zip(ax, ((20, 200), (200, 380))):
            s = (f > lo) & (f < hi)
            a.plot(f[s], r.real[s], "k", lw=0.6, label="real")
            a.plot(f[s], r.imag[s], c="0.7", lw=0.5, label="imag")
            for q in ranges:
                a.axvspan(q[0], q[1], color="#cfe3f5", zorder=0)
            a.set_xlim(lo, hi)
        ax[0].legend()
        ax[0].set_title(f"{args.id}: calibrated phase at the edge {edge * 1e3:.3f} ms; blue: fit ranges")
        fig.tight_layout()
        fig.savefig(out / "phased.png", dpi=90)
    except ImportError:
        pass
    print(json.dumps({"series": str(out / "series.json"), "edge_ms": round(edge * 1e3, 3), "ranges": ranges}))


if __name__ == "__main__":
    main()
