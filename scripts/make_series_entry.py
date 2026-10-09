"""Average FID -> processed complex spectrum + one-entry series file for fit_joint_series / j_tuner.

    python scripts/make_series_entry.py --fid DATA/average_fid.npy --id ethylenediamine --out runs/series/eda
        [--crop 0.1] [--record 8] [--apodization 0.3] [--zero-fill 3] [--ranges 85,119.6;120.4,179.6 | bands]
        [--exclude 81.5,86] [--sampling-rate HZ] [--grid 20,380]

Recipe of the isopropylamine / amine analyses (docs/analysis/2026-10-03_isopropylamine_complex-fit.md): crop from
`--crop` s for `--record` s, exponential window, zero fill, phase from the instrument calibration of
configs/confirmed_samples.json (phase0 at the switching edge + delay offset). Without --ranges the fit range is the
whole grid with the mains harmonics (n x 60.06 Hz) and the instrument lines of the config +-0.4 Hz and the --exclude
intervals left out, so model lines where the data are empty are rejected by the fit (ranges taken from the data
alone let an ethanol fit put strong lines in an unfitted gap, 2026-10-08). `--ranges bands` keeps only the bands
where the smoothed |spectrum| exceeds 5 noise levels (noise from 330-380 Hz), widened by 3 Hz; check the ranges on
OUT/phased.png.
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
from zulf_processing import find_sampling_rate  # noqa: E402
from zulf_processing.series_spectrum import series_spectrum  # noqa: E402


def whole_ranges(lo, hi, lines, exclude=()):
    """The whole grid lo-hi as fit ranges, minus +-0.4 Hz around the power-line harmonics and instrument lines and
    minus the `exclude` intervals: a model line where the data are empty is then seen and rejected by the fit."""
    cuts = sorted([(v - 0.4, v + 0.4) for v in lines if lo < v < hi] + [tuple(e) for e in exclude])
    out, a = [], lo
    for c0, c1 in cuts:
        if c0 > a:
            out.append([round(float(a), 2), round(float(c0), 2)])
        a = max(a, c1)
    if a < hi:
        out.append([round(float(a), 2), round(float(hi), 2)])
    return out


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
    ap.add_argument("--sampling-rate", type=float, default=None,
                    help="Hz; default: scans.json or .ini next to the FID, else 4000")
    ap.add_argument("--grid", default="20,380", help="lo,hi of the processed frequency grid (Hz)")
    ap.add_argument("--crop", type=float, default=0.1, help="record start (s)")
    ap.add_argument("--record", type=float, default=8.0, help="record length (s)")
    ap.add_argument("--apodization", type=float, default=0.3, help="exponential window (1/s)")
    ap.add_argument("--zero-fill", type=int, default=3)
    ap.add_argument("--ranges", default="", help="lo,hi;lo,hi;... fit ranges (Hz); default: the whole grid minus "
                    "the power-line and instrument lines; 'bands': only the bands above 5 noise levels")
    ap.add_argument("--exclude", default="", help="lo,hi;... left out of the default ranges (e.g. 81.5,86)")
    ap.add_argument("--sg-window", type=float, default=0.0, help="drift filter length (s); default the plan's")
    ap.add_argument("--phase0-deg", type=float, default=None, help="zero-order phase; default the calibration")
    ap.add_argument("--delay-ms", type=float, default=None, help="first-order delay; default minus the edge")
    args = ap.parse_args()
    cfg = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]
    fs = args.sampling_rate or find_sampling_rate(args.fid, 4000.0)[0]
    grid = [float(v) for v in args.grid.split(",")]
    y = np.load(args.fid).astype(float)
    sp = series_spectrum(y, fs, args.crop, args.record, args.apodization, args.zero_fill, grid,
                         sg_window_s=args.sg_window or None, phase0_deg=args.phase0_deg, delay_ms=args.delay_ms)
    f, r, edge = sp["f"], sp["spectrum"], sp["edge_s"]
    lines = [60.06 * k for k in range(1, 7)] + [float(v) for v in cfg.get("instrument_lines_hz", []) if v % 60]
    if args.ranges == "bands":
        ranges = auto_ranges(f, r, lines, lo=grid[0], hi=grid[1])
    elif args.ranges:
        ranges = [[float(v) for v in q.split(",")] for q in args.ranges.split(";")]
    else:
        exclude = [[float(v) for v in q.split(",")] for q in args.exclude.split(";")] if args.exclude else []
        ranges = whole_ranges(grid[0], grid[1], lines, exclude)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "frequency.npy", f)
    np.save(out / "amplitude.npy", r)
    entry = {"id": args.id, "x": 1.0, "freq": str(out / "frequency.npy"), "values": str(out / "amplitude.npy"),
             "record": sp["acquisition"], "phasing": sp["phasing"],
             "ranges": ranges, "source_fid": str(Path(args.fid).expanduser().resolve())}
    json.dump([entry], open(out / "series.json", "w"), indent=1)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(16, 8))
        mid = 0.5 * (grid[0] + grid[1])
        for a, (lo, hi) in zip(ax, ((grid[0], mid), (mid, grid[1]))):
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
