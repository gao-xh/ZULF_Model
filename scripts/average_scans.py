"""Average the scans of one NMRduino run, with per-scan screening metrics and two disjoint half averages.

    python scripts/average_scans.py RUN_FOLDER OUT_DIR [--exclude-z 0] [--window 0.1,4.1]

RUN_FOLDER holds <n>.dat / <n>.ini as written by the instrument (read-only; nothing is written there).
Every scan is decoded with zulf_core.io.decode_dat. Per scan, on the window (s, after the record start):
  - deviation: RMS of (scan - mean of all scans) after removing each one's own mean and linear trend;
  - late_noise: RMS of first differences / sqrt(2) over the last 2 s of the record.
Robust z = (value - median) / (1.4826 MAD). With --exclude-z Z > 0, scans whose deviation z exceeds Z are left out
of the averages (and the mean used for the metric is recomputed once without them); default 0 keeps every scan.
Writes OUT_DIR/average_fid.npy (all kept scans), average_even.npy / average_odd.npy (kept scans by position,
disjoint halves for held-out checks), scans.json (metrics, flags, settings, source checksums) and scans.png.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_core.io import DECODER, decode_dat, read_ini   # noqa: E402


def detrended(y, t):
    a = np.polyfit(t, y, 1)
    return y - np.polyval(a, t)


def robust_z(v):
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med))
    return (v - med) / mad if mad > 0 else np.zeros_like(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--exclude-z", type=float, default=0.0)
    ap.add_argument("--window", default="0.1,4.1", help="metric window lo,hi (s)")
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.out)
    files = sorted((int(p.stem), p) for p in run.glob("*.dat") if p.stem.isdecimal())
    ids = np.array([i for i, _ in files])
    info = read_ini(run / f"{ids[0]}.ini")
    fs = float(info["sampling_rate_hz"])
    data = np.stack([decode_dat(p.read_bytes()) for _, p in files]).astype(np.float32)
    n = data.shape[1]
    t = np.arange(n) / fs
    lo, hi = (float(v) for v in args.window.split(","))
    w = (t >= lo) & (t < hi)
    late = t >= t[-1] - 2.0
    late_noise = np.std(np.diff(data[:, late].astype(float), axis=1), axis=1) / np.sqrt(2)

    def deviation(keep):
        mean_w = data[keep][:, w].astype(float).mean(axis=0)
        ref = detrended(mean_w, t[w])
        return np.array([np.sqrt(np.mean((detrended(row.astype(float), t[w]) - ref) ** 2)) for row in data[:, w]])

    keep = np.ones(len(ids), bool)
    dev = deviation(keep)
    z = robust_z(dev)
    if args.exclude_z > 0:
        keep = z <= args.exclude_z
        dev = deviation(keep)
        z = robust_z(dev)
        keep = z <= args.exclude_z
    kept = np.flatnonzero(keep)
    out.mkdir(parents=True, exist_ok=True)
    mean_all = data[kept].astype(float).mean(axis=0)
    np.save(out / "average_fid.npy", mean_all)
    np.save(out / "average_even.npy", data[kept[0::2]].astype(float).mean(axis=0))
    np.save(out / "average_odd.npy", data[kept[1::2]].astype(float).mean(axis=0))
    record = {
        "run": str(run), "decoder": DECODER, "sampling_rate_hz": fs, "points": int(n), "scans_found": int(len(ids)),
        "scans_kept": int(len(kept)), "exclude_z": args.exclude_z, "metric_window_s": [lo, hi],
        "first_ini_sha256": info["sha256"], "excluded_scans": [int(i) for i in ids[~keep]],
        "deviation_median": float(np.median(dev)), "late_noise_median": float(np.median(late_noise)),
        "scans": [{"scan": int(i), "deviation": float(d), "deviation_z": float(q), "late_noise": float(s), "kept": bool(k)}
                  for i, d, q, s, k in zip(ids, dev, z, late_noise, keep)],
        "average_sha256": hashlib.sha256((out / "average_fid.npy").read_bytes()).hexdigest(),
    }
    json.dump(record, open(out / "scans.json", "w"), indent=1)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
        ax[0].plot(ids, dev, ".", ms=2, c="k")
        ax[0].plot(ids[~keep], dev[~keep], "x", c="r", label="excluded")
        ax[0].set_ylabel("deviation (ADC)")
        ax[0].set_yscale("log")
        ax[1].plot(ids, late_noise, ".", ms=2, c="k")
        ax[1].set_ylabel("late noise (ADC)")
        ax[1].set_xlabel("scan")
        ax[0].set_title(f"{run.name}: {len(kept)} of {len(ids)} scans kept")
        fig.tight_layout()
        fig.savefig(out / "scans.png", dpi=90)
    except ImportError:
        pass
    print(json.dumps({k: record[k] for k in ("scans_found", "scans_kept", "deviation_median", "late_noise_median",
                                              "excluded_scans")}))


if __name__ == "__main__":
    main()
