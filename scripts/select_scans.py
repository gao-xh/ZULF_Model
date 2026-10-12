"""Keep list of the scans of a run from their signal band, not from their broadband deviation.

    python scripts/select_scans.py RUN_FOLDER OUT.json [--band 110,270] [--window 0.1,4.1] [--z 5]
        [--rail 28500] [--pinned 0.9]

Per scan (RUN_FOLDER is read only):
  - pinned fraction: share of the samples in 0.1-1.1 s at |x| >= rail. A scan pinned at the ADC rail carries no
    signal; above `--pinned` it is left out;
  - band rms: rms of the scan band-passed to `--band` (Butterworth order 4, zero phase) over `--window`. Offsets,
    drift and the switching transient fall outside the band, so the metric sees the molecular signal and the
    noise there, not the magnetometer offset;
  - robust z of the band rms (median and 1.4826 MAD over the scans that are not pinned); a scan above `--z` is
    an outlier (glitches, interference bursts) and is left out.
Writes {"keep": [...], "scans": [per-scan metrics], settings}; `average_scans.py --keep OUT.json` averages them.
Diethylamine 2025-10-10 (10000 scans): the earlier broadband selection kept 3101 pinned scans and dropped 1434
good ones; this selection (6011 scans) gave higher line heights over noise (124.9 / 251.1 Hz: 15.8 / 9.6
against 14.1 / 6.7). Check the result on line heights before adopting it for another run.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_core.io import decode_dat, read_ini  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--band", default="110,270", help="signal band lo,hi (Hz)")
    ap.add_argument("--window", default="0.1,4.1", help="time window lo,hi (s) of the band rms")
    ap.add_argument("--z", type=float, default=5.0, help="robust z above which a scan is an outlier")
    ap.add_argument("--rail", type=float, default=28500.0)
    ap.add_argument("--pinned", type=float, default=0.9, help="pinned fraction above which a scan is left out")
    args = ap.parse_args()
    run = Path(args.run)
    files = sorted((int(p.stem), p) for p in run.glob("*.dat") if p.stem.isdecimal())
    fs = float(read_ini(run / f"{files[0][0]}.ini")["sampling_rate_hz"])
    lo_hz, hi_hz = (float(v) for v in args.band.split(","))
    sos = butter(4, [lo_hz, hi_hz], btype="band", fs=fs, output="sos")
    w0, w1 = (int(round(float(v) * fs)) for v in args.window.split(","))
    p0, p1 = int(round(0.1 * fs)), int(round(1.1 * fs))
    scans = []
    for scan, path in files:
        x = decode_dat(path.read_bytes()).astype(float)
        band = sosfiltfilt(sos, x)[w0:w1]
        scans.append({"scan": scan, "pinned_fraction": float(np.mean(np.abs(x[p0:p1]) >= args.rail)),
                      "band_rms": float(np.sqrt(np.mean(band ** 2)))})
    free = np.array([s["band_rms"] for s in scans if s["pinned_fraction"] <= args.pinned])
    med = float(np.median(free))
    mad = float(1.4826 * np.median(np.abs(free - med))) or 1.0
    keep = []
    for s in scans:
        s["band_z"] = (s["band_rms"] - med) / mad
        s["kept"] = s["pinned_fraction"] <= args.pinned and s["band_z"] <= args.z
        if s["kept"]:
            keep.append(s["scan"])
    out = {"run": str(run), "settings": vars(args), "scans_found": len(scans), "scans_kept": len(keep),
           "scans_pinned": int(sum(s["pinned_fraction"] > args.pinned for s in scans)),
           "band_rms_median": med, "band_rms_mad": mad, "keep": keep, "scans": scans}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"{len(keep)} of {len(scans)} scans kept ({out['scans_pinned']} pinned, "
          f"{len(scans) - len(keep) - out['scans_pinned']} band outliers) -> {args.out}")


if __name__ == "__main__":
    main()
