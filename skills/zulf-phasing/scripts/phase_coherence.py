"""Model-free phasing of a ZULF spectrum from sign-agnostic peak coherence.

For each trial first-order delay the phases of the strongest peak points
(local maxima of |S| above a noise multiple inside the given bands) are
corrected with the project convention; the doubled angle removes the sign
ambiguity of ZULF lines. The coherence R(delay) = |sum w exp(2 i theta)| / sum w
(w = |S|^2) is maximised; phase0 follows from its angle, with the overall sign
chosen so the strongest peak is positive.

Usage:
    python phase_coherence.py FID.npy [--fs 4000] [--start 200] [--stop 4200] [--sg 201]
        [--bands 120-155,250-265] [--threshold 5] [--delay-range-ms 20] [--plot out.png]
Prints JSON with phase0_rad, delay_s, the coherence and the delay interval over
which the coherence stays above 90 % of its maximum (the delay uncertainty).
"""
import argparse
import json

import numpy as np
from scipy.signal import find_peaks

from zulf_core.render.acquisition import Acquisition, evaluate_spectrum, process_record
from zulf_core.render.phasing import phase_correct, reference_delay_s


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fid")
    p.add_argument("--fs", type=float, default=4000.0)
    p.add_argument("--start", type=int, default=200)
    p.add_argument("--stop", type=int, default=None)
    p.add_argument("--sg", type=int, default=201)
    p.add_argument("--bands", default="120-155,250-265")
    p.add_argument("--noise-band", default="180-250")
    p.add_argument("--threshold", type=float, default=5.0)
    p.add_argument("--delay-range-ms", type=float, default=20.0)
    p.add_argument("--plot", default=None)
    a = p.parse_args()
    x = np.load(a.fid).astype(float)
    acq = Acquisition(a.fs, len(x), start_sample=a.start, stop_sample=a.stop, sg_window=a.sg, sg_order=2,
                      remove_mean=True)
    df = a.fs / acq.n / 4
    f = np.arange(1, int(min(a.fs / 2, 600.0) / df)) * df
    s = evaluate_spectrum(process_record(x, acq), acq, f)
    lo, hi = (float(v) for v in a.noise_band.split("-"))
    noise = np.median(np.abs(s[(f > lo) & (f < hi)]))
    mask = np.zeros(len(f), bool)
    for item in a.bands.split(","):
        b0, b1 = (float(v) for v in item.split("-"))
        mask |= (f > b0) & (f < b1)
    peaks, _ = find_peaks(np.where(mask, np.abs(s), 0.0), height=a.threshold * noise, distance=max(1, int(0.8 / df)))
    if len(peaks) < 2:
        raise SystemExit("Fewer than two peaks above the threshold; lower --threshold or widen --bands.")
    fk, sk = f[peaks], s[peaks]
    w = np.abs(sk) ** 2
    taus = np.linspace(-a.delay_range_ms, a.delay_range_ms, 8001) * 1e-3
    theta = np.angle(sk)[None, :] - 2 * np.pi * fk[None, :] * (taus[:, None] + reference_delay_s(acq))
    coherence = np.abs((w[None, :] * np.exp(2j * theta)).sum(axis=1)) / w.sum()
    i = int(np.argmax(coherence))
    delay = float(taus[i])
    phase0 = float(np.angle((w * np.exp(2j * theta[i])).sum()) / 2)
    if phase_correct(sk, fk, phase0, delay, acq).real[np.argmax(w)] < 0:
        phase0 += np.pi
    phase0 = float((phase0 + np.pi) % (2 * np.pi) - np.pi)
    plateau = taus[coherence >= 0.9 * coherence[i]]
    out = {"phase0_rad": phase0, "delay_s": delay, "coherence": float(coherence[i]),
           "delay_interval_90pct_s": [float(plateau.min()), float(plateau.max())],
           "peaks_hz": [round(float(v), 3) for v in fk]}
    print(json.dumps(out, indent=1))
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(11, 6))
        axes[0].plot(taus * 1e3, coherence)
        axes[0].axvline(delay * 1e3, ls="--", color="k")
        axes[0].set_xlabel("delay (ms)")
        axes[0].set_ylabel("coherence")
        axes[1].plot(f, phase_correct(s, f, phase0, delay, acq).real / noise, lw=0.7)
        axes[1].axhline(0, color="0.7", lw=0.6)
        axes[1].set_xlabel("frequency (Hz)")
        axes[1].set_ylabel("real / noise")
        plt.tight_layout()
        plt.savefig(a.plot, dpi=120)


if __name__ == "__main__":
    main()
