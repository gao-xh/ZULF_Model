"""Refine the two-isotopologue isopropylamine skeleton model on an averaged FID.

Usage:
    python scripts/refine_isopropylamine.py AVERAGE.npy 0.ini OUTPUT_DIR [--budget-s 1800]

The starting J values are the defaults of the earlier ZULF_Analysis_Tools
isopropylamine model; vicinal H-H couplings are tied across isotopologues and
the unresolved methyl-methyl coupling is fixed at zero. Results are
conditional numerical candidates; inspect the figures and flags.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from zulf_model.io import load_average
from zulf_model.render import Acquisition
from zulf_model.solver import ObservedSpectrum, ParameterPolicy, Parameterization, RefineSettings, refine
from zulf_model.spinsystem import Component, Interpretation, SpinSystem
from zulf_model.timing import Timer


def methine(jch=135.0, jhh=6.5, jcm=-3.0):
    g = np.array([[0, jhh, jcm], [jhh, 0, jch], [jcm, jch, 0]])
    return SpinSystem.from_group_couplings(["1H", "1H", "13C"], [6, 1, 1], g)


def methyl(jch=125.0, jhh=6.5, jcmh=-3.0, jcother=2.0):
    g = np.zeros((4, 4))
    g[0, 2] = g[2, 0] = jhh
    g[1, 2] = g[2, 1] = jhh
    g[0, 3] = g[3, 0] = jch
    g[1, 3] = g[3, 1] = jcother
    g[2, 3] = g[3, 2] = jcmh
    return SpinSystem.from_group_couplings(["1H", "1H", "1H", "13C"], [3, 3, 1, 1], g)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("average")
    parser.add_argument("ini")
    parser.add_argument("output")
    parser.add_argument("--start-s", type=float, default=0.1)
    parser.add_argument("--sg-window", type=int, default=801)
    parser.add_argument("--background-order", type=int, default=2)
    parser.add_argument("--budget-s", type=float, default=1800.0)
    parser.add_argument("--ranges", default="110,150;230,275")
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    exp = load_average(args.average, args.ini)
    acq = Acquisition(exp.sampling_rate_hz, exp.points, start_sample=int(round(args.start_s * exp.sampling_rate_hz)),
                      sg_window=args.sg_window, sg_order=2, remove_mean=True)
    ranges = [tuple(float(v) for v in r.split(",")) for r in args.ranges.split(";")]
    obs = ObservedSpectrum.from_fid(exp.fid, acq, ranges)
    candidate = Interpretation((Component(methine(), 1.0, "methine-13C"), Component(methyl(), 2.0, "methyl-13C")))
    policy = ParameterPolicy(coupling_margin_hz=10.0, coupling_margin_relative=0.05, rate_bounds_per_s=(0.05, 20.0))
    param = Parameterization.from_interpretation(candidate, policy)
    param.tie("c0.J0-1", "c1.J0-2", "c1.J1-2").fix("c1.J0-1")
    timer = Timer()
    settings = RefineSettings(starts=1, background_order=args.background_order, max_seconds=args.budget_s,
                              max_evaluations=100000, policy=policy)
    result = refine(candidate, obs, settings, parameterization=param, timer=timer)
    summary = result.summary()
    summary.update(acquisition=acq.to_dict(), source=exp.source, timing=timer.report())
    (out / "result.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    np.savez(out / "spectra.npz", frequency_hz=obs.frequencies_hz, observed=obs.values, model=result.prediction,
             components=np.array(result.component_spectra), band=obs.band_index)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(ranges), 1, figsize=(12, 4 * len(ranges)))
        for axis, band in zip(np.atleast_1d(axes), range(len(ranges))):
            m = obs.band_index == band
            f = obs.frequencies_hz[m]
            axis.plot(f, obs.values[m].real, lw=0.7, label="observed (real)")
            axis.plot(f, result.prediction[m].real, lw=0.7, label="model (real)")
            axis.plot(f, (obs.values[m] - result.prediction[m]).real, lw=0.5, label="residual")
            axis.set_xlabel("Frequency (Hz)")
            axis.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "fit_real.png", dpi=120)
    except ImportError:
        pass
    print(json.dumps({k: summary[k] for k in ("relative_residual", "band_relative_residuals", "flags",
                                               "evaluations", "elapsed_s")}, default=float))
    print(json.dumps({k: round(float(v), 4) for k, v in result.parameters.items()}))


if __name__ == "__main__":
    main()
