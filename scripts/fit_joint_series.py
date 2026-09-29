"""Joint fit of a concentration series: every coupling linear (hence monotonic) in the concentration.

    python scripts/fit_joint_series.py --series series.json --structure '{"motif": "pyridine ring", ...}' \\
        --couplings '{"J(HA2,HA3)": 4.9, ...}' [--prior-sigma-hh 0.5] [--prior-sigma-ch 1.5] [--prior-weight 10]
        [--variant ratios] [--range 140,200] [--starts 4] [--spread 0.5] [--out runs/processed/joint]

series.json: [{"id": ..., "x": 0.5, "freq": F.npy, "values": V.npy, "start": fit.json of a single-spectrum fit
(optional)}, ...]. Model: J_k(x) = a_k + b_k (x - mean x) for every free coupling k (shared by all spectra), and per
spectrum its own decay rates and delay (nonlinear) and its gains, shared phase and background (solved linearly
inside each spectrum's forward model). The residual is the concatenation of the spectra's weighted residuals; the
Jacobian is each spectrum's analytic variable-projection Jacobian mapped by the chain rule (dJ/da = 1,
dJ/db = x - mean x). Optional Gaussian priors on a_k around the starting couplings (as RefineSettings.priors, per
spectrum count: the prior row is scaled by the mean spectrum norm). Writes OUT/fit.json (a, b and their linearised
errors, J at every x), OUT/couplings_vs_x.png and OUT/spectra.png.
"""
import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from zulf_core.physics.protocol import SUDDEN_DROP                     # noqa: E402
from zulf_core.solver import ObservedSpectrum                          # noqa: E402
from zulf_core.solver.forward import MixtureForward                    # noqa: E402
from zulf_core.solver.refine import _signal_kwargs                     # noqa: E402
from zulf_hypothesis.builder import build_model                        # noqa: E402
from zulf_hypothesis.fit import default_fit_base                       # noqa: E402
from zulf_hypothesis.search import _settings_for                       # noqa: E402


class JointSeries:
    """Shared linear couplings over several spectra of one structure."""

    def __init__(self, model, settings, observations, xs):
        self.model, self.settings = model, settings
        self.params = [settings.parameterize(model.interpretation) for _ in observations]
        self.forwards = [MixtureForward(p, o, SUDDEN_DROP, settings.gain_model, settings.background_order,
                                        settings.band_weighting, **_signal_kwargs(settings))
                         for p, o in zip(self.params, observations)]
        self.free = self.params[0].free_names
        self.coupling = [n for n in self.free if self.params[0].parameters[n].kind == "coupling"]
        self.local = [n for n in self.free if n not in self.coupling]
        self.xs = np.asarray(xs, float)
        self.xbar = float(self.xs.mean())
        self.nc, self.nl, self.ns = len(self.coupling), len(self.local), len(observations)
        self.col = {n: i for i, n in enumerate(self.free)}
        self.priors = (np.zeros(0, int), np.zeros(0), np.zeros(0))

    # z = [a (nc), b (nc), local of spectrum 0 (nl), local of spectrum 1, ...]
    def spectrum_vector(self, z, s):
        a, b = z[:self.nc], z[self.nc:2 * self.nc]
        loc = z[2 * self.nc + s * self.nl: 2 * self.nc + (s + 1) * self.nl]
        x = np.empty(len(self.free))
        for i, n in enumerate(self.coupling):
            x[self.col[n]] = a[i] + b[i] * (self.xs[s] - self.xbar)
        for i, n in enumerate(self.local):
            x[self.col[n]] = loc[i]
        return x

    def pack(self, per_spectrum_x):
        """Start vector from per-spectrum free vectors: a = mean of the couplings, b = least-squares slope."""
        X = np.array(per_spectrum_x)
        dx = self.xs - self.xbar
        a = [X[:, self.col[n]].mean() for n in self.coupling]
        b = [float(dx @ (X[:, self.col[n]] - X[:, self.col[n]].mean()) / max(dx @ dx, 1e-12)) for n in self.coupling]
        loc = [X[s, self.col[n]] for s in range(self.ns) for n in self.local]
        return np.r_[a, b, loc]

    def bounds(self, slope_bound):
        lo, hi = self.params[0].bounds()
        la = [lo[self.col[n]] for n in self.coupling]
        ha = [hi[self.col[n]] for n in self.coupling]
        ll = [lo[self.col[n]] for _ in range(self.ns) for n in self.local]
        hl = [hi[self.col[n]] for _ in range(self.ns) for n in self.local]
        return (np.r_[la, [-slope_bound] * self.nc, ll], np.r_[ha, [slope_bound] * self.nc, hl])

    def set_priors(self, index, mean, sigma, weight):
        norm = float(np.mean([f.norm for f in self.forwards]))
        self.priors = (np.asarray(index, int), np.asarray(mean, float),
                       np.sqrt(weight) / np.asarray(sigma, float) / norm)

    def residual(self, z):
        parts = [f.predict(self.spectrum_vector(z, s)).residual for s, f in enumerate(self.forwards)]
        idx, mean, scale = self.priors
        parts.append((z[idx] - mean) * scale)
        return np.concatenate(parts)

    def jacobian(self, z):
        blocks = []
        for s, f in enumerate(self.forwards):
            js = f.jacobian(self.spectrum_vector(z, s))
            out = np.zeros((js.shape[0], len(z)))
            for i, n in enumerate(self.coupling):
                out[:, i] = js[:, self.col[n]]
                out[:, self.nc + i] = js[:, self.col[n]] * (self.xs[s] - self.xbar)
            for i, n in enumerate(self.local):
                out[:, 2 * self.nc + s * self.nl + i] = js[:, self.col[n]]
            blocks.append(out)
        idx, _, scale = self.priors
        rows = np.zeros((len(idx), len(z)))
        rows[np.arange(len(idx)), idx] = scale
        blocks.append(rows)
        return np.vstack(blocks)

    def data_residuals(self, z):
        out = []
        for s, f in enumerate(self.forwards):
            pred = f.predict(self.spectrum_vector(z, s))
            m = f.data_signal_mask if f.data_signal_mask is not None else np.ones(len(f.y), bool)
            out.append(float(np.linalg.norm(f.mismatch(pred.model, m)) /
                             max(np.linalg.norm(f.mismatch(np.zeros_like(f.y), m)), 1e-30)))
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--prior-sigma-hh", type=float, default=0.0)
    ap.add_argument("--prior-sigma-ch", type=float, default=0.0)
    ap.add_argument("--prior-weight", type=float, default=1.0)
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--starts", type=int, default=4)
    ap.add_argument("--spread", type=float, default=0.5)
    ap.add_argument("--slope-bound", type=float, default=20.0, help="bound on |dJ/dx| in Hz per mole fraction")
    ap.add_argument("--max-nfev", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/processed/joint")
    args = ap.parse_args()
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    series = json.load(open(args.series))
    lo, hi = (float(v) for v in args.range.split(","))
    obs = [ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), np.load(e["values"]).astype(float),
                                          [(lo, hi)], record=None, real_only=True, label=e["id"]) for e in series]
    spec = json.loads(args.structure)
    spec.setdefault("compound", "series")
    fragment = override_couplings(reg.structure_for(spec), json.loads(args.couplings))
    model = build_model(fragment, ranges=[(lo, hi)])
    settings = _settings_for(model, args.variant, default_fit_base())
    joint = JointSeries(model, settings, obs, [e["x"] for e in series])
    centre = joint.params[0].values()
    # per-spectrum starts: a single-spectrum fit when given, else the structure's values
    xs0 = []
    for e, p in zip(series, joint.params):
        x = p.vector()
        if e.get("start"):
            fit = json.load(open(e["start"]))
            names = {**fit.get("stage2", {}).get("couplings", {})}
            for key, value in names.items():
                for n in model.coupling_names.get(key, []):
                    if n in joint.col:
                        x[joint.col[n]] = value
            for n, rate in fit.get("stage2", {}).get("rates_per_s", {}).items():
                if n in joint.col:
                    x[joint.col[n]] = np.log(rate)
            if "phase_delay" in joint.col and "delay_ms" in fit.get("stage2", {}):
                x[joint.col["phase_delay"]] = fit["stage2"]["delay_ms"] * 1e-3
        xs0.append(x)
    z0 = joint.pack(xs0)
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = joint.params[0].ties.get(n, n)
            if leader in joint.coupling and leader not in key_of:
                key_of[leader] = key
    idx, mean, sigma = [], [], []
    for i, n in enumerate(joint.coupling):
        value, key = centre[n], key_of.get(n, n)
        if abs(value) >= 50.0:
            continue
        s = args.prior_sigma_hh if key.startswith("J(H") else args.prior_sigma_ch
        if s > 0:
            idx.append(i)
            mean.append(value)
            sigma.append(s)
    joint.set_priors(idx, mean, sigma, args.prior_weight)
    lower, upper = joint.bounds(args.slope_bound)
    z0 = np.clip(z0, lower + 1e-9, upper - 1e-9)
    rng = np.random.default_rng(args.seed)
    starts = [z0] + [np.clip(z0 + np.r_[rng.normal(0, args.spread, joint.nc), np.zeros(len(z0) - joint.nc)],
                             lower + 1e-9, upper - 1e-9) for _ in range(max(args.starts - 1, 0))]
    t0 = time.time()
    solutions = []
    for k, z in enumerate(starts):
        sol = least_squares(joint.residual, z, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                            max_nfev=args.max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        solutions.append((float(2 * sol.cost), sol.x))
        print(f"start {k}: score {2 * sol.cost:.5f} ({time.time() - t0:.0f} s)", flush=True)
    solutions.sort(key=lambda s: s[0])
    score, z = solutions[0]
    r = joint.residual(z)
    jac = joint.jacobian(z)
    dof = max(len(r) - len(z), 1)
    cov = float(r @ r) / dof * np.linalg.pinv(jac.T @ jac)
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result = {"x": joint.xs.tolist(), "xbar": joint.xbar, "prior": {"sigma_hh": args.prior_sigma_hh,
              "sigma_ch": args.prior_sigma_ch, "weight": args.prior_weight}, "scores": [s for s, _ in solutions],
              "data_region_residuals": dict(zip([e["id"] for e in series], joint.data_residuals(z))),
              "seconds": round(time.time() - t0), "couplings": {}}
    for i, n in enumerate(joint.coupling):
        key = key_of.get(n, n)
        a, b = z[i], z[joint.nc + i]
        result["couplings"][key] = {"a_hz": a, "a_std": se[i], "slope_hz_per_x": b, "slope_std": se[joint.nc + i],
                                    "prior_centre": centre[n] if abs(centre[n]) < 50 else None,
                                    "J_at_x": [a + b * (x - joint.xbar) for x in joint.xs]}
    json.dump(result, open(out / "fit.json", "w"), indent=1)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = list(result["couplings"])
        cols = 5
        rows = int(np.ceil(len(keys) / cols))
        fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.4 * rows))
        for ax, key in zip(np.ravel(axes), keys):
            c = result["couplings"][key]
            ax.plot(joint.xs, c["J_at_x"], "o-", color="#d1495b", ms=3)
            if c["prior_centre"] is not None:
                ax.axhline(c["prior_centre"], color="#999999", ls=":", lw=0.8)
            ax.set_title(f"{key}: {c['slope_hz_per_x']:+.2f} +- {c['slope_std']:.2f} Hz/x", fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in np.ravel(axes)[len(keys):]:
            ax.axis("off")
        fig.supxlabel("pyridine mole fraction")
        fig.tight_layout()
        fig.savefig(out / "couplings_vs_x.png", dpi=100)
        plt.close(fig)
        fig, axes = plt.subplots(len(obs), 1, figsize=(12, 2.3 * len(obs)), sharex=True)
        for s, (ax, o, f) in enumerate(zip(np.atleast_1d(axes), obs, joint.forwards)):
            pred = f.predict(joint.spectrum_vector(z, s)).model.real
            y = o.values.real
            m = np.abs(y).max()
            ax.plot(o.frequencies_hz, y / m, color="#222222", lw=0.8, label="experiment")
            ax.plot(o.frequencies_hz, pred / m, color="#d1495b", lw=0.9, label="simulation")
            ax.plot(o.frequencies_hz, (y - pred) / m - 0.35, color="#999999", lw=0.6, label="residual")
            ax.set_title(f"{series[s]['id']} (x = {joint.xs[s]:.2f}): relative residual "
                         f"{result['data_region_residuals'][series[s]['id']]:.3f}", fontsize=9, loc="left")
            ax.set_yticks([])
        np.atleast_1d(axes)[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
        np.atleast_1d(axes)[-1].set_xlabel("frequency (Hz)")
        fig.tight_layout()
        fig.savefig(out / "spectra.png", dpi=100)
        plt.close(fig)
    except ImportError:
        pass
    print(json.dumps({"score": score, "residuals": result["data_region_residuals"], "seconds": result["seconds"]}))


if __name__ == "__main__":
    main()
