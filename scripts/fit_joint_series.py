"""Joint fit of a concentration series: every coupling linear, or only monotonic, in the concentration.

    python scripts/fit_joint_series.py --series series.json --structure '{"motif": "pyridine ring", ...}' \\
        --couplings '{"J(HA2,HA3)": 4.9, ...}' [--prior-sigma-hh 0.5] [--prior-sigma-ch 1.5] [--prior-weight 10]
        [--variant ratios] [--range 140,200] [--starts 4] [--spread 0.5] [--out runs/processed/joint]

--shape monotone --from-joint LINEAR/fit.json: each coupling takes any value at every concentration but moves in one
direction only (steps >= 0 times the direction of the linear fit's slope); priors act on the series average.
series.json (increasing x): [{"id": ..., "x": 0.5, "freq": F.npy, "values": V.npy, "start": fit.json of a single-spectrum fit
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
    """Shared couplings over several spectra of one structure: coupling k at spectrum s is G_k[s] @ theta_k (a
    per-coupling design matrix: linear [1, x - mean x], or monotone [1, sign * step indicators]); every spectrum keeps
    its own decay rates and delay."""

    def __init__(self, model, settings, observations, xs, shape="linear", signs=None):
        self.model, self.settings = model, settings
        self.params = [settings.parameterize(model.interpretation) for _ in observations]
        self.forwards = [MixtureForward(p, o, SUDDEN_DROP, settings.gain_model, settings.background_order,
                                        settings.band_weighting, **_signal_kwargs(settings))
                         for p, o in zip(self.params, observations)]
        self.free = self.params[0].free_names
        self.coupling = [n for n in self.free if self.params[0].parameters[n].kind == "coupling"]
        self.local = [n for n in self.free if n not in self.coupling]
        self.xs = np.asarray(xs, float)
        if np.any(np.diff(self.xs) <= 0):
            raise ValueError("Give the series in increasing concentration.")
        self.xbar = float(self.xs.mean())
        self.nc, self.nl, self.ns = len(self.coupling), len(self.local), len(observations)
        self.col = {n: i for i, n in enumerate(self.free)}
        self.shape = shape
        signs = np.ones(self.nc) if signs is None else np.asarray(signs, float)
        if shape == "linear":
            self.G = [np.column_stack([np.ones(self.ns), self.xs - self.xbar]) for _ in range(self.nc)]
        elif shape == "monotone":
            steps = (np.arange(self.ns)[:, None] > np.arange(self.ns - 1)[None, :]).astype(float)
            self.G = [np.column_stack([np.ones(self.ns), sg * steps]) for sg in signs]
        else:
            raise ValueError("shape must be 'linear' or 'monotone'.")
        self.m = self.G[0].shape[1]
        self.nt = self.nc * self.m
        self.priors = (np.zeros((0, 0)), np.zeros(0), np.zeros(0))

    def theta(self, z, k):
        return z[k * self.m:(k + 1) * self.m]

    def coupling_values(self, z, k):
        return self.G[k] @ self.theta(z, k)

    # z = [theta_0, ..., theta_{nc-1}, local of spectrum 0 (nl), local of spectrum 1, ...]
    def spectrum_vector(self, z, s):
        x = np.empty(len(self.free))
        for k, n in enumerate(self.coupling):
            x[self.col[n]] = self.G[k][s] @ self.theta(z, k)
        loc = z[self.nt + s * self.nl: self.nt + (s + 1) * self.nl]
        for i, n in enumerate(self.local):
            x[self.col[n]] = loc[i]
        return x

    def pack(self, coupling_tables, per_spectrum_x):
        """Start vector: theta_k by least squares from a (spectra x couplings) table (monotone: steps clipped at 0),
        locals from the per-spectrum free vectors."""
        z = []
        for k in range(self.nc):
            th = np.linalg.lstsq(self.G[k], coupling_tables[:, k], rcond=None)[0]
            if self.shape == "monotone":
                th[1:] = np.maximum(th[1:], 0.0)
            z.extend(th)
        X = np.array(per_spectrum_x)
        z.extend(X[s, self.col[n]] for s in range(self.ns) for n in self.local)
        return np.array(z, float)

    def bounds(self, slope_bound, step_bound):
        lo, hi = self.params[0].bounds()
        lz, hz = [], []
        for n in self.coupling:
            lz.append(lo[self.col[n]])
            hz.append(hi[self.col[n]])
            if self.shape == "linear":
                lz.append(-slope_bound)
                hz.append(slope_bound)
            else:
                lz.extend([0.0] * (self.m - 1))
                hz.extend([step_bound] * (self.m - 1))
        lz.extend(lo[self.col[n]] for _ in range(self.ns) for n in self.local)
        hz.extend(hi[self.col[n]] for _ in range(self.ns) for n in self.local)
        return np.array(lz), np.array(hz)

    def set_priors(self, couplings, mean, sigma, weight):
        """Gaussian priors on the series-average value of couplings (indices into self.coupling)."""
        norm = float(np.mean([f.norm for f in self.forwards]))
        rows = np.zeros((len(couplings), len(self.coupling) * self.m + self.ns * self.nl))
        for r, k in enumerate(couplings):
            rows[r, k * self.m:(k + 1) * self.m] = self.G[k].mean(axis=0)
        self.priors = (rows, np.asarray(mean, float), np.sqrt(weight) / np.asarray(sigma, float) / norm)

    def residual(self, z):
        parts = [f.predict(self.spectrum_vector(z, s)).residual for s, f in enumerate(self.forwards)]
        rows, mean, scale = self.priors
        if len(mean):
            parts.append((rows @ z - mean) * scale)
        return np.concatenate(parts)

    def jacobian(self, z):
        blocks = []
        for s, f in enumerate(self.forwards):
            js = f.jacobian(self.spectrum_vector(z, s))
            out = np.zeros((js.shape[0], len(z)))
            for k, n in enumerate(self.coupling):
                out[:, k * self.m:(k + 1) * self.m] = np.outer(js[:, self.col[n]], self.G[k][s])
            for i, n in enumerate(self.local):
                out[:, self.nt + s * self.nl + i] = js[:, self.col[n]]
            blocks.append(out)
        rows, mean, scale = self.priors
        if len(mean):
            blocks.append(rows * scale[:, None])
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
    ap.add_argument("--slope-bound", type=float, default=20.0, help="linear: bound on |dJ/dx| (Hz per x)")
    ap.add_argument("--shape", default="linear", choices=["linear", "monotone"])
    ap.add_argument("--step-bound", type=float, default=10.0, help="monotone: bound on each step between spectra")
    ap.add_argument("--start-couplings", default="", help="JSON {key: Hz}: starting couplings for every spectrum "
                    "(overrides the per-spectrum starts; the prior centres stay the --couplings values)")
    ap.add_argument("--from-joint", default="", help="fit.json of an earlier joint fit: starting couplings and, "
                    "for the monotone shape, the direction of every coupling (sign of its linear slope)")
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
    previous = json.load(open(args.from_joint)) if args.from_joint else None
    probe = settings.parameterize(model.interpretation)
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = probe.ties.get(n, n)
            if leader in probe.free_names and probe.parameters[leader].kind == "coupling" and leader not in key_of:
                key_of[leader] = key
    coupling_names = [n for n in probe.free_names if probe.parameters[n].kind == "coupling"]
    signs = None
    if args.shape == "monotone":
        if previous is None:
            raise SystemExit("--shape monotone needs --from-joint (directions from a linear joint fit).")
        signs = [1.0 if previous["couplings"][key_of[n]]["slope_hz_per_x"] >= 0 else -1.0 for n in coupling_names]
    joint = JointSeries(model, settings, obs, [e["x"] for e in series], shape=args.shape, signs=signs)
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
    if args.start_couplings:
        for key, value in json.loads(args.start_couplings).items():
            for n in model.coupling_names.get(key, []):
                if n in joint.col:
                    for x in xs0:
                        x[joint.col[n]] = value
    table = np.array([[x[joint.col[n]] for n in joint.coupling] for x in xs0])
    if previous is not None:
        table = np.array([[previous["couplings"][key_of[n]]["J_at_x"][s] for n in joint.coupling]
                          for s in range(joint.ns)])
    z0 = joint.pack(table, xs0)
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
    lower, upper = joint.bounds(args.slope_bound, args.step_bound)
    z0 = np.clip(z0, lower + 1e-9, upper - 1e-9)
    rng = np.random.default_rng(args.seed)
    level = np.zeros(len(z0))
    level[:joint.nt:joint.m] = 1.0                        # perturb the level of every coupling, not its shape
    starts = [z0] + [np.clip(z0 + level * rng.normal(0, args.spread, len(z0)), lower + 1e-9, upper - 1e-9)
                     for _ in range(max(args.starts - 1, 0))]
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
    result = {"shape": args.shape, "x": joint.xs.tolist(), "xbar": joint.xbar, "prior": {"sigma_hh": args.prior_sigma_hh,
              "sigma_ch": args.prior_sigma_ch, "weight": args.prior_weight}, "scores": [s for s, _ in solutions],
              "data_region_residuals": dict(zip([e["id"] for e in series], joint.data_residuals(z))),
              "seconds": round(time.time() - t0), "couplings": {}}
    for k, n in enumerate(joint.coupling):
        key = key_of.get(n, n)
        block = slice(k * joint.m, (k + 1) * joint.m)
        values = joint.coupling_values(z, k)
        value_cov = joint.G[k] @ cov[block, block] @ joint.G[k].T
        entry = {"J_at_x": values.tolist(), "J_std_at_x": np.sqrt(np.maximum(np.diag(value_cov), 0)).tolist(),
                 "prior_centre": centre[n] if abs(centre[n]) < 50 else None}
        if joint.shape == "linear":
            entry.update(a_hz=z[block][0], a_std=se[block][0], slope_hz_per_x=z[block][1], slope_std=se[block][1])
        else:
            entry.update(direction="increasing" if joint.G[k][-1, 1:].sum() > 0 else "decreasing",
                         steps_hz=z[block][1:].tolist(), change_hz=float(values[-1] - values[0]))
        result["couplings"][key] = entry
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
            ax.errorbar(joint.xs, c["J_at_x"], yerr=c["J_std_at_x"], fmt="o-", color="#d1495b", ms=3, lw=1,
                        capsize=2)
            if c["prior_centre"] is not None:
                ax.axhline(c["prior_centre"], color="#999999", ls=":", lw=0.8)
            label = (f"{c['slope_hz_per_x']:+.2f} +- {c['slope_std']:.2f} Hz/x" if joint.shape == "linear"
                     else f"change {c['change_hz']:+.2f} Hz")
            ax.set_title(f"{key}: {label}", fontsize=8)
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
