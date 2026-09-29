"""Joint fit of a concentration series: every coupling monotonic in the concentration, direction and shape free.

    python scripts/fit_joint_series.py --series series.json --structure '{"motif": "pyridine ring", ...}' \\
        --couplings '{"J(HA2,HA3)": 4.9, ...}' [--start-couplings '{...}' | --from-joint PREVIOUS/fit.json]
        [--prior-sigma-hh 0.5] [--prior-sigma-ch 1.5] [--prior-weight 10] [--variant ratios] [--range 140,200]
        [--starts 6] [--spread 0.5] [--out runs/processed/joint]

series.json (increasing x): [{"id": ..., "x": 0.5, "freq": F.npy, "values": V.npy, "start": fit.json of a
single-spectrum fit (optional)}, ...].

Model (user: monotonic, never linear): for every free coupling k and spectrum i (x_1 < ... < x_n)

    J_k(x_i) = v_k + A_k c_k,i,   c_k,i = (e^w_k,1 + ... + e^w_k,i-1) / (e^w_k,1 + ... + e^w_k,n-1)

so c rises from 0 to 1 in steps of any size (shape free), and the sign of A_k (fitted) is the direction: J_k is
monotonic by construction and the data choose direction, size and shape. Every spectrum keeps its own decay rates
and delay (nonlinear) and its gains, shared phase and background (solved linearly in its forward model). Residual:
the spectra's weighted residuals, plus optional Gaussian priors on each coupling's series average around the
--couplings values (weight as RefineSettings.priors, scaled by the mean spectrum norm). Jacobian: each spectrum's
analytic variable-projection Jacobian, chain rule through J_k(x_i). Outputs OUT/fit.json (J at every x with
linearised errors, direction, total change), OUT/J_table.csv, OUT/couplings_vs_x.png and OUT/spectra.png.
"""
import argparse
import csv
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

W_BOUND = 8.0          # |w| bound: the smallest step share is about e^-16 of the largest, i.e. effectively zero


def monotone_profile(w: np.ndarray):
    """c (n,) rising from 0 to 1 and dc/dw (n, n-1) for the softmax step shares e^w / sum e^w."""
    e = np.exp(w - w.max())
    total = e.sum()
    share = e / total
    n = len(w) + 1
    below = (np.arange(n)[:, None] > np.arange(n - 1)[None, :]).astype(float)      # step j lies below point i
    c = below @ share
    dc = below * share[None, :] - c[:, None] * share[None, :]
    return c, dc


class JointSeries:
    """Shared monotonic couplings over several spectra of one structure (see the module docstring)."""

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
        if len(self.xs) < 2 or np.any(np.diff(self.xs) <= 0):
            raise ValueError("Give at least two spectra in increasing concentration.")
        self.nc, self.nl, self.ns = len(self.coupling), len(self.local), len(observations)
        self.m = self.ns + 1                      # v, A, w_1 .. w_{n-1}
        self.nt = self.nc * self.m
        self.col = {n: i for i, n in enumerate(self.free)}
        self.priors = ([], np.zeros(0), np.zeros(0))

    def theta(self, z, k):
        return z[k * self.m:(k + 1) * self.m]

    def coupling_values(self, z, k):
        th = self.theta(z, k)
        c, _ = monotone_profile(th[2:])
        return th[0] + th[1] * c

    def coupling_jacobian(self, z, k):
        """d J_k(x_i) / d theta_k: (n spectra, m)."""
        th = self.theta(z, k)
        c, dc = monotone_profile(th[2:])
        return np.column_stack([np.ones(self.ns), c, th[1] * dc])

    # z = [theta_0, ..., theta_{nc-1}, local of spectrum 0 (nl), local of spectrum 1, ...]
    def spectrum_vector(self, z, s, values=None):
        values = values if values is not None else [self.coupling_values(z, k) for k in range(self.nc)]
        x = np.empty(len(self.free))
        for k, n in enumerate(self.coupling):
            x[self.col[n]] = values[k][s]
        loc = z[self.nt + s * self.nl: self.nt + (s + 1) * self.nl]
        for i, n in enumerate(self.local):
            x[self.col[n]] = loc[i]
        return x

    def pack(self, table, per_spectrum_x):
        """Start vector from a (spectra x couplings) table: v = first value, A = last - first, step shares from the
        table's steps in the direction of A (equal shares where the table is flat or not monotonic)."""
        z = []
        for k in range(self.nc):
            col = table[:, k]
            a = col[-1] - col[0]
            steps = np.diff(col) * (1.0 if a >= 0 else -1.0)
            steps = np.maximum(steps, 0.0)
            w = np.log(steps + 1e-3 * max(steps.max(), 1e-3)) if steps.sum() > 0 else np.zeros(self.ns - 1)
            z.extend([col[0], a] + list(np.clip(w - w.max(), -W_BOUND, W_BOUND)))
        X = np.array(per_spectrum_x)
        z.extend(X[s, self.col[n]] for s in range(self.ns) for n in self.local)
        return np.array(z, float)

    def bounds(self, change_bound):
        lo, hi = self.params[0].bounds()
        lz, hz = [], []
        for n in self.coupling:
            lz += [lo[self.col[n]], -change_bound] + [-W_BOUND] * (self.ns - 1)
            hz += [hi[self.col[n]], change_bound] + [W_BOUND] * (self.ns - 1)
        lz.extend(lo[self.col[n]] for _ in range(self.ns) for n in self.local)
        hz.extend(hi[self.col[n]] for _ in range(self.ns) for n in self.local)
        return np.array(lz), np.array(hz)

    def set_priors(self, couplings, mean, sigma, weight):
        """Gaussian priors on the series-average value of couplings (indices into self.coupling)."""
        norm = float(np.mean([f.norm for f in self.forwards]))
        self.priors = (list(couplings), np.asarray(mean, float), np.sqrt(weight) / np.asarray(sigma, float) / norm)

    def residual(self, z):
        values = [self.coupling_values(z, k) for k in range(self.nc)]
        parts = [f.predict(self.spectrum_vector(z, s, values)).residual for s, f in enumerate(self.forwards)]
        ks, mean, scale = self.priors
        if ks:
            parts.append((np.array([values[k].mean() for k in ks]) - mean) * scale)
        return np.concatenate(parts)

    def jacobian(self, z):
        values = [self.coupling_values(z, k) for k in range(self.nc)]
        dvals = [self.coupling_jacobian(z, k) for k in range(self.nc)]
        blocks = []
        for s, f in enumerate(self.forwards):
            js = f.jacobian(self.spectrum_vector(z, s, values))
            out = np.zeros((js.shape[0], len(z)))
            for k, n in enumerate(self.coupling):
                out[:, k * self.m:(k + 1) * self.m] = np.outer(js[:, self.col[n]], dvals[k][s])
            for i, n in enumerate(self.local):
                out[:, self.nt + s * self.nl + i] = js[:, self.col[n]]
            blocks.append(out)
        ks, _, scale = self.priors
        if ks:
            rows = np.zeros((len(ks), len(z)))
            for r, k in enumerate(ks):
                rows[r, k * self.m:(k + 1) * self.m] = dvals[k].mean(axis=0) * scale[r]
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
    ap.add_argument("--couplings", default="{}", help="JSON {key: Hz}: structure values and prior centres")
    ap.add_argument("--start-couplings", default="", help="JSON {key: Hz}: common starting couplings for every "
                    "spectrum (the prior centres stay the --couplings values)")
    ap.add_argument("--from-joint", default="", help="fit.json of an earlier joint fit: its J at every x as start")
    ap.add_argument("--prior-sigma-hh", type=float, default=0.0)
    ap.add_argument("--prior-sigma-ch", type=float, default=0.0)
    ap.add_argument("--prior-weight", type=float, default=1.0)
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--starts", type=int, default=6)
    ap.add_argument("--spread", type=float, default=0.5, help="perturbation of every coupling level (Hz)")
    ap.add_argument("--change-bound", type=float, default=20.0, help="bound on |total change| of a coupling (Hz)")
    ap.add_argument("--signal-threshold", type=float, default=0.0,
                    help="peak-core threshold of the signal weighting in noise sigma (default: fit base, 4)")
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
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    settings = _settings_for(model, args.variant, base)
    joint = JointSeries(model, settings, obs, [e["x"] for e in series])
    centre = joint.params[0].values()
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = joint.params[0].ties.get(n, n)
            if leader in joint.coupling and leader not in key_of:
                key_of[leader] = key
    # per-spectrum starts: a single-spectrum fit when given, else the structure's values
    xs0 = []
    for e, p in zip(series, joint.params):
        x = p.vector()
        if e.get("start"):
            fit = json.load(open(e["start"]))
            for key, value in fit.get("stage2", {}).get("couplings", {}).items():
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
    if args.from_joint:
        previous = json.load(open(args.from_joint))
        table = np.array([[previous["couplings"][key_of[n]]["J_at_x"][s] for n in joint.coupling]
                          for s in range(joint.ns)])
    z0 = joint.pack(table, xs0)
    ks, mean, sigma = [], [], []
    for k, n in enumerate(joint.coupling):
        value, key = centre[n], key_of.get(n, n)
        if abs(value) >= 50.0:
            continue
        s = args.prior_sigma_hh if key.startswith("J(H") else args.prior_sigma_ch
        if s > 0:
            ks.append(k)
            mean.append(value)
            sigma.append(s)
    joint.set_priors(ks, mean, sigma, args.prior_weight)
    lower, upper = joint.bounds(args.change_bound)
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
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result = {"shape": "monotone (direction and shape free)", "x": joint.xs.tolist(),
              "signal_threshold": settings.signal_threshold,
              "prior": {"sigma_hh": args.prior_sigma_hh, "sigma_ch": args.prior_sigma_ch, "weight": args.prior_weight},
              "scores": [s for s, _ in solutions],
              "data_region_residuals": dict(zip([e["id"] for e in series], joint.data_residuals(z))),
              "seconds": round(time.time() - t0), "couplings": {}}
    for k, n in enumerate(joint.coupling):
        block = slice(k * joint.m, (k + 1) * joint.m)
        values = joint.coupling_values(z, k)
        g = joint.coupling_jacobian(z, k)
        value_cov = g @ cov[block, block] @ g.T
        result["couplings"][key_of.get(n, n)] = {
            "J_at_x": values.tolist(), "J_std_at_x": np.sqrt(np.maximum(np.diag(value_cov), 0)).tolist(),
            "direction": "increasing" if z[block][1] > 0 else "decreasing", "change_hz": float(z[block][1]),
            "prior_centre": centre[n] if abs(centre[n]) < 50 else None}
    json.dump(result, open(out / "fit.json", "w"), indent=1)
    with open(out / "J_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["coupling", "direction"] + [f"J at x={x:.2f} (Hz)" for x in joint.xs] +
                   [f"std at x={x:.2f} (Hz)" for x in joint.xs] + ["prior centre (Hz)"])
        for key, c in result["couplings"].items():
            w.writerow([key, c["direction"]] + [round(v, 3) for v in c["J_at_x"]] +
                       [round(v, 3) for v in c["J_std_at_x"]] + [c["prior_centre"]])
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
            ax.set_title(f"{key}: change {c['change_hz']:+.2f} Hz", fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in np.ravel(axes)[len(keys):]:
            ax.axis("off")
        fig.supxlabel("mole fraction")
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
