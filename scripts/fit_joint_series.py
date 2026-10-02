"""Joint fit of a concentration series: every coupling monotonic in the concentration, direction and shape free.

    python scripts/fit_joint_series.py --series series.json --structure '{"motif": "pyridine ring", ...}' \\
        --couplings '{"J(HA2,HA3)": 4.9, ...}' [--start-couplings '{...}' | --from-joint PREVIOUS/fit.json]
        [--prior-sigma-hh 0.5] [--prior-sigma-ch 1.5] [--prior-weight 10] [--variant ratios] [--range 140,200]
        [--starts 6] [--spread 0.5] [--out runs/processed/joint]

series.json (increasing x): [{"id": ..., "x": 0.5, "freq": F.npy, "values": V.npy, "start": fit.json of a
single-spectrum fit (optional)}, ...]; an entry may instead give "fid": raw averaged FID (.npy; processed with the
confirmed-sample recipe of configs/confirmed_samples.json and fitted as complex data) and "ranges".

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
from scipy.special import expit
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
        if len(self.xs) < 2 or np.any(np.diff(self.xs) < 0):
            raise ValueError("Give the spectra in non-decreasing concentration.")
        # concentration nodes: spectra at the same concentration share one set of couplings
        self.nodes = np.unique(self.xs)
        if len(self.nodes) < 2:
            raise ValueError("Need at least two different concentrations.")
        self.node_of = [int(np.searchsorted(self.nodes, x)) for x in self.xs]
        self.nc, self.nl, self.ns = len(self.coupling), len(self.local), len(observations)
        self.nn = len(self.nodes)
        self.m = self.nn + 1                      # v, A, w_1 .. w_{nodes-1}
        self.nt = self.nc * self.m
        self.col = {n: i for i, n in enumerate(self.free)}
        self.priors = ([], np.zeros(0), np.zeros(0))
        self.smoothing = None                     # per spectrum: residual transform S = W K W^-1 (or None)
        self.peaks = None                         # per spectrum: index windows around data peak tops (or None)
        self.peak_strength = 0.0

    def set_peak_penalty(self, strength, prominence=0.12, tolerance_hz=0.15, min_sigma=2.0, smooth=0.0):
        """Missing-peak rows: for every data peak top (prominence >= `prominence` times the largest abs value and
        >= `min_sigma` noise sigma) one extra residual row strength * min(max over +-tolerance_hz of r, 0),
        r = W (model - data) / norm the ordinary residual. The row is zero while the model reaches the data peak
        somewhere near the top and grows when the model is too low there (peak missing or too weak); the
        ordinary residual and its weights are unchanged. Data only (peak tops are fixed). strength 0 = off.

        smooth > 0 replaces max and min(., 0) by smooth versions for the optimiser (the hard rows have kinks,
        which keep least squares from converging): soft max tau log sum exp(r / tau) and soft hinge
        -s log(1 + exp(-m / s)), tau = s = smooth times the peak's own height in residual units. Scores for
        ranking should use the hard rows (smooth 0)."""
        from scipy.signal import find_peaks
        from zulf_core.solver.forward import signal_regions
        self.peak_strength = float(strength)
        self.peak_smooth = float(smooth)
        if not strength:
            self.peaks = None
            return
        self.peaks, self.peak_heights = [], []
        for f in self.forwards:
            y = f.y.real
            _, sigma = signal_regions(f.f, f.y, f.band, 4.0)
            tops, _ = find_peaks(y, prominence=max(prominence * float(np.abs(f.y).max()), min_sigma * sigma))
            spacing = float(np.median(np.diff(f.f)))
            half = max(int(round(tolerance_hz / spacing)), 0)
            self.peaks.append([np.arange(max(t - half, 0), min(t + half + 1, len(y))) for t in tops])
            self.peak_heights.append(np.abs(y[tops]) * f.weight[tops] / f.norm)

    def _peak_rows(self, s, r):
        """Rows of spectrum s from its ordinary residual r, and per row (indices, coefficients) of
        d row / d r (the Jacobian row is the matching combination of ordinary rows)."""
        rows, grads = [], []
        for w, h in zip(self.peaks[s], self.peak_heights[s]):
            if self.peak_smooth > 0:
                t = self.peak_smooth * h
                a = r[w] / t
                p = np.exp(a - a.max())
                m = t * (a.max() + np.log(p.sum()))
                p /= p.sum()
                rows.append(-self.peak_strength * t * np.logaddexp(0.0, -m / t))
                grads.append((w, self.peak_strength * expit(-m / t) * p))
            else:
                j = w[int(np.argmax(r[w]))]
                v = min(float(r[j]), 0.0)
                rows.append(self.peak_strength * v)
                grads.append((np.array([j]), np.array([self.peak_strength])) if v < 0 else (np.zeros(0, int), np.zeros(0)))
        return np.array(rows), grads

    def set_smoothing(self, sigma_hz):
        """Coarse-to-fine continuation: compare Gaussian-smoothed spectra. The same kernel K (sigma_hz, over the
        points' frequencies, rows normalised) acts on data and model, so K (model - data) is compared; in residual
        units r = W (model - data) / norm this is S r with S = W K W^-1 (block-diagonal for complex data).
        sigma_hz <= 0 switches smoothing off."""
        if not sigma_hz or sigma_hz <= 0:
            self.smoothing = None
            return
        mats = []
        for f in self.forwards:
            d = f.f[:, None] - f.f[None, :]
            k = np.exp(-0.5 * (d / sigma_hz) ** 2)
            k /= k.sum(axis=1, keepdims=True)
            w = f.weight / f.norm
            s = (w[:, None] * k) / w[None, :]
            mats.append(s if f.real_only else np.block([[s, np.zeros_like(s)], [np.zeros_like(s), s]]))
        self.smoothing = mats

    def _smooth(self, s, vec_or_mat):
        return vec_or_mat if self.smoothing is None else self.smoothing[s] @ vec_or_mat

    def theta(self, z, k):
        return z[k * self.m:(k + 1) * self.m]

    def coupling_values(self, z, k):
        th = self.theta(z, k)
        c, _ = monotone_profile(th[2:])
        return th[0] + th[1] * c

    def coupling_jacobian(self, z, k):
        """d J_k(node) / d theta_k: (nodes, m)."""
        th = self.theta(z, k)
        c, dc = monotone_profile(th[2:])
        return np.column_stack([np.ones(self.nn), c, th[1] * dc])

    # z = [theta_0, ..., theta_{nc-1}, local of spectrum 0 (nl), local of spectrum 1, ...]
    def spectrum_vector(self, z, s, values=None):
        values = values if values is not None else [self.coupling_values(z, k) for k in range(self.nc)]
        x = np.empty(len(self.free))
        for k, n in enumerate(self.coupling):
            x[self.col[n]] = values[k][self.node_of[s]]
        loc = z[self.nt + s * self.nl: self.nt + (s + 1) * self.nl]
        for i, n in enumerate(self.local):
            x[self.col[n]] = loc[i]
        return x

    def pack(self, table, per_spectrum_x):
        """Start vector from a (nodes x couplings) table: v = first value, A = last - first, step shares from the
        table's steps in the direction of A (equal shares where the table is flat or not monotonic)."""
        z = []
        for k in range(self.nc):
            col = table[:, k]
            a = col[-1] - col[0]
            steps = np.diff(col) * (1.0 if a >= 0 else -1.0)
            steps = np.maximum(steps, 0.0)
            w = np.log(steps + 1e-3 * max(steps.max(), 1e-3)) if steps.sum() > 0 else np.zeros(self.nn - 1)
            z.extend([col[0], a] + list(np.clip(w - w.max(), -W_BOUND, W_BOUND)))
        X = np.array(per_spectrum_x)
        z.extend(X[s, self.col[n]] for s in range(self.ns) for n in self.local)
        return np.array(z, float)

    def bounds(self, change_bound):
        lo, hi = self.params[0].bounds()
        lz, hz = [], []
        for n in self.coupling:
            lz += [lo[self.col[n]], -change_bound] + [-W_BOUND] * (self.nn - 1)
            hz += [hi[self.col[n]], change_bound] + [W_BOUND] * (self.nn - 1)
        lz.extend(lo[self.col[n]] for _ in range(self.ns) for n in self.local)
        hz.extend(hi[self.col[n]] for _ in range(self.ns) for n in self.local)
        return np.array(lz), np.array(hz)

    def set_priors(self, couplings, mean, sigma, weight):
        """Gaussian priors on the series-average value of couplings (indices into self.coupling)."""
        norm = float(np.mean([f.norm for f in self.forwards]))
        self.priors = (list(couplings), np.asarray(mean, float), np.sqrt(weight) / np.asarray(sigma, float) / norm)

    def residual(self, z):
        values = [self.coupling_values(z, k) for k in range(self.nc)]
        parts = []
        for s, f in enumerate(self.forwards):
            r = f.predict(self.spectrum_vector(z, s, values)).residual
            parts.append(self._smooth(s, r))
            if self.peaks is not None:
                parts.append(self._peak_rows(s, r)[0])
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
                out[:, k * self.m:(k + 1) * self.m] = np.outer(js[:, self.col[n]], dvals[k][self.node_of[s]])
            for i, n in enumerate(self.local):
                out[:, self.nt + s * self.nl + i] = js[:, self.col[n]]
            blocks.append(self._smooth(s, out))
            if self.peaks is not None:
                _, grads = self._peak_rows(s, f.predict(self.spectrum_vector(z, s, values)).residual)
                blocks.append(np.array([c @ out[idx] for idx, c in grads]).reshape(len(grads), len(z)))
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


def predict_left_out(joint, z, entry, band, key_of, settings, max_nfev=200):
    """Fit the left-out spectrum with every coupling held between its values at the neighbouring concentrations
    (monotonicity allows nothing else): J_k = J_k(left) + u_k (J_k(right) - J_k(left)), u_k in [0, 1]; rates, delay,
    gains, phase and background free. Returns the relative residual on the data cores and the u_k (0 or 1 means
    the data push the coupling to a bracket end). At the ends of the series there is only one neighbour: no
    prediction."""
    x = float(entry["x"])
    left = [i for i in range(joint.nn) if joint.nodes[i] < x]
    right = [i for i in range(joint.nn) if joint.nodes[i] > x]
    if not left or not right:
        return {"prediction": "none (end of the series: one neighbour only)"}
    i, j = left[-1], right[0]
    obs = ObservedSpectrum.from_spectrum(np.load(entry["freq"]).astype(float), np.load(entry["values"]).astype(float),
                                         [band], record=entry.get("record"), phasing=entry.get("phasing"), real_only=True,
                                         label=entry["id"])
    param = joint.settings.parameterize(joint.model.interpretation)
    fw = MixtureForward(param, obs, SUDDEN_DROP, settings.gain_model, settings.background_order,
                        settings.band_weighting, **_signal_kwargs(settings))
    values = [joint.coupling_values(z, k) for k in range(joint.nc)]
    lo_v = np.array([values[k][i] for k in range(joint.nc)])
    hi_v = np.array([values[k][j] for k in range(joint.nc)])
    si, sj = joint.node_of.index(i), joint.node_of.index(j)
    local0 = 0.5 * (z[joint.nt + si * joint.nl: joint.nt + (si + 1) * joint.nl] +
                    z[joint.nt + sj * joint.nl: joint.nt + (sj + 1) * joint.nl])
    lo_b, hi_b = param.bounds()
    col = joint.col

    def vector(p):
        u, loc = p[:joint.nc], p[joint.nc:]
        xv = np.empty(len(joint.free))
        for k, n in enumerate(joint.coupling):
            xv[col[n]] = lo_v[k] + u[k] * (hi_v[k] - lo_v[k])
        for m, n in enumerate(joint.local):
            xv[col[n]] = loc[m]
        return xv

    def residual(p):
        return fw.predict(vector(p)).residual

    def jacobian(p):
        js = fw.jacobian(vector(p))
        out = np.zeros((js.shape[0], len(p)))
        for k, n in enumerate(joint.coupling):
            out[:, k] = js[:, col[n]] * (hi_v[k] - lo_v[k])
        for m, n in enumerate(joint.local):
            out[:, joint.nc + m] = js[:, col[n]]
        return out

    lower = np.r_[np.zeros(joint.nc), [lo_b[col[n]] for n in joint.local]]
    upper = np.r_[np.ones(joint.nc), [hi_b[col[n]] for n in joint.local]]
    p0 = np.clip(np.r_[np.full(joint.nc, 0.5), local0], lower + 1e-9, upper - 1e-9)
    sol = least_squares(residual, p0, jac=jacobian, bounds=(lower, upper), x_scale="jac", max_nfev=max_nfev,
                        ftol=1e-10, xtol=1e-10, gtol=1e-10)
    pred = fw.predict(vector(sol.x))
    m = fw.data_signal_mask if fw.data_signal_mask is not None else np.ones(len(fw.y), bool)
    rel = float(np.linalg.norm(fw.mismatch(pred.model, m)) / max(np.linalg.norm(fw.mismatch(np.zeros_like(fw.y), m)), 1e-30))
    return {"relative_residual": rel, "neighbours_x": [float(joint.nodes[i]), float(joint.nodes[j])],
            "u": {key_of.get(n, n): float(u) for n, u in zip(joint.coupling, sol.x[:joint.nc])},
            "prediction_spectrum": np.asarray(pred.model.real).tolist(),
            "frequencies_hz": obs.frequencies_hz.tolist(), "data": obs.values.real.tolist()}


_JOINT_TASK = None


def _coordinate_scan(joint, z, lower, upper, max_nfev, cycles=2, half_width=4.0, step=0.5, seed=0):
    """Global moves along single couplings: in every cycle each small coupling's level v_k (in random order) is
    set to the best of a grid (current value +- half_width, `step`) with everything else held (objective:
    residual norm, priors included), then all parameters are refined by least squares. Crosses barriers along
    one coordinate that a local fit cannot."""
    rng = np.random.default_rng(seed)
    small = [k for k in range(joint.nc) if abs(z[k * joint.m]) < 50.0]
    for _ in range(cycles):
        for k in rng.permutation(small):
            i = k * joint.m
            grid = np.arange(z[i] - half_width, z[i] + half_width + 1e-9, step)
            grid = grid[(grid > lower[i]) & (grid < upper[i])]
            best_v, best_c = z[i], None
            for v in grid:
                zz = z.copy()
                zz[i] = v
                r = joint.residual(zz)
                c = float(r @ r)
                if best_c is None or c < best_c:
                    best_v, best_c = v, c
            z = z.copy()
            z[i] = best_v
        z = least_squares(joint.residual, z, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                          max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10).x
    return z


def _solve_start(z):
    """One start (module level so worker processes can run it; fork start method): least squares at every level
    of the smoothing schedule in turn (coarse to fine; the last level is always unsmoothed); returns the score of
    the unsmoothed objective."""
    joint, lower, upper, max_nfev, schedule, scan, hold_small = _JOINT_TASK
    if hold_small:
        # stage 1: every small coupling (level, change and shape) held at its start; 1J, rates and delays free
        lo_h, hi_h = lower.copy(), upper.copy()
        for k in range(joint.nc):
            if abs(z[k * joint.m]) < 50.0:
                blk = slice(k * joint.m, (k + 1) * joint.m)
                lo_h[blk] = z[blk] - 1e-9
                hi_h[blk] = z[blk] + 1e-9
        z = least_squares(joint.residual, np.clip(z, lo_h, hi_h), jac=joint.jacobian, bounds=(lo_h, hi_h),
                          x_scale="jac", max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10).x
        z = np.clip(z, lower + 1e-9, upper - 1e-9)
    if scan:
        z = _coordinate_scan(joint, z, lower, upper, max_nfev, **scan)
    for sigma in list(schedule) + [0.0]:
        joint.set_smoothing(sigma)
        sol = least_squares(joint.residual, z, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                            max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        z = sol.x
    return float(2 * sol.cost), sol.x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}", help="JSON {key: Hz}: structure values and prior centres")
    ap.add_argument("--start-couplings", default="", help="JSON {key: Hz}: common starting couplings for every "
                    "spectrum (the prior centres stay the --couplings values)")
    ap.add_argument("--from-joint", default="", help="fit.json of an earlier joint fit: its J at every x as start")
    ap.add_argument("--seeds", default="", help="JSON {\"x\": [...], \"solutions\": [{key: [J at every x]}, ...]}: "
                    "several start tables (e.g. the near-best solutions of earlier fits); each is a start and the "
                    "perturbed starts are spread over them in turn (spectrum parameters as for --from-joint)")
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
    ap.add_argument("--signal-taper", type=float, default=0.0,
                    help="Gaussian fall-off (Hz) of the weight around peak cores (default: fit base, 2 Hz)")
    ap.add_argument("--signal-height-power", type=float, default=0.0,
                    help="extra weight (local peak height / max)^-power: small peaks count more (0 = off)")
    ap.add_argument("--peak-penalty", type=float, default=0.0,
                    help="missing-peak rows: strength per data peak top whose model stays below the data (0 = off)")
    ap.add_argument("--peak-prominence", type=float, default=0.12,
                    help="peak tops: prominence as a fraction of the largest value (and at least 2 noise sigma)")
    ap.add_argument("--peak-tolerance", type=float, default=0.15, help="peak tops: position tolerance (Hz)")
    ap.add_argument("--peak-min-sigma", type=float, default=2.0, help="peak tops: prominence at least this many noise sigma")
    ap.add_argument("--peak-smooth", type=float, default=0.0,
                    help="optimise with smooth missing-peak rows (width = this fraction of each peak's height); "
                         "solutions are then rescored and ranked with the hard rows")
    ap.add_argument("--leave-out", type=int, default=-1, help="index of a spectrum to leave out (leave-one-out)")
    ap.add_argument("--prior-starts", type=int, default=0, help="extra starts drawn from the priors")
    ap.add_argument("--hold-small-first", action="store_true",
                    help="stage 1 of every start: small couplings held, 1J, rates and delays fitted; then all free")
    ap.add_argument("--scan-cycles", type=int, default=0, help="coordinate grid-scan cycles before each fit")
    ap.add_argument("--scan-half-width", type=float, default=4.0)
    ap.add_argument("--scan-step", type=float, default=0.5)
    ap.add_argument("--smoothing", default="", help="coarse-to-fine schedule of Gaussian smoothing widths (Hz), "
                    "e.g. 1.5,0.8,0.4; the unsmoothed fit always ends every start")
    ap.add_argument("--workers", type=int, default=1, help="processes for the starts (set OMP_NUM_THREADS=1)")
    ap.add_argument("--max-nfev", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/processed/joint")
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out), "fit_joint_series")
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    series = json.load(open(args.series))
    left_out = None
    if args.leave_out >= 0:
        left_out = series[args.leave_out]
        series = [e for i, e in enumerate(series) if i != args.leave_out]
    lo, hi = (float(v) for v in args.range.split(","))
    config = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]

    def observation(e):
        """A processed real spectrum (freq, values), or a raw FID ("fid") processed with the confirmed-sample
        recipe and compared as complex data (phase and delay fitted), restricted to e["ranges"] or the band."""
        if "fid" in e:
            path = Path(e["fid"])
            o = reg.observed_for({"file": path.name}, config, str(path.parent))
            return o.restricted([tuple(r) for r in e.get("ranges", [(lo, hi)])])
        return ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), np.load(e["values"]).astype(float),
                                              [(lo, hi)], record=e.get("record"), phasing=e.get("phasing"), real_only=True,
                                              label=e["id"])

    obs = [observation(e) for e in series]
    spec = json.loads(args.structure)
    spec.setdefault("compound", "series")
    fragment = override_couplings(reg.structure_for(spec), json.loads(args.couplings))
    model = build_model(fragment, ranges=[(lo, hi)])
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    if getattr(args, "signal_height_power", 0.0):
        base = dataclasses.replace(base, signal_height_power=args.signal_height_power)
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
        if "delay_ms" in e and "phase_delay" in joint.col:       # starting delay (e.g. the switching edge)
            x[joint.col["phase_delay"]] = float(e["delay_ms"]) * 1e-3
        xs0.append(x)
    if args.start_couplings:
        for key, value in json.loads(args.start_couplings).items():
            for n in model.coupling_names.get(key, []):
                if n in joint.col:
                    for x in xs0:
                        x[joint.col[n]] = value
    first = [joint.node_of.index(i) for i in range(joint.nn)]          # first spectrum of every node
    table = np.array([[xs0[s][joint.col[n]] for n in joint.coupling] for s in first])
    if args.from_joint:
        previous = json.load(open(args.from_joint))
        # the previous fit's row at the same concentration, else the nearest one (a concentration new to the series)
        rows = [int(np.argmin(np.abs(np.asarray(previous["x"]) - x))) for x in joint.nodes]
        table = np.array([[previous["couplings"][key_of[n]]["J_at_x"][s] for n in joint.coupling] for s in rows])
    if args.from_joint:
        # decay rates and delays of spectra already in the previous fit (same id)
        for e, x in zip(series, xs0):
            for n, v in previous.get("spectrum_parameters", {}).get(e["id"], {}).items():
                if n in joint.col:
                    x[joint.col[n]] = v
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
    joint.set_peak_penalty(args.peak_penalty, args.peak_prominence, args.peak_tolerance, min_sigma=args.peak_min_sigma,
                           smooth=args.peak_smooth)
    lower, upper = joint.bounds(args.change_bound)
    z0 = np.clip(z0, lower + 1e-9, upper - 1e-9)
    rng = np.random.default_rng(args.seed)
    level = np.zeros(len(z0))
    level[:joint.nt:joint.m] = 1.0                        # perturb the level of every coupling, not its shape
    centres = [z0]
    if args.seeds:
        seeds = json.load(open(args.seeds))
        rows = [int(np.argmin(np.abs(np.asarray(seeds["x"]) - x))) for x in joint.nodes]
        centres = [np.clip(joint.pack(np.array([[sol[key_of[n]][s] for n in joint.coupling] for s in rows]), xs0),
                           lower + 1e-9, upper - 1e-9) for sol in seeds["solutions"]]
    starts = list(centres) + [
        np.clip(centres[i % len(centres)] + level * rng.normal(0, args.spread, len(z0)), lower + 1e-9, upper - 1e-9)
        for i in range(max(args.starts - len(centres), 0))]
    # random starts drawn from the priors: level of every small coupling ~ N(centre, 2 sigma), direction and size
    # of its change random, shape uniform; 1J levels and the spectrum parameters from the first start
    for _ in range(args.prior_starts):
        z = z0.copy()
        for k, n in enumerate(joint.coupling):
            base_k = k * joint.m
            if k in ks:
                sd = 2.0 * sigma[ks.index(k)]
                z[base_k] = rng.normal(mean[ks.index(k)], sd)
                z[base_k + 1] = rng.normal(0.0, sd)
                z[base_k + 2: base_k + joint.m] = 0.0
        starts.append(np.clip(z, lower + 1e-9, upper - 1e-9))
    t0 = time.time()
    global _JOINT_TASK
    schedule = [float(v) for v in args.smoothing.split(",") if v.strip()] if args.smoothing else []
    scan = ({"cycles": args.scan_cycles, "half_width": args.scan_half_width, "step": args.scan_step,
             "seed": args.seed} if args.scan_cycles else None)
    _JOINT_TASK = (joint, lower, upper, args.max_nfev, schedule, scan, args.hold_small_first)
    if args.workers > 1:
        import multiprocessing
        with multiprocessing.get_context("fork").Pool(args.workers) as pool:
            solved = pool.map(_solve_start, starts, chunksize=1)
    else:
        solved = [_solve_start(z) for z in starts]
    solutions = []
    for k, (score, x) in enumerate(solved):
        solutions.append((score, x))
        print(f"start {k}: score {score:.5f}", flush=True)
    print(f"{len(starts)} starts in {time.time() - t0:.0f} s", flush=True)
    joint.set_smoothing(0.0)
    smooth_scores = None
    if joint.peaks is not None and joint.peak_smooth > 0:
        # rank with the hard missing-peak rows; keep the smooth objective values for the record
        smooth_scores = [sc for sc, _ in solutions]
        joint.peak_smooth = 0.0
        solutions = [(float(np.sum(joint.residual(zz) ** 2)), zz) for _, zz in solutions]
        order = np.argsort([sc for sc, _ in solutions])
        smooth_scores = [smooth_scores[i] for i in order]
    solutions.sort(key=lambda s: s[0])
    score, z = solutions[0]
    r = joint.residual(z)
    jac = joint.jacobian(z)
    dof = max(len(r) - len(z), 1)
    cov = float(r @ r) / dof * np.linalg.pinv(jac.T @ jac)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    start_tables = [{"score": sc, "J_at_x": {key_of.get(n, n): joint.coupling_values(zz, k).tolist()
                                              for k, n in enumerate(joint.coupling)}} for sc, zz in solutions]
    result = {"shape": "monotone (direction and shape free)", "x": joint.nodes.tolist(),
              "spectra": [{"id": e["id"], "x": e["x"]} for e in series],
              "start_solutions": start_tables,
              "signal_threshold": settings.signal_threshold, "signal_taper_hz": settings.signal_taper_hz,
              "signal_height_power": settings.signal_height_power,
              "peak_penalty": {"strength": args.peak_penalty, "prominence": args.peak_prominence,
                               "tolerance_hz": args.peak_tolerance, "min_sigma": args.peak_min_sigma,
                               "tops": [len(t) for t in joint.peaks] if joint.peaks else [],
                               "smooth": args.peak_smooth, "smooth_scores_sorted_by_hard": smooth_scores},
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
    if left_out is not None:
        result["left_out"] = {"id": left_out["id"], "x": left_out["x"],
                              **predict_left_out(joint, z, left_out, (lo, hi), key_of, settings)}
    result["spectrum_parameters"] = {series[si]["id"]: {n: float(z[joint.nt + si * joint.nl + i])
                                                       for i, n in enumerate(joint.local)} for si in range(joint.ns)}
    json.dump(result, open(out / "fit.json", "w"), indent=1)
    with open(out / "J_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["coupling", "direction"] + [f"J at x={x:.3f} (Hz)" for x in joint.nodes] +
                   [f"std at x={x:.3f} (Hz)" for x in joint.nodes] + ["prior centre (Hz)"])
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
            ax.errorbar(joint.nodes, c["J_at_x"], yerr=c["J_std_at_x"], fmt="o-", color="#d1495b", ms=3, lw=1,
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
            pred = f.predict(joint.spectrum_vector(z, s)).model
            y = o.values
            if o.real_only:
                pred, y = pred.real, y.real
            else:                              # complex data (unphased): compare magnitudes in the figure
                pred, y = np.abs(pred), np.abs(y)
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
    run_log.finish({"score": score, "residuals": result["data_region_residuals"],
                    "J_at_x": {k: v["J_at_x"] for k, v in result["couplings"].items()}, "x": result["x"]},
                   figures=["spectra.png", "couplings_vs_x.png", "J_table.csv"])
    print(json.dumps({"score": score, "residuals": result["data_region_residuals"], "seconds": result["seconds"]}))


if __name__ == "__main__":
    main()
