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
and delay (nonlinear) and its gains, shared phase and background (solved linearly in its forward model);
--shared phase_delay makes the delay (or any other listed spectrum parameter) one value for the series;
--shape free drops the monotone constraint (independent J at every concentration; also the mode for one
spectrum); a series entry may give its own fit "ranges" (e.g. to leave out mains harmonics). Residual:
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
from types import SimpleNamespace

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



def _values(entry, real_only=True):
    """Spectrum values of a series entry: real (real_only, the default) or complex."""
    v = np.load(entry["values"])
    return v.real.astype(float) if real_only else v.astype(complex)


def _flag(text):
    """--real-only true|false"""
    t = str(text).strip().lower()
    if t in ("true", "1", "yes"):
        return True
    if t in ("false", "0", "no"):
        return False
    raise ValueError(f"expected true or false, got {text!r}")

def _rate_policy(base, args):
    """Decay-rate families and bounds from the command line: --family-edges splits every isotopologue's
    transitions by frequency (one rate per family, D31 derivatives), --rate-bounds sets the rate range,
    --phase-delay-bounds the range of the fitted delay (ms; an instrument prior, e.g. -4.6,-2.6 for the NMRduino
    setup: near 200 Hz a delay is ambiguous by about 1 / f within one band)."""
    policy = base.policy
    edges = [float(v) for v in getattr(args, "family_edges", "").split(",") if v.strip()]
    if edges:
        if np.any(np.diff(edges) <= 0):
            raise SystemExit("--family-edges must increase")
        policy = dataclasses.replace(policy, family_edges_hz=tuple(edges))
    if getattr(args, "rate_bounds", ""):
        lo, hi = (float(v) for v in args.rate_bounds.split(","))
        policy = dataclasses.replace(policy, rate_bounds_per_s=(lo, hi),
                                     initial_rate_per_s=float(np.clip(policy.initial_rate_per_s, lo, hi)))
    if getattr(args, "phase_delay_bounds", ""):
        lo, hi = (1e-3 * float(v) for v in args.phase_delay_bounds.split(","))
        if not lo < hi:
            raise SystemExit("--phase-delay-bounds: lo < hi (ms)")
        policy = dataclasses.replace(policy, phase_delay_bounds_s=(lo, hi))
    return dataclasses.replace(base, policy=policy)


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


def add_nh_exchange(param, model, labels, start, bounds):
    """One exchange rate k (1/s, log parameter "log_kex") for every listed proton group in every component that
    contains it, all tied to one leader: the N-H / O-H protons exchange with the solvent at one rate
    (physics.exchange; slow exchange model, the 15N isotopologue included). Returns the leader name."""
    from fit_staged import group_index
    leader = None
    for label in labels:
        for c in range(len(model.component_labels)):
            try:
                g = group_index(model, label, c)
            except ValueError:
                continue
            name = param.add_exchange(c, g, start, bounds, name=f"c{c}.log_kex_{label}")
            if leader is None:
                leader = name
            else:
                param.tie(leader, name)
    if leader is None:
        raise ValueError(f"No component contains the exchanging groups {labels}.")
    return leader


def exchange_variants_for(fragment, exchange, overrides):
    from zulf_hypothesis import exchange_variants
    from fit_processed_spectrum import override_couplings
    return exchange_variants(override_couplings(fragment, overrides), exchange)[0]


def exchangeable_labels(fragment):
    """Proton groups on N or O sites (candidates for --nh-exchange)."""
    element = {site.label: site.element for site in fragment.sites}
    return [g.label for g in fragment.protons if element.get(g.site) in ("N", "O")]


class JointSeries:
    """Shared monotonic couplings over several spectra of one structure (see the module docstring)."""

    def __init__(self, model, settings, observations, xs, shared=(), shape="monotone", exchange=None, fixed=()):
        self.model, self.settings = model, settings
        self.params = [settings.parameterize(model.interpretation) for _ in observations]
        for p in self.params:                      # couplings held at their structure values (--free-couplings)
            for key in fixed:
                for n in model.coupling_names.get(key, []):
                    if n in p.parameters:
                        p.parameters[n].free = False
        if exchange:
            for p in self.params:
                add_nh_exchange(p, model, exchange["labels"], exchange["start"], exchange["bounds"])
        self.forwards = [MixtureForward(p, o, SUDDEN_DROP, settings.gain_model, settings.background_order,
                                        settings.band_weighting, **_signal_kwargs(settings))
                         for p, o in zip(self.params, observations)]
        self.free = self.params[0].free_names
        self.coupling = [n for n in self.free if self.params[0].parameters[n].kind == "coupling"]
        # spectrum parameters shared by all spectra (e.g. one delay for the series); the others are per spectrum
        self.shared = [n for n in self.free if n in set(shared) and n not in self.coupling]
        self.local = [n for n in self.free if n not in self.coupling and n not in self.shared]
        self.xs = np.asarray(xs, float)
        if len(self.xs) < 1 or np.any(np.diff(self.xs) < 0):
            raise ValueError("Give the spectra in non-decreasing concentration.")
        # concentration nodes: spectra at the same concentration share one set of couplings
        self.nodes = np.unique(self.xs)
        if len(self.nodes) < 2 and shape == "monotone":
            raise ValueError("A monotone series needs at least two concentrations (one spectrum: shape free).")
        self.node_of = [int(np.searchsorted(self.nodes, x)) for x in self.xs]
        self.nc, self.nl, self.ns = len(self.coupling), len(self.local), len(observations)
        self.nn = len(self.nodes)
        if shape not in ("monotone", "free"):
            raise ValueError(f"shape must be monotone or free, not {shape!r}")
        self.shape = shape
        # monotone: v, A, w_1 .. w_{nodes-1}; free: J at every node (no constraint between concentrations)
        self.m = self.nn + 1 if shape == "monotone" else self.nn
        self.ntheta = self.nc * self.m           # coupling blocks
        self.nt = self.ntheta + len(self.shared)  # then the shared spectrum parameters
        self.col = {n: i for i, n in enumerate(self.free)}
        self.priors = ([], np.zeros(0), np.zeros(0))
        self.smoothing = None                     # per spectrum: residual transform S = W K W^-1 (or None)
        self.peaks = None                         # per spectrum: index windows around data peak tops (or None)
        self.peak_strength = 0.0
        self.res_windows = None                   # per spectrum: residual-vector indices of residual-peak windows
        self.res_strength = 0.0
        self.trace = None                         # list collecting (stage, cost, z) at every residual evaluation
        self.trace_label = ""
        self.monitor = None                       # callable (stage, cost, z) at every evaluation (fit_monitor)

    def set_peak_penalty(self, strength, prominence=0.12, tolerance_hz=0.15, min_sigma=2.0, smooth=0.0, dips=0.0,
                         max_width_hz=None):
        """Missing-peak rows: for every data peak top (prominence >= `prominence` times the largest abs value and
        >= `min_sigma` noise sigma) one extra residual row strength * min(max over +-tolerance_hz of r, 0),
        r = W (model - data) / norm the ordinary residual. The row is zero while the model reaches the data peak
        somewhere near the top and grows when the model is too low there (peak missing or too weak); the
        ordinary residual and its weights are unchanged. Data only (peak tops are fixed). strength 0 = off.

        smooth > 0 replaces max and min(., 0) by smooth versions for the optimiser (the hard rows have kinks,
        which keep least squares from converging): soft max tau log sum exp(r / tau) and soft hinge
        -s log(1 + exp(-m / s)), tau = s = smooth times the peak's own height in residual units. Scores for
        ranking should use the hard rows (smooth 0).

        dips > 0 adds the mirror rows for every data valley (a local minimum of the real part with the same
        prominence rules): dips * max(min over +-tolerance_hz of r, 0), zero while the model comes down to the data
        somewhere near the valley and growing when it fills the valley (one broad line over a doublet) or misses a
        negative line. Together the two kinds penalise structure of the data that the model smooths over; structure
        of the model finer than the data (a broad data line refined into close lines) costs nothing extra.

        max_width_hz: only sharp extrema get rows (width at half prominence below this): the rows then act on the
        sharp lines and narrow valleys of the data (a missed line is a sharp negative residual peak), not on broad
        humps, which the model may refine into close lines without penalty."""
        from scipy.signal import find_peaks
        from zulf_core.solver.forward import signal_regions
        self.peak_strength = float(strength)
        self.peak_smooth = float(smooth)
        self.dip_strength = float(dips)
        if not strength and not dips:
            self.peaks = None
            return
        self.peaks, self.peak_heights, self.peak_signs = [], [], []
        for f in self.forwards:
            y = f.y.real
            _, sigma = signal_regions(f.f, f.y, f.band, 4.0)
            prom = max(prominence * float(np.abs(f.y).max()), min_sigma * sigma)
            spacing = float(np.median(np.diff(f.f)))
            half = max(int(round(tolerance_hz / spacing)), 0)
            wins, heights, signs = [], [], []
            for sign, on in ((1.0, bool(strength)), (-1.0, bool(dips))):
                if not on:
                    continue
                width = (None, max_width_hz / spacing) if max_width_hz else None
                ext, props = find_peaks(sign * y, prominence=prom, width=width)
                wins += [np.arange(max(t - half, 0), min(t + half + 1, len(y))) for t in ext]
                # smoothing scale of a dip row: its prominence (the depth of the valley), like a peak's height
                scale = np.abs(y[ext]) if sign > 0 else props["prominences"]
                heights += list(scale * f.weight[ext] / f.norm)
                signs += [sign] * len(ext)
            self.peaks.append(wins)
            self.peak_heights.append(np.asarray(heights))
            self.peak_signs.append(np.asarray(signs))

    def _peak_rows(self, s, r):
        """Rows of spectrum s from its ordinary residual r, and per row (indices, coefficients) of
        d row / d r (the Jacobian row is the matching combination of ordinary rows)."""
        rows, grads = [], []
        signs = getattr(self, "peak_signs", None)
        for i, (w, h) in enumerate(zip(self.peaks[s], self.peak_heights[s])):
            sg = 1.0 if signs is None else float(signs[s][i])
            st = self.peak_strength if sg > 0 else self.dip_strength
            rw = sg * r[w]                       # a dip row is a peak row of the mirrored residual
            if self.peak_smooth > 0:
                t = self.peak_smooth * h
                a = rw / t
                p = np.exp(a - a.max())
                m = t * (a.max() + np.log(p.sum()))
                p /= p.sum()
                rows.append(-st * t * np.logaddexp(0.0, -m / t))
                grads.append((w, sg * st * expit(-m / t) * p))
            else:
                j = int(np.argmax(rw))
                v = min(float(rw[j]), 0.0)
                rows.append(st * v)
                grads.append((np.array([w[j]]), np.array([sg * st])) if v < 0 else (np.zeros(0, int), np.zeros(0)))
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
        if self.shape == "free":
            return th.copy()
        c, _ = monotone_profile(th[2:])
        return th[0] + th[1] * c

    def coupling_jacobian(self, z, k):
        """d J_k(node) / d theta_k: (nodes, m)."""
        if self.shape == "free":
            return np.eye(self.nn)
        th = self.theta(z, k)
        c, dc = monotone_profile(th[2:])
        return np.column_stack([np.ones(self.nn), c, th[1] * dc])

    # z = [theta_0, ..., theta_{nc-1}, shared, local of spectrum 0 (nl), local of spectrum 1, ...]
    def spectrum_vector(self, z, s, values=None):
        values = values if values is not None else [self.coupling_values(z, k) for k in range(self.nc)]
        x = np.empty(len(self.free))
        for k, n in enumerate(self.coupling):
            x[self.col[n]] = values[k][self.node_of[s]]
        for i, n in enumerate(self.shared):
            x[self.col[n]] = z[self.ntheta + i]
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
            if self.shape == "free":
                z.extend(col)
                continue
            a = col[-1] - col[0]
            steps = np.diff(col) * (1.0 if a >= 0 else -1.0)
            steps = np.maximum(steps, 0.0)
            w = np.log(steps + 1e-3 * max(steps.max(), 1e-3)) if steps.sum() > 0 else np.zeros(self.nn - 1)
            z.extend([col[0], a] + list(np.clip(w - w.max(), -W_BOUND, W_BOUND)))
        X = np.array(per_spectrum_x)
        z.extend(float(np.mean(X[:, self.col[n]])) for n in self.shared)
        z.extend(X[s, self.col[n]] for s in range(self.ns) for n in self.local)
        return np.array(z, float)

    def level_shift(self, draws):
        """Start perturbation from random draws (len(z)): the level of every coupling moves, its shape is kept
        (monotone: the v entries; free: one offset for all nodes of a coupling)."""
        out = np.zeros(len(draws))
        for k in range(self.nc):
            if self.shape == "free":
                out[k * self.m:(k + 1) * self.m] = draws[k * self.m]
            else:
                out[k * self.m] = draws[k * self.m]
        return out

    def bounds(self, change_bound):
        lo, hi = self.params[0].bounds()
        lz, hz = [], []
        for n in self.coupling:
            if self.shape == "free":
                lz += [lo[self.col[n]]] * self.nn
                hz += [hi[self.col[n]]] * self.nn
                continue
            lz += [lo[self.col[n]], -change_bound] + [-W_BOUND] * (self.nn - 1)
            hz += [hi[self.col[n]], change_bound] + [W_BOUND] * (self.nn - 1)
        lz.extend(lo[self.col[n]] for n in self.shared)
        hz.extend(hi[self.col[n]] for n in self.shared)
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
            if self.res_windows is not None and len(self.res_windows[s]):
                parts.append(self.res_strength * r[self.res_windows[s]])
        ks, mean, scale = self.priors
        if ks:
            parts.append((np.array([values[k].mean() for k in ks]) - mean) * scale)
        out = np.concatenate(parts)
        if self.trace is not None:
            self.trace.append((self.trace_label, float(out @ out), np.array(z, float)))
        if self.monitor is not None:                 # observes only: the returned residual is unchanged
            self.monitor(self.trace_label, float(out @ out), z)
        return out

    def jacobian(self, z):
        values = [self.coupling_values(z, k) for k in range(self.nc)]
        dvals = [self.coupling_jacobian(z, k) for k in range(self.nc)]
        blocks = []
        for s, f in enumerate(self.forwards):
            js = f.jacobian(self.spectrum_vector(z, s, values))
            out = np.zeros((js.shape[0], len(z)))
            for k, n in enumerate(self.coupling):
                out[:, k * self.m:(k + 1) * self.m] = np.outer(js[:, self.col[n]], dvals[k][self.node_of[s]])
            for i, n in enumerate(self.shared):
                out[:, self.ntheta + i] = js[:, self.col[n]]
            for i, n in enumerate(self.local):
                out[:, self.nt + s * self.nl + i] = js[:, self.col[n]]
            blocks.append(self._smooth(s, out))
            if self.peaks is not None:
                _, grads = self._peak_rows(s, f.predict(self.spectrum_vector(z, s, values)).residual)
                blocks.append(np.array([c @ out[idx] for idx, c in grads]).reshape(len(grads), len(z)))
            if self.res_windows is not None and len(self.res_windows[s]):
                blocks.append(self.res_strength * out[self.res_windows[s]])
        ks, _, scale = self.priors
        if ks:
            rows = np.zeros((len(ks), len(z)))
            for r, k in enumerate(ks):
                rows[r, k * self.m:(k + 1) * self.m] = dvals[k].mean(axis=0) * scale[r]
            blocks.append(rows)
        return np.vstack(blocks)

    def find_residual_peaks(self, z, k_sigma=4.0, half_width_hz=0.3, assign_hz=0.6, noise_window_hz=5.0,
                            max_width_hz=1.5, min_relative_amplitude=0.02):
        """Localized residual peaks of either sign: |model - data| (complex magnitude of the weighted residual, so a
        model line too tall or a filled data dip counts as much as a missing line) standing out by
        k_sigma (prominence) above its surroundings, in units of a robust local noise level (1.4826 MAD of the real and imaginary residual within +-noise_window_hz), narrower
        than max_width_hz at half height. Broad, spread-out residual is not a peak. A peak is assignable when a
        model transition (|amplitude| >= min_relative_amplitude of its component's largest) lies within assign_hz:
        then a coupling can move a line onto it; otherwise (an unexplained feature: impurity, instrument line,
        missing species) it is reported only. Returns per spectrum the residual-vector indices of the windows
        (+-half_width_hz) around assignable peaks, and the report."""
        from scipy.signal import find_peaks
        windows, report = [], []
        for s, f in enumerate(self.forwards):
            x = self.spectrum_vector(z, s)
            r = f.predict(x).residual
            nf = len(f.f)
            rc = r[:nf] + 1j * r[nf:] if len(r) == 2 * nf else r.astype(complex)
            step = float(np.median(np.diff(f.f)))
            half = max(int(noise_window_hz / step), 5)
            centres = np.arange(0, nf, max(half // 10, 1))
            parts = np.c_[rc.real, rc.imag]
            sig = []
            for c in centres:
                block = parts[max(c - half, 0):c + half + 1]
                sig.append(1.4826 * np.median(np.abs(block - np.median(block, axis=0))))
            sigma = np.maximum(np.interp(np.arange(nf), centres, sig), 1e-30)
            snr = np.abs(rc) / sigma
            # prominence (height above the surrounding residual), not height: ripples on a broad misfit are not peaks,
            # and the broad misfit itself is wider than max_width_hz at half its prominence
            tops, props = find_peaks(snr, prominence=k_sigma, width=(None, max_width_hz / step), rel_height=0.5)
            values = self.params[s].values(x)
            lines = []
            for c, system in enumerate(self.params[s].systems(values)):
                tl = f.transitions(values, c, system)
                a = np.abs(np.asarray(tl.amplitudes))
                if len(a):
                    lines.extend(np.asarray(tl.frequencies_hz)[a >= min_relative_amplitude * a.max()])
            lines = np.asarray(lines)
            idx = []
            for t, h in zip(tops, props["prominences"]):
                fp = float(f.f[t])
                near = lines[np.abs(lines - fp) <= assign_hz] if len(lines) else lines
                # sign of the real (absorption) residual at the top: model above the data (+, e.g. a dip the model
                # fills, a line too tall) or below it (-, e.g. a data line the model misses); both are peaks of |r|
                report.append({"spectrum": s, "frequency_hz": round(fp, 3), "height_sigma": round(float(h), 1),
                               "sign": "model above data" if rc[t].real > 0 else "model below data",
                               "assignable": bool(len(near)), "nearest_transition_hz":
                               round(float(near[np.argmin(np.abs(near - fp))]), 3) if len(near) else None})
                if len(near):
                    idx.append(np.flatnonzero(np.abs(f.f - fp) <= half_width_hz))
            idx = np.unique(np.concatenate(idx)) if idx else np.zeros(0, int)
            windows.append(np.r_[idx, idx + nf] if len(r) == 2 * nf else idx)
        return windows, report

    def peak_sources(self, z, s, frequency_hz, assign_hz=0.6, bands=None, min_relative_amplitude=0.02):
        """Where a feature of spectrum s comes from, from the analytic transitions: every transition within assign_hz
        of frequency_hz (component, frequency, amplitude relative to the component's largest, df/dJ for every free
        coupling, tied names summed into their leader), and per coupling its local effect (sum over those lines of
        |A| |df/dJ|) and its selectivity (local effect / the same sum over all lines in `bands`, default the
        spectrum's fit bands). A coupling with a large local effect and a high selectivity moves this feature and
        little else: the one to adjust first. `shift` / `split`: weighted mean / spread of df/dJ over the nearby lines
        (a split-type coupling changes the spacing inside a feature, e.g. a dip between close lines)."""
        from zulf_core.physics.derivatives import transition_derivatives
        f = self.forwards[s]
        P = self.params[s]
        x = self.spectrum_vector(z, s)
        values = P.values(x)
        bands = bands or [(float(f.f.min()), float(f.f.max()))]
        lines, local, total = [], {n: 0.0 for n in self.coupling}, {n: 0.0 for n in self.coupling}
        near_cols = {n: [] for n in self.coupling}
        near_amp = []
        for c, system in enumerate(P.systems(values)):
            names = P.coupling_names(c)
            pairs = [P.parameters[n].detail for n in names]
            d = transition_derivatives(system, pairs, f.protocol)
            if not len(d):
                continue
            amp = np.abs(d.amplitudes)
            rel = amp / amp.max()
            dfdj = {}
            for j, n in enumerate(names):
                leader = P.ties.get(n, n)
                if leader not in self.col or leader not in local:
                    continue
                with np.errstate(divide="ignore", invalid="ignore"):
                    col = np.where(amp > 0, (d.t_weights[j] / d.amplitudes).real, 0.0)
                dfdj[leader] = dfdj.get(leader, 0.0) + col
            inband = np.zeros(len(d), bool)
            for lo, hi in bands:
                inband |= (d.frequencies_hz >= lo) & (d.frequencies_hz <= hi)
            near = (np.abs(d.frequencies_hz - frequency_hz) <= assign_hz) & (rel >= min_relative_amplitude)
            near_amp.extend((amp[near] / amp.max()).tolist())
            for n in self.coupling:
                near_cols[n].extend((dfdj[n][near] if n in dfdj else np.zeros(int(near.sum()))).tolist())
            for n, col in dfdj.items():
                local[n] += float(np.sum(amp[near] * np.abs(col[near])))
                total[n] += float(np.sum(amp[inband] * np.abs(col[inband])))
            for i in np.flatnonzero(near):
                lines.append({"component": c, "frequency_hz": float(d.frequencies_hz[i]), "relative_amplitude":
                              float(rel[i]), "df_dJ": {n: float(col[i]) for n, col in dfdj.items()}})
        w = np.asarray(near_amp)
        effect = {}
        for n in self.coupling:
            col = np.asarray(near_cols[n])
            # shift: amplitude-weighted mean df/dJ of the nearby lines (moves the feature as a whole);
            # split: their weighted spread (changes the spacing inside the feature, e.g. opens or fills a dip)
            shift = float(np.sum(w * col) / w.sum()) if len(w) else 0.0
            split = float(np.sqrt(np.sum(w * (col - shift) ** 2) / w.sum())) if len(w) > 1 else 0.0
            effect[n] = {"local": local[n], "selectivity": local[n] / total[n] if total[n] > 0 else 0.0,
                         "shift": shift, "split": split}
        return {"frequency_hz": float(frequency_hz), "lines": lines, "couplings": effect}

    def local_fit(self, z, s, window_hz, couplings, locals_=(), lower=None, upper=None, starts=None, max_nfev=100,
                  margin_hz=2.0):
        """Targeted local refinement: fit only spectrum s inside window_hz = (lo, hi), with only the named couplings
        (leaders, e.g. the split-type couplings of peak_sources) and spectrum parameters (e.g. the rate family holding
        those lines) free; everything else held at z. Fast: a forward model of the window points alone that renders
        and differentiates only the transitions within window +- margin_hz (farther lines add only their tails;
        gains, phase and background are solved again on the window). `starts`: list of dicts {name: start value}
        (default: z itself). Returns [(local cost, z)] sorted, each z a full vector for a global refit."""
        f = self.forwards[s]
        lo_w, hi_w = float(window_hz[0]), float(window_hz[1])
        obs_w = f.obs.restricted([(lo_w, hi_w)])
        fw = MixtureForward(self.params[s], obs_w, f.protocol, self.settings.gain_model,
                            self.settings.background_order, self.settings.band_weighting,
                            **_signal_kwargs(self.settings))
        fw.line_band_hz = (lo_w - margin_hz, hi_w + margin_hz)
        fw.jacobian_only = set(couplings) | set(locals_)
        zi, xcols = [], []
        for n in couplings:
            k = self.coupling.index(n)
            zi.append(k * self.m + (self.node_of[s] if self.shape == "free" else 0))
            xcols.append(self.col[n])
        for n in locals_:
            if n in self.shared:
                zi.append(self.ntheta + self.shared.index(n))
            else:
                zi.append(self.nt + s * self.nl + self.local.index(n))
            xcols.append(self.col[n])
        zi = np.asarray(zi, int)
        lo = lower[zi] if lower is not None else np.full(len(zi), -np.inf)
        hi = upper[zi] if upper is not None else np.full(len(zi), np.inf)

        def full(p):
            zz = z.copy()
            zz[zi] = p
            return zz

        def res(p):
            return fw.predict(self.spectrum_vector(full(p), s)).residual

        def jac(p):
            return fw.jacobian(self.spectrum_vector(full(p), s))[:, xcols]
        out = []
        names = list(couplings) + list(locals_)
        for st in (starts or [{}]):
            p0 = z[zi].copy()
            for n, v in st.items():
                p0[names.index(n)] = v
            p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)
            sol = least_squares(res, p0, jac=jac, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev,
                                ftol=1e-10, xtol=1e-10, gtol=1e-10)
            out.append((float(2 * sol.cost), full(sol.x)))
        out.sort(key=lambda t: t[0])
        return out

    def data_residuals(self, z):
        out = []
        for s, f in enumerate(self.forwards):
            pred = f.predict(self.spectrum_vector(z, s))
            m = f.data_signal_mask if f.data_signal_mask is not None else np.ones(len(f.y), bool)
            out.append(float(np.linalg.norm(f.mismatch(pred.model, m)) /
                             max(np.linalg.norm(f.mismatch(np.zeros_like(f.y), m)), 1e-30)))
        return out


def predict_left_out(joint, z, entry, band, key_of, settings, max_nfev=200, real_only=True):
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
    obs = ObservedSpectrum.from_spectrum(np.load(entry["freq"]).astype(float), _values(entry, real_only),
                                         [band], record=entry.get("record"), phasing=entry.get("phasing"), real_only=real_only,
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
        for m, n in enumerate(joint.shared):           # held at the series value
            xv[col[n]] = z[joint.ntheta + m]
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


def _write_trace(out, joint, series, trace, z_final, frames, key_of, origin):
    """Fit trace for a viewer (scripts/trace_view.py): up to `frames` evaluations of the recorded path (evenly
    spread, first and last included, plus the final solution), each with the model spectra, the couplings and the
    objective. Saved as OUT/trace.npz (arrays) and OUT/trace.json (labels, couplings, objective)."""
    if not trace:
        return
    out.mkdir(parents=True, exist_ok=True)
    picks = sorted(set(np.linspace(0, len(trace) - 1, min(frames, len(trace))).round().astype(int).tolist()))
    entries = [trace[i] for i in picks] + [("final", float(np.sum(joint.residual(z_final) ** 2)), z_final)]
    data, meta = {}, {"origin": origin, "evaluations": len(trace), "frames": []}
    for s, (e, f) in enumerate(zip(series, joint.forwards)):
        data[f"f{s}"] = f.f
        data[f"y{s}"] = f.y
    models = [[] for _ in joint.forwards]
    for i, (label, cost, zz) in zip(picks + [len(trace)], entries):
        values = [joint.coupling_values(zz, k) for k in range(joint.nc)]
        for s, f in enumerate(joint.forwards):
            models[s].append(f.predict(joint.spectrum_vector(zz, s, values)).model)
        meta["frames"].append({"evaluation": int(i), "stage": label, "objective": cost,
                               "J": {key_of.get(n, n): values[k].tolist() for k, n in enumerate(joint.coupling)}})
    for s, m in enumerate(models):
        data[f"model{s}"] = np.array(m)
    meta["spectra"] = [{"id": e["id"], "x": e["x"]} for e in series]
    np.savez_compressed(out / "trace.npz", **data)
    (out / "trace.json").write_text(json.dumps(meta, indent=1))


def _residual_peak_stage(joint, solutions, lower, upper, args):
    """Outer loop on the best candidates: find localized, assignable residual peaks, refit with those windows
    weighted (rows res_strength * r there; windows fixed within a refit), repeat. A broad residual hardly affects
    the couplings, a localized one marks a structure the model does not match. All candidates, before and after,
    are then scored on one common window set (the union of every detection), so the ranking does not depend on
    which windows a candidate happened to see; the plain objective is kept for the record."""
    opts = dict(k_sigma=args.residual_peak_sigma, half_width_hz=args.residual_peak_halfwidth,
                assign_hz=args.residual_peak_assign)
    joint.res_strength = float(args.residual_peaks)
    hard_smooth = joint.peak_smooth if joint.peaks is not None else 0.0
    candidates = [zz for _, zz in solutions[:max(args.residual_peak_candidates, 1)]]
    union = [set() for _ in joint.forwards]
    rounds, refined = [], []
    for ci, zz in enumerate(candidates):
        for rnd in range(args.residual_peak_rounds):
            windows, report = joint.find_residual_peaks(zz, **opts)
            for s, w in enumerate(windows):
                union[s].update(w.tolist())
            rounds.append({"candidate": ci, "round": rnd, "peaks": report})
            if not any(len(w) for w in windows):
                break
            joint.res_windows = windows
            joint.peak_smooth = args.peak_smooth if joint.peaks is not None else 0.0
            zz = least_squares(joint.residual, zz, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                               max_nfev=args.max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10).x
            print(f"residual peaks: candidate {ci} round {rnd}: {sum(p['assignable'] for p in report)} assignable "
                  f"of {len(report)}", flush=True)
        windows, report = joint.find_residual_peaks(zz, **opts)
        for s, w in enumerate(windows):
            union[s].update(w.tolist())
        rounds.append({"candidate": ci, "round": "final", "peaks": report})
        refined.append(zz)
    joint.peak_smooth = 0.0 if joint.peaks is not None else hard_smooth
    common = [np.array(sorted(u), int) for u in union]
    scored = []
    for kind, pool in (("before", candidates), ("after", refined)):
        for ci, zz in enumerate(pool):
            joint.res_windows = common
            penalized = float(np.sum(joint.residual(zz) ** 2))
            joint.res_windows = None
            plain = float(np.sum(joint.residual(zz) ** 2))
            scored.append({"candidate": ci, "stage": kind, "penalized": penalized, "plain": plain, "z": zz})
    scored.sort(key=lambda d: d["penalized"])
    best = scored[0]
    rest = [(sc, zz) for sc, zz in solutions if not np.array_equal(zz, best["z"])]
    record = {"strength": args.residual_peaks, "sigma": args.residual_peak_sigma,
              "half_width_hz": args.residual_peak_halfwidth, "assign_hz": args.residual_peak_assign,
              "common_window_points": [int(len(c)) for c in common], "rounds": rounds,
              "scores": [{k: v for k, v in d.items() if k != "z"} for d in scored],
              "chosen": {"candidate": best["candidate"], "stage": best["stage"]}}
    return [(best["plain"], best["z"])] + rest, record


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
            blk = slice(i, i + joint.m) if joint.shape == "free" else slice(i, i + 1)   # free: shift all nodes
            for v in grid:
                zz = z.copy()
                zz[blk] += v - z[i]
                r = joint.residual(zz)
                c = float(r @ r)
                if best_c is None or c < best_c:
                    best_v, best_c = v, c
            z = z.copy()
            z[blk] += best_v - z[i]
        z = least_squares(joint.residual, z, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                          max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10).x
    return z


_MONITOR = None        # (directory, coupling keys) when fit_monitor records the starts



# ---------------------------------------------------------------------------------------------------------------
# Component search: every isotopologue on the window where it dominates, its own couplings from random starts
# ---------------------------------------------------------------------------------------------------------------

_CS_TASK = None


def _component_windows(joint, z, s, min_line=0.05, gap_hz=3.0, dominance=0.8):
    """[(component, lo, hi)] for spectrum s: the lines of each component (relative amplitude >= min_line) clustered
    at gaps > gap_hz, each cluster +- 1 Hz, merged with the next cluster of the same component whenever the merged
    window stays dominated by that component (>= dominance of the summed component magnitudes on its points), so
    a component alone in a region gets one window over all of it."""
    f = joint.forwards[s]
    P = joint.params[s]
    values = P.values(joint.spectrum_vector(z, s))
    pred = f.predict(joint.spectrum_vector(z, s))
    nf = len(f.f)
    comps = [np.abs(np.asarray(c)[:nf] if len(c) == nf else np.asarray(c)) for c in pred.component_spectra]
    f_lo, f_hi = float(f.f.min()), float(f.f.max())

    def share(c, lo, hi):
        w = (f.f >= lo) & (f.f <= hi)
        if not w.any():
            return 0.0
        tot = sum(float(np.linalg.norm(x[w])) for x in comps)
        return float(np.linalg.norm(comps[c][w])) / tot if tot > 0 else 0.0

    out = []
    for c, system in enumerate(P.systems(values)):
        tl = f.transitions(values, c, system)
        a = np.abs(np.asarray(tl.amplitudes))
        fr = np.asarray(tl.frequencies_hz)
        if not len(a):
            continue
        fr = np.sort(fr[(a >= min_line * a.max()) & (fr >= f_lo) & (fr <= f_hi)])
        if not len(fr):
            continue
        clusters, start = [], fr[0]
        for x0, x1 in zip(fr[:-1], fr[1:]):
            if x1 - x0 > gap_hz:
                clusters.append([start - 1.0, x0 + 1.0])
                start = x1
        clusters.append([start - 1.0, fr[-1] + 1.0])
        merged = []
        for cl in clusters:
            if merged and share(c, merged[-1][0], cl[1]) >= dominance:
                merged[-1][1] = cl[1]
            else:
                merged.append(cl)
        out += [(c, max(lo, f_lo), min(hi, f_hi)) for lo, hi in merged if share(c, lo, hi) >= dominance]
    return out


def _cs_local(start):
    joint, s, window, names, rates, lower, upper, max_nfev, z = _CS_TASK
    return joint.local_fit(z, s, window, names, locals_=rates, lower=lower, upper=upper, starts=[start],
                           max_nfev=max_nfev)[0]


def component_search(joint, prob, z, args, stage):
    """For every spectrum and every component, on each window it dominates (_component_windows): its own couplings
    (J(site, ...) of its labelled nucleus; shared H-H couplings held) and its rate families in the window are refit
    on the window alone (JointSeries.local_fit, the window's own gains) from the current values and from
    --component-search-starts random starts (one-bond couplings +- --component-search-box[2] Hz, the others uniform
    in the box). The best distinct window solutions are scored on the whole objective with everything else held;
    the best improving one is accepted (greedy, worst bands first). Returns (z, record)."""
    global _CS_TASK
    lo_box, hi_box, one_bond = (float(v) for v in args.component_search_box.split(","))
    rng = np.random.default_rng(args.seed + (0 if stage == "start" else 1))
    key_of = prob.key_of
    labels = joint.model.component_labels
    edges = np.asarray(joint.params[0].policy.family_edges_hz or [], float)
    z = z.copy()
    whole = float(np.sum(joint.residual(z) ** 2))
    record = {"stage": stage, "objective_before": whole, "windows": []}
    t0 = time.time()
    for s in range(len(joint.forwards)):
        windows = _component_windows(joint, z, s, args.component_search_min_line, args.component_search_gap,
                                     args.component_search_dominance)
        for c, lo, hi in windows:
            site = labels[c].split("@")[-1].split(" ")[0] if "@" in labels[c] else None
            names = [n for n in joint.coupling if site and key_of.get(n, n).startswith(f"J({site},")]
            if not names:
                continue
            fams = set(np.searchsorted(edges, np.linspace(lo, hi, 200), side="right").tolist()) if len(edges) \
                else {0}
            rates = [n for n in joint.local if n.startswith(f"c{c}.log_rate") and int(n.split("log_rate")[1]) in fams]
            current = {n: float(joint.coupling_values(z, joint.coupling.index(n))[joint.node_of[s]]) for n in names}
            starts = [{}]
            for _ in range(args.component_search_starts):
                starts.append({n: float(rng.uniform(v - one_bond, v + one_bond) if abs(v) >= 50.0
                                        else rng.uniform(lo_box, hi_box)) for n, v in current.items()})
            _CS_TASK = (joint, s, (lo, hi), names, rates, prob.lower, prob.upper, args.component_search_nfev, z)
            if args.workers > 1:
                import multiprocessing
                with multiprocessing.get_context("fork").Pool(args.workers) as pool:
                    results = pool.map(_cs_local, starts, chunksize=1)
            else:
                results = [_cs_local(st) for st in starts]
            here = results[0][0]
            results.sort(key=lambda t: t[0])
            distinct = []
            for cost, zz in results:
                vals = {n: float(joint.coupling_values(zz, joint.coupling.index(n))[joint.node_of[s]]) for n in names}
                if all(max(abs(vals[n] - d[2][n]) for n in names) > 0.15 for d in distinct):
                    distinct.append((cost, zz, vals))
                if len(distinct) == args.component_search_keep:
                    break
            scored = [(float(np.sum(joint.residual(zz) ** 2)), cost, zz, vals) for cost, zz, vals in distinct]
            scored.sort(key=lambda t: t[0])
            best = scored[0]
            accepted = best[0] < whole * (1 - 1e-4)
            entry = {"spectrum": s, "component": labels[c], "window_hz": [lo, hi],
                     "couplings": [key_of.get(n, n) for n in names], "rates": rates,
                     "window_cost_here": here, "window_cost_best": best[1], "objective_before": whole,
                     "objective_best": best[0], "accepted": bool(accepted),
                     "values_before": {key_of.get(n, n): v for n, v in current.items()},
                     "values_best": {key_of.get(n, n): v for n, v in best[3].items()}}
            record["windows"].append(entry)
            print(f"component search ({stage}): {labels[c]} {lo:.1f}-{hi:.1f} Hz, {len(starts)} starts: window cost "
                  f"{here:.4f} -> {best[1]:.4f}; whole objective {whole:.5f} -> {best[0]:.5f}"
                  f"{' accepted' if accepted else ''}  ({time.time() - t0:.0f} s)", flush=True)
            if accepted:
                z, whole = best[2].copy(), best[0]
    record["objective_after"] = whole
    record["seconds"] = round(time.time() - t0)
    return z, record


def _solve_indexed(item):
    """_solve_start for start k, recorded by a fit_monitor.MonitorWriter when monitoring is on (observation only)."""
    k, z = item
    if _MONITOR is None:
        return _solve_start(z)
    from fit_monitor import MonitorWriter
    joint = _JOINT_TASK[0]
    joint.monitor = MonitorWriter(_MONITOR[0], f"start_{k:03d}", joint, _MONITOR[1])
    score = None
    try:
        out = _solve_start(z)
        score = out[0]
        return out
    finally:
        joint.monitor.close(score)
        joint.monitor = None


def _solve_start(z):
    """One start (module level so worker processes can run it; fork start method): least squares at every level
    of the smoothing schedule in turn (coarse to fine; the last level is always unsmoothed); returns the score of
    the unsmoothed objective."""
    joint, lower, upper, max_nfev, schedule, scan, hold_small = _JOINT_TASK
    tracing = joint.trace is not None
    if tracing:
        joint.trace = []
    joint.trace_label = "hold small couplings" if hold_small else "start"      # stage name for trace and monitor
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
        joint.trace_label = "coordinate scan"
        z = _coordinate_scan(joint, z, lower, upper, max_nfev, **scan)
    for sigma in list(schedule) + [0.0]:
        joint.set_smoothing(sigma)
        joint.trace_label = f"smoothing {sigma:g} Hz" if sigma else "fit"
        sol = least_squares(joint.residual, z, jac=joint.jacobian, bounds=(lower, upper), x_scale="jac",
                            max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        z = sol.x
    if tracing:
        return float(2 * sol.cost), sol.x, joint.trace
    return float(2 * sol.cost), sol.x


def build_problem(args):
    """The fitting problem of a command line: observations, model, JointSeries (weights, priors, missing-peak
    rows), the start vector z0 and the bounds. Shared by main and scripts/j_tuner.py."""
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
        return ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), _values(e, args.real_only),
                                              [tuple(r) for r in e.get("ranges", [(lo, hi)])], record=e.get("record"),
                                              phasing=e.get("phasing"), real_only=args.real_only, label=e["id"])

    obs = [observation(e) for e in series]
    spec = json.loads(args.structure)
    if isinstance(spec, list):
        model = _combined_model(spec, args, lo, hi)
        fragment, fixed_keys = model.fragment, []
        return _finish_problem(args, series, left_out, lo, hi, obs, model, fragment, fixed_keys)
    spec.setdefault("compound", "series")
    overrides = json.loads(args.couplings)
    fixed_keys = []
    if getattr(args, "free_couplings", ""):
        # every coupling not matching the pattern is held, at its value in --from-joint (else the structure's)
        import re
        pattern = re.compile(args.free_couplings)
        previous = json.load(open(args.from_joint))["couplings"] if args.from_joint else {}
        probe = exchange_variants_for(reg.structure_for(spec), args.exchange, overrides)
        for key in build_model(probe, ranges=[(lo, hi)]).coupling_names:
            if not pattern.search(key):
                fixed_keys.append(key)
                if key in previous:
                    overrides[key] = float(previous[key]["J_at_x"][0])
    fragment = override_couplings(reg.structure_for(spec), overrides)
    from zulf_hypothesis import exchange_variants
    fragment = exchange_variants(fragment, args.exchange)[0]     # fast: protons on N/O/S dropped (decoupled)
    model = build_model(fragment, ranges=[(lo, hi)])
    return _finish_problem(args, series, left_out, lo, hi, obs, model, fragment, fixed_keys)


def _combined_model(specs, args, lo, hi):
    """Several molecules (a JSON list of structure specs) as one model: combine_models keeps the ratios fixed only
    within each molecule (free between them) and prefixes the coupling keys P1:, P2:, ... (--couplings takes the
    prefixed keys)."""
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    from zulf_hypothesis import exchange_variants
    from zulf_hypothesis.builder import combine_models
    overrides = json.loads(args.couplings)
    models = []
    for k, sp in enumerate(specs):
        sp = dict(sp)
        sp.setdefault("compound", f"part {k + 1}")
        own = {key.split(":", 1)[1]: v for key, v in overrides.items() if key.startswith(f"P{k + 1}:")}
        frag = exchange_variants(override_couplings(reg.structure_for(sp), own), args.exchange)[0]
        models.append(build_model(frag, ranges=[(lo, hi)]))
    return combine_models(models)


def _finish_problem(args, series, left_out, lo, hi, obs, model, fragment, fixed_keys):
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    if getattr(args, "signal_height_power", 0.0):
        base = dataclasses.replace(base, signal_height_power=args.signal_height_power)
    base = _rate_policy(base, args)
    settings = _settings_for(model, args.variant, base)
    exchange = None
    if getattr(args, "nh_exchange", 0.0) > 0:
        if args.exchange != "slow":
            raise SystemExit("--nh-exchange needs --exchange slow (the N-H protons kept in the spin system).")
        lo_k, hi_k = (float(v) for v in args.nh_exchange_bounds.split(","))
        labels = ([v.strip() for v in args.nh_exchange_labels.split(",") if v.strip()] if args.nh_exchange_labels
                  else exchangeable_labels(fragment))
        exchange = {"labels": labels, "start": args.nh_exchange, "bounds": (lo_k, hi_k)}
    joint = JointSeries(model, settings, obs, [e["x"] for e in series],
                        shared=[n.strip() for n in args.shared.split(",") if n.strip()], shape=args.shape,
                        exchange=exchange, fixed=fixed_keys)
    missing = [n for n in args.shared.split(",") if n.strip() and n.strip() not in joint.shared]
    if missing:
        raise SystemExit(f"--shared: not a spectrum parameter of this model: {missing} (have {joint.local})")
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
        # couplings the previous fit does not have (a model extended since) keep their structure values
        table = np.array([[previous["couplings"][key_of[n]]["J_at_x"][s] if key_of[n] in previous["couplings"]
                           else table[i, k] for k, n in enumerate(joint.coupling)] for i, s in enumerate(rows)])
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
                           smooth=args.peak_smooth, dips=args.dip_penalty, max_width_hz=args.peak_max_width or None)
    lower, upper = joint.bounds(args.change_bound)
    z0 = np.clip(z0, lower + 1e-9, upper - 1e-9)
    return SimpleNamespace(series=series, left_out=left_out, lo=lo, hi=hi, obs=obs, model=model,
                           settings=settings, joint=joint, key_of=key_of, centre=centre, xs0=xs0, z0=z0,
                           prior=(ks, mean, sigma), lower=lower, upper=upper)


def make_parser():
    """Command line of the series fitter (also used by scripts/j_tuner.py)."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--real-only", type=_flag, default=True,
                    help="true (default): fit the real part of the phased spectra; false: fit complex spectra "
                         "(values stored complex, model through the same record and phasing)")
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
    ap.add_argument("--peak-max-width", type=float, default=0.0,
                    help="peak tops and dips: only extrema narrower than this (Hz, width at half prominence; 0 = any)")
    ap.add_argument("--dip-penalty", type=float, default=0.0,
                    help="filled-dip rows: strength per data valley (or negative line) the model stays above (0 = off); "
                         "with --peak-penalty it stops one broad line from covering a resolved splitting")
    ap.add_argument("--peak-min-sigma", type=float, default=2.0, help="peak tops: prominence at least this many noise sigma")
    ap.add_argument("--peak-smooth", type=float, default=0.0,
                    help="optimise with smooth missing-peak rows (width = this fraction of each peak's height); "
                         "solutions are then rescored and ranked with the hard rows")
    ap.add_argument("--family-edges", default="", help="comma-separated transition frequencies (Hz) splitting "
                    "every isotopologue's lines into decay-rate families (default: one rate per isotopologue)")
    ap.add_argument("--phase-delay-bounds", default="", help="lo,hi of the fitted delay in ms (instrument prior)")
    ap.add_argument("--rate-bounds", default="", help="lo,hi decay-rate bounds in 1/s (default: the policy's)")
    ap.add_argument("--exchange", default="slow", choices=["slow", "fast"],
                    help="N-H / O-H protons: slow (default, kept in the spin system) or fast (dropped: decoupled)")
    ap.add_argument("--nh-exchange", type=float, default=0.0,
                    help="start value (1/s) of one fitted exchange rate k of the N-H / O-H protons with the solvent "
                         "(exact Liouville model, physics.exchange; needs --exchange slow; 0 = static, off)")
    ap.add_argument("--nh-exchange-bounds", default="0.05,5000", help="bounds of k (1/s)")
    ap.add_argument("--component-search", default="both", choices=["off", "start", "end", "both"],
                    help="component search: every isotopologue on the window it dominates, its own couplings from "
                         "random starts (window fit, shared H-H held); 'start' guides the multi-start (its result is "
                         "the first centre), 'end' searches from the best solution and refits globally if it wins")
    ap.add_argument("--component-search-starts", type=int, default=30, help="random starts per window")
    ap.add_argument("--component-search-nfev", type=int, default=30, help="max_nfev of each window fit")
    ap.add_argument("--component-search-keep", type=int, default=6, help="distinct window solutions scored")
    ap.add_argument("--component-search-box", default="-8,13,7",
                    help="small couplings uniform in lo,hi; one-bond couplings current +- the third value (Hz)")
    ap.add_argument("--component-search-dominance", type=float, default=0.8)
    ap.add_argument("--component-search-min-line", type=float, default=0.05)
    ap.add_argument("--component-search-gap", type=float, default=3.0)
    ap.add_argument("--free-couplings", default="",
                    help="regex: only couplings whose key matches are fitted; the others are held at their "
                         "--from-joint values (e.g. 'HN|N1' for a staged N-H fit)")
    ap.add_argument("--nh-exchange-labels", default="", help="exchanging proton groups (default: all on N / O)")
    ap.add_argument("--shape", default="monotone", choices=["monotone", "free"],
                    help="monotone (default): every coupling monotonic in x; free: independent J at every x")
    ap.add_argument("--shared", default="", help="comma-separated spectrum parameters shared by all spectra "
                    "(e.g. phase_delay: one model delay for the series)")
    ap.add_argument("--residual-peaks", type=float, default=0.0,
                    help="strength of the residual-peak rows (0 = off): after the starts, the best candidates are "
                         "refitted with localized, assignable residual peaks weighted, then ranked on one common "
                         "window set")
    ap.add_argument("--residual-peak-sigma", type=float, default=4.0, help="residual peaks: height in local noise sigma")
    ap.add_argument("--residual-peak-halfwidth", type=float, default=0.3, help="residual peaks: window half width (Hz)")
    ap.add_argument("--residual-peak-assign", type=float, default=0.6,
                    help="residual peaks: a model transition within this distance (Hz) makes a peak assignable")
    ap.add_argument("--residual-peak-rounds", type=int, default=3)
    ap.add_argument("--residual-peak-candidates", type=int, default=3)
    ap.add_argument("--trace", type=int, default=0,
                    help="record the fit path of the best start (and the residual-peak stage) as up to this many "
                         "frames: OUT/trace.npz + trace.json for scripts/trace_view.py (0 = off)")
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
    ap.add_argument("--model-line-passes", type=int, default=0,
                    help="after the fit: refits with the best model's own lines (envelope above "
                         "--model-line-threshold noise sigma) added to the signal-weighting cores, so model lines "
                         "where the data show none are fully weighted (0 = off)")
    ap.add_argument("--model-line-threshold", type=float, default=2.0)
    ap.add_argument("--max-nfev", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/processed/joint")
    ap.add_argument("--monitor", default="on", choices=["on", "off"],
                    help="record the objective of every start live in OUT/monitor (scripts/fit_monitor.py shows it)")
    return ap


def model_line_passes(joint, z, solutions, passes, threshold, index):
    """Model-line passes of signal weighting (as refine's signal_model_passes): the grid points where the best
    model's own line envelope exceeds `threshold` noise sigma join the peak cores of every spectrum, so a line the
    model puts where the data show none (or a dip) is weighted fully, and the best solution is refit; repeated
    until no point is added. The objective changes with the weights, so every solution is rescored under the final
    weights; the record keeps the objective of the final vector under the original weights for comparison with
    runs without passes."""
    saved = [(f.signal_cores.copy(), f.weight.copy(), f.norm, f.signal_mask) for f in joint.forwards]
    record = {"threshold_sigma": threshold, "passes": []}
    for k in range(passes):
        added = []
        for s, f in enumerate(joint.forwards):
            x = joint.spectrum_vector(z, s)
            pred = f.predict(x)
            points = f.model_line_points(f.p.values(x), np.asarray(pred.gains), threshold)
            new = points & ~f.signal_cores
            added.append({"points": int(new.sum()),
                          "hz": [round(float(v), 2) for v in f.f[new][:: max(1, int(new.sum()) // 40)]]})
            f.add_signal_cores(new)
        if not any(a["points"] for a in added):
            break
        before = float(np.sum(joint.residual(z) ** 2))
        refit = _solve_indexed((index + k, z))
        after = float(np.sum(joint.residual(refit[1]) ** 2))
        if after < before:
            z = refit[1]
        record["passes"].append({"added": added, "objective_before_refit": before, "objective_after_refit": after})
        print(f"model-line pass {k + 1}: +{sum(a['points'] for a in added)} core points; objective (new weights) "
              f"{before:.5f} -> {after:.5f}", flush=True)
    final = [(f.signal_cores, f.weight, f.norm, f.signal_mask) for f in joint.forwards]
    for f, (cores, weight, norm, mask) in zip(joint.forwards, saved):
        f.signal_cores, f.weight, f.norm, f.signal_mask = cores, weight, norm, mask
    record["objective_original_weights"] = float(np.sum(joint.residual(z) ** 2))
    for f, (cores, weight, norm, mask) in zip(joint.forwards, final):
        f.signal_cores, f.weight, f.norm, f.signal_mask = cores, weight, norm, mask
    rescored = [(float(np.sum(joint.residual(zz) ** 2)), zz) for _, zz in solutions]
    if not any(np.array_equal(z, zz) for _, zz in solutions):
        rescored.append((float(np.sum(joint.residual(z) ** 2)), z))
    rescored.sort(key=lambda t: t[0])
    record["objective"] = rescored[0][0]
    return rescored[0][1], rescored, record


def main():
    ap = make_parser()
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out), "fit_joint_series")
    prob = build_problem(args)
    series, left_out, lo, hi, obs, model = prob.series, prob.left_out, prob.lo, prob.hi, prob.obs, prob.model
    settings, joint, key_of, centre, xs0, z0 = (prob.settings, prob.joint, prob.key_of, prob.centre, prob.xs0,
                                                prob.z0)
    (ks, mean, sigma), lower, upper = prob.prior, prob.lower, prob.upper
    rng = np.random.default_rng(args.seed)
    centres = [z0]
    if args.seeds:
        seeds = json.load(open(args.seeds))
        rows = [int(np.argmin(np.abs(np.asarray(seeds["x"]) - x))) for x in joint.nodes]
        centres = [np.clip(joint.pack(np.array([[sol[key_of[n]][s] for n in joint.coupling] for s in rows]), xs0),
                           lower + 1e-9, upper - 1e-9) for sol in seeds["solutions"]]
    component_record = {}
    if args.component_search in ("start", "both"):
        # one pass from the first centre guides the multi-start: its result becomes the first centre
        z_cs, component_record["start"] = component_search(joint, prob, centres[0], args, "start")
        if component_record["start"]["objective_after"] < component_record["start"]["objective_before"]:
            centres = [z_cs] + centres
    starts = list(centres) + [
        np.clip(centres[i % len(centres)] + joint.level_shift(rng.normal(0, args.spread, len(z0))), lower + 1e-9, upper - 1e-9)
        for i in range(max(args.starts - len(centres), 0))]
    # random starts drawn from the priors: level of every small coupling ~ N(centre, 2 sigma), direction and size
    # of its change random, shape uniform; 1J levels and the spectrum parameters from the first start
    for _ in range(args.prior_starts):
        z = z0.copy()
        for k, n in enumerate(joint.coupling):
            base_k = k * joint.m
            if k in ks:
                sd = 2.0 * sigma[ks.index(k)]
                level, change = rng.normal(mean[ks.index(k)], sd), rng.normal(0.0, sd)
                if joint.shape == "free":
                    z[base_k: base_k + joint.m] = level + change * np.linspace(0.0, 1.0, joint.nn)
                else:
                    z[base_k], z[base_k + 1] = level, change
                    z[base_k + 2: base_k + joint.m] = 0.0
        starts.append(np.clip(z, lower + 1e-9, upper - 1e-9))
    t0 = time.time()
    global _JOINT_TASK
    schedule = [float(v) for v in args.smoothing.split(",") if v.strip()] if args.smoothing else []
    scan = ({"cycles": args.scan_cycles, "half_width": args.scan_half_width, "step": args.scan_step,
             "seed": args.seed} if args.scan_cycles else None)
    _JOINT_TASK = (joint, lower, upper, args.max_nfev, schedule, scan, args.hold_small_first)
    if args.trace:
        joint.trace = []                          # each start collects its own (worker processes fork this)
    global _MONITOR
    status = None
    if args.monitor == "on":
        from fit_monitor import RunStatus, Tee
        keys = [key_of.get(n, n) for n in joint.coupling]
        status = RunStatus(Path(args.out) / "monitor", sys.argv, len(starts), keys, joint.nodes.tolist(),
                           {key: joint.coupling_values(z0, k).tolist() for k, key in enumerate(keys)})
        sys.stdout = Tee(sys.stdout, status.dir / "console.log")
        _MONITOR = (str(status.dir), keys)
        status.set(phase=f"starts ({args.workers} worker{'s' if args.workers > 1 else ''})")
        print(f"monitor: python scripts/fit_monitor.py {args.out}", flush=True)
    solved = []
    if args.workers > 1:
        import multiprocessing
        with multiprocessing.get_context("fork").Pool(args.workers) as pool:
            for k, out in enumerate(pool.imap(_solve_indexed, list(enumerate(starts)), chunksize=1)):
                solved.append(out)
                print(f"start {k} finished: objective {out[0]:.5f}", flush=True)
                if status is not None:
                    status.finished(k, out[0])
    else:
        for k, z in enumerate(starts):
            solved.append(_solve_indexed((k, z)))
            print(f"start {k} finished: objective {solved[-1][0]:.5f}", flush=True)
            if status is not None:
                status.finished(k, solved[-1][0])
    traces = [t[2] for t in solved] if args.trace else None
    solved = [(t[0], t[1]) for t in solved]
    solutions = []
    for k, (score, x) in enumerate(solved):
        solutions.append((score, x))
        print(f"start {k}: score {score:.5f}", flush=True)
    print(f"{len(starts)} starts in {time.time() - t0:.0f} s", flush=True)
    if status is not None:
        status.set(phase="ranking and writing results")
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
    residual_peak_record = None
    stage_trace = []
    if args.residual_peaks > 0:
        if args.trace:
            joint.trace = stage_trace
        if status is not None:
            from fit_monitor import MonitorWriter
            status.set(phase="residual-peak stage")
            joint.monitor = MonitorWriter(status.dir, "residual_peak_stage", joint, _MONITOR[1])
            joint.trace_label = "residual-peak stage"
        solutions, residual_peak_record = _residual_peak_stage(joint, solutions, lower, upper, args)
        if status is not None:
            joint.monitor.close(solutions[0][0])
            joint.monitor = None
        score, z = solutions[0]
        joint.trace = None
    if args.component_search in ("end", "both"):
        # the best solution again, component by component; a better basin is refit globally and kept if it wins
        if status is not None:
            status.set(phase="component search (end)")
        z_cs, rec = component_search(joint, prob, z, args, "end")
        if rec["objective_after"] < rec["objective_before"]:
            if joint.peaks is not None:
                joint.peak_smooth = args.peak_smooth
            refit = _solve_indexed((len(solved), z_cs))
            if joint.peaks is not None:
                joint.peak_smooth = 0.0
            z_new = refit[1]
            score_new = float(np.sum(joint.residual(z_new) ** 2))
            rec["refit_objective"] = score_new
            rec["refit_kept"] = bool(score_new < score)
            print(f"component search (end): global refit {score_new:.5f} against {score:.5f}"
                  f"{' kept' if score_new < score else ''}", flush=True)
            if score_new < score:
                solutions.insert(0, (score_new, z_new))
                score, z = score_new, z_new
        component_record["end"] = rec
    model_line_record = None
    if args.model_line_passes > 0:
        if status is not None:
            status.set(phase="model-line passes")
        if joint.peaks is not None:
            joint.peak_smooth = args.peak_smooth
        z, solutions, model_line_record = model_line_passes(joint, z, solutions, args.model_line_passes,
                                                            args.model_line_threshold, len(solved))
        if joint.peaks is not None:
            joint.peak_smooth = 0.0
        score = solutions[0][0]
    if args.trace:
        best_start = next((k for k, (_, x) in enumerate(solved)
                           if any(np.array_equal(x, zz) for _, zz in solutions[:1])), None)
        if best_start is None:          # chosen in the residual-peak stage: show its refit after the best start
            best_start = int(np.argmin([sc for sc, _ in solved]))
        _write_trace(Path(args.out), joint, series, traces[best_start] + stage_trace, z, args.trace, key_of,
                     f"start {best_start}")
    r = joint.residual(z)
    jac = joint.jacobian(z)
    dof = max(len(r) - len(z), 1)
    cov = float(r @ r) / dof * np.linalg.pinv(jac.T @ jac)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    start_tables = [{"score": sc, "J_at_x": {key_of.get(n, n): joint.coupling_values(zz, k).tolist()
                                              for k, n in enumerate(joint.coupling)}} for sc, zz in solutions]
    result = {"shape": "monotone (direction and shape free)" if joint.shape == "monotone" else
              "free (J at every concentration independent)", "x": joint.nodes.tolist(),
              "spectra": [{"id": e["id"], "x": e["x"]} for e in series],
              "start_solutions": start_tables,
              "signal_threshold": settings.signal_threshold, "signal_taper_hz": settings.signal_taper_hz,
              "signal_height_power": settings.signal_height_power,
              "peak_penalty": {"strength": args.peak_penalty, "prominence": args.peak_prominence,
                               "tolerance_hz": args.peak_tolerance, "min_sigma": args.peak_min_sigma,
                               "tops": [len(t) for t in joint.peaks] if joint.peaks else [],
                               "dip_strength": args.dip_penalty, "max_width_hz": args.peak_max_width,
                               "smooth": args.peak_smooth, "smooth_scores_sorted_by_hard": smooth_scores},
              "prior": {"sigma_hh": args.prior_sigma_hh, "sigma_ch": args.prior_sigma_ch, "weight": args.prior_weight},
              "scores": [s for s, _ in solutions],
              "data_region_residuals": dict(zip([e["id"] for e in series], joint.data_residuals(z))),
              "seconds": round(time.time() - t0), "couplings": {}, "residual_peaks": residual_peak_record,
              "component_search": component_record, "model_line_passes": model_line_record}
    for k, n in enumerate(joint.coupling):
        block = slice(k * joint.m, (k + 1) * joint.m)
        values = joint.coupling_values(z, k)
        g = joint.coupling_jacobian(z, k)
        value_cov = g @ cov[block, block] @ g.T
        result["couplings"][key_of.get(n, n)] = {
            "J_at_x": values.tolist(), "J_std_at_x": np.sqrt(np.maximum(np.diag(value_cov), 0)).tolist(),
            "direction": ("increasing" if values[-1] > values[0] else "decreasing") if joint.shape == "free" else
                         ("increasing" if z[block][1] > 0 else "decreasing"),
            "change_hz": float(values[-1] - values[0]) if joint.shape == "free" else float(z[block][1]),
            "prior_centre": centre[n] if abs(centre[n]) < 50 else None}
    if left_out is not None:
        result["left_out"] = {"id": left_out["id"], "x": left_out["x"],
                              **predict_left_out(joint, z, left_out, (lo, hi), key_of, settings, real_only=args.real_only)}
    result["shared"] = joint.shared
    shared_values = {n: float(z[joint.ntheta + i]) for i, n in enumerate(joint.shared)}
    result["spectrum_parameters"] = {series[si]["id"]: {**shared_values,
                                                       **{n: float(z[joint.nt + si * joint.nl + i])
                                                          for i, n in enumerate(joint.local)}}
                                     for si in range(joint.ns)}
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
    if status is not None:
        status.set(phase="finished", score=score, finished_at=time.time())


if __name__ == "__main__":
    main()
