"""Line-shape probe: where the data want the model lines in one window, without the couplings (D69).

    python scripts/line_shape_probe.py RUN_FOLDER --window 129.8,133.8 [--components 0,2] [--groups 8]
        [--max-shift 0.6] [--rate-min 0.1] [--rate-up 1.5] [--merge 0.06] [--spectrum 0] [--max-nfev 300]
        [--out OUT.json] [--figure OUT.png]

A diagnostic, not a fit of the structure. From a finished fit (RUN_FOLDER/fit.json and the command recorded in
RUN_FOLDER/monitor/status.json) the transitions of the chosen components are edited inside the forward model:
the `--groups` strongest line groups in the window (lines closer than `--merge` Hz form one group) each get a free
frequency shift within +-`--max-shift` Hz, and the decay-rate families of those components that hold lines in the
window are free between `--rate-min` 1/s and `--rate-up` times their fitted value. Everything else (couplings,
the other lines, phase, field) stays at the fit; the gains are solved again as in the fit. The shifts and rates
are fitted (least squares) to the residual data - model in the window (real and imaginary parts unless the series
is real only). Reported: every moved group (component, position, summed relative amplitude, shift), the rates
before and after, the rms of the real residual in the window (units of the largest real data value there) and the
whole objective of the fit (hard missing-peak rows) before and after.

The result says which lines would have to sit where, and how wide, for the window to match: e.g. DMEDA
130.9 / 131.75 / 132.75 Hz (2026-10-16), where the 13C@C1 lines had to gather into three narrow clusters that no
coupling set of the structure produces. The moved lines do not correspond to coupling values; conditional
numerical result of the given structure and data.
"""
import argparse
import contextlib
import json
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def line_groups(frequencies_hz, amplitudes, window, merge_hz=0.06):
    """Lines inside window (lo, hi) grouped where neighbours are closer than merge_hz: list of index arrays,
    in frequency order."""
    lo, hi = window
    idx = np.flatnonzero((frequencies_hz >= lo) & (frequencies_hz <= hi))
    idx = idx[np.argsort(frequencies_hz[idx])]
    groups = []
    for i in idx:
        if groups and frequencies_hz[i] - frequencies_hz[groups[-1][-1]] < merge_hz:
            groups[-1].append(int(i))
        else:
            groups.append([int(i)])
    return [np.array(g, int) for g in groups]


@contextlib.contextmanager
def edited_lines(forward, edits):
    """Inside the block, forward.transitions returns the transitions of component c with the lines edits[c] =
    [(indices, shift_hz), ...] moved by their shifts (amplitudes and families unchanged). The shifts are read at
    call time, so a caller may change them in place between predictions."""
    from zulf_core.physics.transitions import TransitionList
    original = forward.transitions

    def transitions(self, values, c, system):
        tl = original(values, c, system)
        if c not in edits:
            return tl
        fr = tl.frequencies_hz.copy()
        for idx, shift in edits[c]:
            fr[idx] += float(shift[0]) if np.ndim(shift) else float(shift)
        return TransitionList(fr, tl.amplitudes, tl.dc, tl.families, dict(tl.metadata), tl.line_rates)

    forward.transitions = types.MethodType(transitions, forward)
    try:
        yield
    finally:
        forward.transitions = original


def probe(joint, z, window, s=0, components=None, groups=8, max_shift=0.6, rate_min=0.1, rate_up=1.5,
          merge_hz=0.06, max_nfev=300, min_relative=0.05):
    """Fit free shifts of the strongest line groups and free rate families in window (see the module docstring).
    Returns a dict with the moved groups, rates, rms and objective before and after, and the model curves."""
    from scipy.optimize import least_squares
    fw = joint.forwards[s]
    P = joint.params[s]
    values = P.values(joint.spectrum_vector(z, s))
    systems = P.systems(values)
    lo, hi = window
    lists = {c: fw.transitions(values, c, system) for c, system in enumerate(systems)}
    if components is None:
        components = []
        for c, tl in lists.items():
            amp = np.abs(tl.amplitudes)
            inside = (tl.frequencies_hz >= lo) & (tl.frequencies_hz <= hi)
            if amp.size and np.any(inside & (amp >= min_relative * amp.max())):
                components.append(c)
    candidates = []
    for c in components:
        tl = lists[c]
        amp = np.abs(tl.amplitudes)
        for g in line_groups(tl.frequencies_hz, amp, window, merge_hz):
            candidates.append((float(amp[g].sum() / amp.max()), c, g, float(np.average(tl.frequencies_hz[g],
                                                                                       weights=amp[g]))))
    chosen = sorted(sorted(candidates, key=lambda t: -t[0])[:groups], key=lambda t: t[3])
    shifts = [np.zeros(1) for _ in chosen]
    edits = {}
    for (_, c, g, _), sh in zip(chosen, shifts):
        edits.setdefault(c, []).append((g, sh))
    edges = np.asarray(P.policy.family_edges_hz, float)
    rate_names = []
    for c in components:
        fr = lists[c].frequencies_hz
        for q in sorted({int(np.searchsorted(edges, v, side="right")) for v in fr[(fr >= lo) & (fr <= hi)]}):
            n = f"c{c}.log_rate{q}"
            if n in joint.local or n in joint.shared:
                rate_names.append(n)
    ridx = np.array([joint.ntheta + joint.shared.index(n) if n in joint.shared else
                     joint.nt + s * joint.nl + joint.local.index(n) for n in rate_names], int)
    f = np.asarray(fw.f, float)
    y = np.asarray(fw.y)
    win = (f >= lo) & (f <= hi)
    scale = float(np.max(np.abs(np.real(y[win])))) or 1.0
    ng = len(chosen)

    def z_of(p):
        zz = z.copy()
        zz[ridx] = p[ng:]
        return zz

    def model(zz):
        return np.asarray(fw.predict(joint.spectrum_vector(zz, s)).model)

    def residual(p):
        for k in range(ng):
            shifts[k][0] = p[k]
        r = (y - model(z_of(p)))[win] / scale
        return r.real if fw.real_only or not np.iscomplexobj(r) else np.r_[r.real, r.imag]

    rms = lambda m: float(np.sqrt(np.mean(np.real(y - m)[win] ** 2)) / scale)
    obj = lambda zz: float(np.sum(joint.residual(zz) ** 2))
    with edited_lines(fw, edits):
        p0 = np.r_[np.zeros(ng), z[ridx]]
        lower = np.r_[np.full(ng, -max_shift), np.full(len(ridx), np.log(rate_min))]
        upper = np.r_[np.full(ng, max_shift), z[ridx] + np.log(rate_up)]
        p0 = np.clip(p0, lower + 1e-9, upper - 1e-9)
        m0 = model(z)
        before = (rms(m0), obj(z))
        sol = least_squares(residual, p0, bounds=(lower, upper), diff_step=0.01, max_nfev=max_nfev)
        for k in range(ng):
            shifts[k][0] = sol.x[k]
        z1 = z_of(sol.x)
        m1 = model(z1)
        after = (rms(m1), obj(z1))
        moved = [{"component": c, "frequency_hz": fq, "relative_amplitude": a, "lines": int(len(g)),
                  "shift_hz": float(sol.x[k]), "new_frequency_hz": fq + float(sol.x[k])}
                 for k, (a, c, g, fq) in enumerate(chosen)]
        sticks_after = {c: lists[c].frequencies_hz.copy() for c in components}
        for (_, c, g, _), sh in zip(chosen, shifts):
            sticks_after[c][g] += sh[0]
    return {"window_hz": [lo, hi], "spectrum": s, "components": components, "groups": moved,
            "rates": {n: [float(np.exp(z[i])), float(np.exp(sol.x[ng + k]))] for k, (n, i) in
                      enumerate(zip(rate_names, ridx))},
            "rms_window": [before[0], after[0]], "objective": [before[1], after[1]],
            "settings": {"groups": groups, "max_shift_hz": max_shift, "rate_min": rate_min, "rate_up": rate_up,
                         "merge_hz": merge_hz, "max_nfev": max_nfev},
            "_curves": {"f": f, "data": y, "before": m0, "after": m1, "scale": scale,
                        "sticks_before": {c: (lists[c].frequencies_hz, np.abs(lists[c].amplitudes))
                                          for c in components},
                        "sticks_after": {c: (sticks_after[c], np.abs(lists[c].amplitudes)) for c in components}}}


def figure(result, path, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cv = result["_curves"]
    lo, hi = result["window_hz"]
    pad = max(2.0, 0.5 * (hi - lo))
    f, y, sc = cv["f"], cv["data"], cv["scale"]
    k = (f >= lo - pad) & (f <= hi + pad)
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [1.4, 1]})
    ax[0].plot(f[k], y.real[k] / sc, color="black", lw=1.0, label="experiment")
    ax[0].plot(f[k], cv["before"].real[k] / sc, color="#d4682e", lw=1.3,
               label=f"fit (window rms {result['rms_window'][0]:.3f})")
    ax[0].plot(f[k], cv["after"].real[k] / sc, color="#3d7cc9", lw=1.3,
               label=f"lines moved, rates free (window rms {result['rms_window'][1]:.3f})")
    ax[0].plot(f[k], (y.real[k] - cv["after"].real[k]) / sc - 0.95, color="0.55", lw=0.7)
    ax[0].axvspan(lo, hi, color="#24476e", alpha=0.05, lw=0)
    ax[0].legend(fontsize=8, frameon=False, loc="upper left")
    colours = ["#d4682e", "#2e8b57", "#9a6fb0", "#c8553d"]
    for n, c in enumerate(result["components"]):
        fb, ab = cv["sticks_before"][c]
        fa, aa = cv["sticks_after"][c]
        top = ab.max() or 1.0
        sel = (fb >= lo - pad) & (fb <= hi + pad)
        ax[1].vlines(fb[sel], 0, ab[sel] / top, color=colours[n % 4], lw=1.1, label=f"component {c}, fit")
        sel = (fa >= lo - pad) & (fa <= hi + pad)
        ax[1].vlines(fa[sel], -aa[sel] / top, 0, color="#3d7cc9", lw=1.1)
    ax[1].axhline(0, color="0.7", lw=0.6)
    ax[1].set_yticks([])
    ax[1].set_xlabel("frequency (Hz)")
    ax[1].text(0.01, 0.92, "fit lines (up), moved lines (down)", transform=ax[1].transAxes, fontsize=9, color="0.3")
    if result["components"]:
        ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    for a in ax:
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.suptitle((title + ": " if title else "") + "line-shape probe (diagnostic, not a fit of the structure); "
                 f"objective {result['objective'][0]:.4f} -> {result['objective'][1]:.4f}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--window", required=True, help="lo,hi (Hz)")
    ap.add_argument("--components", default="", help="component indices (default: those with lines in the window)")
    ap.add_argument("--groups", type=int, default=8)
    ap.add_argument("--max-shift", type=float, default=0.6)
    ap.add_argument("--rate-min", type=float, default=0.1)
    ap.add_argument("--rate-up", type=float, default=1.5)
    ap.add_argument("--merge", type=float, default=0.06)
    ap.add_argument("--spectrum", type=int, default=0)
    ap.add_argument("--max-nfev", type=int, default=300)
    ap.add_argument("--out", default="")
    ap.add_argument("--figure", default="")
    args = ap.parse_args()
    from fit_monitor import SpectrumCache
    joint = SpectrumCache().get(args.run).joint
    z = np.array(json.load(open(Path(args.run) / "fit.json"))["z_final"], float)
    lo, hi = (float(v) for v in args.window.split(","))
    comps = [int(c) for c in args.components.split(",")] if args.components else None
    res = probe(joint, z, (lo, hi), s=args.spectrum, components=comps, groups=args.groups, max_shift=args.max_shift,
                rate_min=args.rate_min, rate_up=args.rate_up, merge_hz=args.merge, max_nfev=args.max_nfev)
    print(f"window {lo}-{hi} Hz, components {res['components']}: rms {res['rms_window'][0]:.4f} -> "
          f"{res['rms_window'][1]:.4f}, objective {res['objective'][0]:.5f} -> {res['objective'][1]:.5f}")
    for g in res["groups"]:
        print(f"  c{g['component']} {g['frequency_hz']:8.3f} Hz (amp {g['relative_amplitude']:.2f}, {g['lines']} lines)"
              f" -> {g['new_frequency_hz']:8.3f} Hz ({g['shift_hz']:+.3f})")
    for n, (a, b) in res["rates"].items():
        print(f"  {n}: {a:.2f} -> {b:.2f} 1/s")
    if args.figure:
        figure(res, args.figure, title=Path(args.run).name)
    if args.out:
        out = {k: v for k, v in res.items() if not k.startswith("_")}
        out["run"] = args.run
        json.dump(out, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
