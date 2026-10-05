"""Which couplings of joint concentration-series fits are determined: compare every stored solution on the data alone.

    python scripts/reliability_series.py --series series.json --structure '{"motif": "pyridine ring", ...}' \\
        --couplings "$(cat centres.json)" --signal-threshold 2.5 --signal-taper 4.0 \\
        --fit 48-start=runs/processed/joint_new_ms48/fit.json --fit staged=runs/processed/joint_real_staged/fit.json \\
        [--tolerance 0.03] [--reliable-hz 1.6] [--trend-hz 1.0] [--literature configs/literature/....json] \\
        [--out runs/processed/reliability]

The fit score of fit_joint_series.py is the weighted data residual sum of squares plus the prior term, and the prior
centres differ between fits, so scores of different fits are not comparable. Here every stored solution (all starts of
every fit) gets its data-only score: score minus its prior term, the prior term being recomputed from the solution's
couplings and the fit's own prior centres, sigmas and weight (JointSeries.set_priors). The data residual depends only
on the data, range, variant and signal weighting, which must be the same for every fit (checked). For each fit's best
solution the data score is also recomputed directly from its couplings and spectrum parameters as a check.

The near-equivalent set pools the solutions of all fits whose data score is within --tolerance of the best data score.
Per coupling: spread = largest (max - min) over the set at any concentration; reliable if spread <= --reliable-hz;
trend only if every solution in the set changes in the same direction by more than --trend-hz from the first to the
last concentration; otherwise not determined. Outputs OUT/reliability.json and OUT/J_trends_reliability.png.

--spectra: the stored start solutions keep only their couplings, so for the best data-score solution the spectrum
parameters (decay rates, delay) are refitted with its couplings held (started from every fit's best spectrum
parameters, the lowest kept; gains, phase and background are linear as in the fit). The refitted data score must
not exceed the stored one. Writes OUT/best_data_spectra.png (experiment, simulation, residual per spectrum).
"""
import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from fit_joint_series import JointSeries                               # noqa: E402
from zulf_core.solver import ObservedSpectrum                          # noqa: E402
from zulf_hypothesis.builder import build_model                        # noqa: E402
from zulf_hypothesis.fit import default_fit_base                       # noqa: E402
from zulf_hypothesis.search import _settings_for                       # noqa: E402

CLASS_COLOR = {"reliable": "#1b7f3b", "trend": "#b8860b", "undetermined": "#777777"}
CLASS_LABEL = {"reliable": "RELIABLE", "trend": "trend only", "undetermined": "not determined"}



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

def build_joint(args):
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    series = json.load(open(args.series))
    lo, hi = (float(v) for v in args.range.split(","))
    obs = [ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), _values(e, getattr(args, "real_only", True)),
                                          [tuple(r) for r in e.get("ranges", [(lo, hi)])], record=e.get("record"),
                                          phasing=e.get("phasing"), real_only=getattr(args, "real_only", True),
                                              label=e["id"]) for e in series]
    spec = json.loads(args.structure)
    spec.setdefault("compound", "series")
    from zulf_hypothesis import exchange_variants
    fragment = exchange_variants(override_couplings(reg.structure_for(spec), json.loads(args.couplings)),
                                 getattr(args, "exchange", "slow"))[0]
    model = build_model(fragment, ranges=[(lo, hi)])
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    if getattr(args, "signal_height_power", 0.0):
        base = dataclasses.replace(base, signal_height_power=args.signal_height_power)
    from fit_joint_series import _rate_policy
    base = _rate_policy(base, args)
    settings = _settings_for(model, args.variant, base)
    shared = [n.strip() for n in getattr(args, "shared", "").split(",") if n.strip()]
    joint = JointSeries(model, settings, obs, [e["x"] for e in series], shared=shared)   # shared: held per refit
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = joint.params[0].ties.get(n, n)
            if leader in joint.coupling and leader not in key_of:
                key_of[leader] = key
    return series, obs, joint, settings, key_of


def prior_term(fit, keys, norm, J):
    """The prior part of a stored score (JointSeries.residual): Gaussian priors on the series average."""
    p = fit["prior"]
    total = 0.0
    for k in keys:
        centre = fit["couplings"][k]["prior_centre"]
        sigma = p["sigma_hh"] if k.startswith("J(H") else p["sigma_ch"]
        if centre is None or sigma <= 0:
            continue
        total += ((np.mean(J[k]) - centre) * np.sqrt(p["weight"]) / sigma / norm) ** 2
    return float(total)


def direct_vector(joint, series, key_of, fit, i):
    """Parameter vector of spectrum i for a fit's stored best (couplings and spectrum parameters)."""
    x = joint.params[i].vector().copy()
    for n in joint.coupling:
        x[joint.col[n]] = fit["couplings"][key_of[n]]["J_at_x"][joint.node_of[i]]
    for n, v in fit["spectrum_parameters"][series[i]["id"]].items():
        x[joint.col[n]] = v
    return x


def direct_data_score(joint, series, key_of, fit):
    """Data residual sum of squares of a fit's stored best, plus its hard missing-peak rows when the fits use them."""
    total = 0.0
    for i, f in enumerate(joint.forwards):
        r = f.predict(direct_vector(joint, series, key_of, fit, i)).residual
        total += float(r @ r)
        if joint.peaks is not None:
            total += float(np.sum(joint._peak_rows(i, r)[0] ** 2))
    return total


def peak_settings(fit):
    """A fit's missing-peak settings (strength 0 = none), without the per-run record fields."""
    p = fit.get("peak_penalty") or {}
    if not p.get("strength") and not p.get("dip_strength"):
        return {"strength": 0.0}
    return {**{k: p[k] for k in ("strength", "prominence", "tolerance_hz", "min_sigma")},
            "dip_strength": p.get("dip_strength", 0.0), "max_width_hz": p.get("max_width_hz", 0.0)}


def refit_spectrum_parameters(joint, series, key_of, J_at_x, fit_files):
    """Spectrum parameters of every spectrum for couplings held at J_at_x; returns per-spectrum parameter vectors."""
    from scipy.optimize import least_squares
    local = [joint.col[n] for n in joint.local]
    vectors = []
    for i, (e, f) in enumerate(zip(series, joint.forwards)):
        lo, hi = (np.asarray(b)[local] for b in joint.params[i].bounds())
        best = None
        for fit in fit_files:
            x = joint.params[i].vector().copy()
            for n in joint.coupling:
                x[joint.col[n]] = J_at_x[key_of[n]][joint.node_of[i]]
            for n, v in fit["spectrum_parameters"][e["id"]].items():
                x[joint.col[n]] = v

            def residual(p, x=x):
                y = x.copy()
                y[local] = p
                return f.predict(y).residual
            sol = least_squares(residual, np.clip(x[local], lo + 1e-9, hi - 1e-9), bounds=(lo, hi), x_scale="jac",
                                max_nfev=200, ftol=1e-12, xtol=1e-12, gtol=1e-12)
            if best is None or sol.cost < best[0]:
                y = x.copy()
                y[local] = sol.x
                best = (sol.cost, y)
        vectors.append(best[1])
    return vectors


def plot_spectra(path, joint, series, obs, vectors, title, reference=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def relative(f, x):
        m = f.data_signal_mask if f.data_signal_mask is not None else np.ones(len(f.y), bool)
        return float(np.linalg.norm(f.mismatch(f.predict(x).model, m)) /
                     max(np.linalg.norm(f.mismatch(np.zeros_like(f.y), m)), 1e-30))
    fig, axes = plt.subplots(len(obs), 1, figsize=(12, 2.3 * len(obs)), sharex=True)
    rows = []
    for s, (ax, o, f) in enumerate(zip(np.atleast_1d(axes), obs, joint.forwards)):
        pred, y = f.predict(vectors[s]).model.real, o.values.real
        m = np.abs(y).max()
        ax.plot(o.frequencies_hz, y / m, color="#222222", lw=0.8, label="experiment")
        ax.plot(o.frequencies_hz, pred / m, color="#c0182a", lw=0.9, label="simulation")
        ax.plot(o.frequencies_hz, (y - pred) / m - 0.35, color="#999999", lw=0.6, label="residual (offset)")
        r = relative(f, vectors[s])
        ref = relative(f, reference[s]) if reference is not None else None
        rows.append((series[s]["id"], r, ref))
        note = f" (lowest fit score solution {ref:.3f})" if ref is not None else ""
        ax.set_title(f"{series[s]['id']} (x = {joint.xs[s]:.2f}): relative residual {r:.3f}{note}", fontsize=9,
                     loc="left")
        ax.set_yticks([])
    np.atleast_1d(axes)[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
    np.atleast_1d(axes)[-1].set_xlabel("frequency (Hz)")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return rows


def classify(values, reliable_hz, trend_hz):
    a = np.asarray(values)
    spread = float((a.max(axis=0) - a.min(axis=0)).max())
    change = a[:, -1] - a[:, 0]
    if spread <= reliable_hz:
        return "reliable", spread
    if len({np.sign(c) for c in change}) == 1 and np.abs(change).min() > trend_hz:
        return "trend", spread
    return "undetermined", spread


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--real-only", type=_flag, default=True, help="true (default) or false: complex spectra")
    ap.add_argument("--shared", default="", help="spectrum parameters shared by all spectra, as in fit_joint_series")
    ap.add_argument("--exchange", default="slow", choices=["slow", "fast"], help="as in fit_joint_series")
    ap.add_argument("--family-edges", default="", help="as in fit_joint_series")
    ap.add_argument("--rate-bounds", default="", help="as in fit_joint_series")
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--signal-threshold", type=float, default=0.0)
    ap.add_argument("--signal-taper", type=float, default=0.0)
    ap.add_argument("--signal-height-power", type=float, default=0.0,
                    help="extra weight (local peak height / max)^-power: small peaks count more (0 = off)")
    ap.add_argument("--fit", action="append", required=True, help="LABEL=path/fit.json (repeat)")
    ap.add_argument("--tolerance", type=float, default=0.03)
    ap.add_argument("--reliable-hz", type=float, default=1.6)
    ap.add_argument("--trend-hz", type=float, default=1.0)
    ap.add_argument("--literature", default="")
    ap.add_argument("--spectra", action="store_true", help="refit and plot the spectra of the best data-score solution")
    ap.add_argument("--out", default="runs/processed/reliability")
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out), "reliability_series")
    series, obs, joint, settings, key_of = build_joint(args)
    keys = [key_of[n] for n in joint.coupling]
    norm = float(np.mean([f.norm for f in joint.forwards]))
    # fits with missing-peak rows: their stored scores contain the hard rows, so the compared score is data plus
    # those rows (the same settings in every fit; checked below)
    peaks = peak_settings(json.load(open(args.fit[0].split("=", 1)[1])))
    if peaks["strength"] or peaks.get("dip_strength"):
        joint.set_peak_penalty(peaks["strength"], peaks["prominence"], peaks["tolerance_hz"],
                               min_sigma=peaks["min_sigma"], smooth=0.0, dips=peaks.get("dip_strength", 0.0),
                               max_width_hz=peaks.get("max_width_hz") or None)
    fits, pool = {}, []
    for item in args.fit:
        label, path = item.split("=", 1)
        fit = json.load(open(path))
        if ([s["id"] for s in fit["spectra"]] != [e["id"] for e in series]
                or fit["signal_threshold"] != settings.signal_threshold
                or fit["signal_taper_hz"] != settings.signal_taper_hz
                or fit.get("signal_height_power", 0.0) != settings.signal_height_power):
            raise ValueError(f"{label}: other spectra or signal weighting than this comparison")
        if peak_settings(fit) != peaks:
            raise ValueError(f"{label}: other missing-peak settings than the first fit")
        sols = []
        for i, s in enumerate(fit["start_solutions"]):
            prior = prior_term(fit, keys, norm, s["J_at_x"])
            sols.append({"fit": label, "rank_in_fit": i, "score": s["score"], "prior": prior,
                         "data_score": s["score"] - prior, "J_at_x": s["J_at_x"]})
        best = min(sols, key=lambda s: s["score"])
        direct = direct_data_score(joint, series, key_of, fit)
        if abs(direct - best["data_score"]) > 1e-6 * max(direct, 1e-12):
            raise RuntimeError(f"{label}: data score {best['data_score']} != direct {direct}")
        fits[label] = {"path": path, "starts": len(sols), "best_score": best["score"], "best_prior": best["prior"],
                       "best_data_score": best["data_score"], "min_data_score": min(s["data_score"] for s in sols)}
        pool += sols
    best_data = min(s["data_score"] for s in pool)
    chosen = sorted([s for s in pool if s["data_score"] <= (1 + args.tolerance) * best_data],
                    key=lambda s: s["data_score"])
    for f in fits.values():
        f["relative_to_best_data"] = f["min_data_score"] / best_data - 1
        f["in_set"] = sum(1 for s in chosen if fits.get(s["fit"]) is f)
    x = json.load(open(fits[next(iter(fits))]["path"]))["x"]
    status = {k: classify([s["J_at_x"][k] for s in chosen], args.reliable_hz, args.trend_hz) for k in keys}
    result = {"criterion": {"score": "weighted data residual sum of squares (prior term removed)" +
                            (" plus the hard missing-peak rows" if peaks["strength"] else ""),
                            "peak_penalty": peaks,
                            "tolerance": args.tolerance, "reliable_hz": args.reliable_hz, "trend_hz": args.trend_hz},
              "x": x, "best_data_score": best_data, "fits": fits,
              "set": [{k: s[k] for k in ("fit", "rank_in_fit", "score", "data_score")} for s in chosen],
              "couplings": {k: {"class": c, "spread_hz": sp, "best_J_at_x": chosen[0]["J_at_x"][k]}
                            for k, (c, sp) in status.items()}}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reliability.json").write_text(json.dumps(result, indent=1))
    figure = plot(out / "J_trends_reliability.png", x, keys, status, chosen, pool, fits, best_data, args)
    run_log.finish({"best_data_score": best_data, "set_size": len(chosen),
                    "fits": {k: round(v["relative_to_best_data"], 4) for k, v in fits.items()},
                    "classes": {k: c for k, (c, _) in status.items()}}, figures=[figure])
    if args.spectra:
        fit_files = [json.load(open(f["path"])) for f in fits.values()]
        vectors = refit_spectrum_parameters(joint, series, key_of, chosen[0]["J_at_x"], fit_files)
        refit = sum(float(r @ r) for r in (f.predict(v).residual for f, v in zip(joint.forwards, vectors)))
        if refit > chosen[0]["data_score"] * (1 + 1e-6):
            raise RuntimeError(f"refitted data score {refit} above the stored {chosen[0]['data_score']}")
        lowest = min((s for s in pool), key=lambda s: s["score"])
        owner = next(f for f in fit_files if min(s["score"] for s in f["start_solutions"]) == lowest["score"])
        reference = [direct_vector(joint, series, key_of, owner, i) for i in range(len(series))]
        rows = plot_spectra(out / "best_data_spectra.png", joint, series, obs, vectors,
                            f"Best data-score solution ({chosen[0]['fit']}, data score {refit:.4f}); couplings as "
                            f"stored, decay rates and delays refitted", reference)
        result["best_data_spectra"] = {"data_score_refitted": refit,
                                       "relative_residuals": {i: r for i, r, _ in rows},
                                       "lowest_fit_score_relative_residuals": {i: r for i, _, r in rows}}
        (out / "reliability.json").write_text(json.dumps(result, indent=1))
        print(f"best data-score solution: refitted data score {refit:.5f} (stored {chosen[0]['data_score']:.5f})")
        for i, r, ref in rows:
            print(f"  {i}: relative residual {r:.3f} (lowest fit score solution {ref:.3f})")
    for label, f in fits.items():
        print(f"{label}: {f['starts']} solutions, best score {f['best_score']:.4f} = data {f['best_data_score']:.4f} "
              f"+ prior {f['best_prior']:.4f}; best data score {f['min_data_score']:.4f} "
              f"({100 * f['relative_to_best_data']:+.1f} %), in set {f['in_set']}")
    for k in sorted(keys, key=lambda k: (list(CLASS_LABEL).index(status[k][0]), status[k][1])):
        print(f"{k:12s} {CLASS_LABEL[status[k][0]]:15s} spread {status[k][1]:5.2f} Hz")


def plot(path, x, keys, status, chosen, pool, fits, best_data, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    lit = {}
    if args.literature:
        L = json.load(open(args.literature))
        lit = dict(L["carbon_proton"]["absolute_values"])
        lit.update(L["proton_proton"]["values"])
    order = sorted(keys, key=lambda k: (list(CLASS_LABEL).index(status[k][0]), status[k][1]))
    outside = {}
    for s in pool:                                   # best solution of every fit with none in the set
        if not any(c["fit"] == s["fit"] for c in chosen):
            if s["fit"] not in outside or s["data_score"] < outside[s["fit"]]["data_score"]:
                outside[s["fit"]] = s
    palette = ["#1f5fa8", "#2a9d8f", "#8a5cb8", "#e07b39", "#6b7280"]
    colour = {label: palette[i % len(palette)] for i, label in enumerate(outside)}
    nc = 5
    nr = int(np.ceil((len(order) + 1) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(17, 3.2 * nr))
    axs = axs.ravel()
    for ax, k in zip(axs, order):
        st, spread = status[k]
        for label, s in outside.items():
            ax.plot(x, s["J_at_x"][k], color=colour[label], lw=1.1, ls="--", zorder=1)
        for s in chosen[1:]:
            ax.plot(x, s["J_at_x"][k], color="#f2a7b0", lw=0.9, zorder=2)
        ax.plot(x, chosen[0]["J_at_x"][k], "o-", color="#c0182a", lw=2.4, ms=4, zorder=3)
        if k in lit:
            ax.axhline(lit[k], color="k", ls=":", lw=1)
        ax.set_title(f"{k}: {CLASS_LABEL[st]} (spread {spread:.1f} Hz)", fontsize=9, color=CLASS_COLOR[st],
                     fontweight="bold" if st != "undetermined" else "normal")
        if st != "undetermined":
            ax.set_facecolor("#eaf6ec" if st == "reliable" else "#fbf4e2")
        for sp in ax.spines.values():
            sp.set_edgecolor(CLASS_COLOR[st])
            sp.set_linewidth(2 if st != "undetermined" else 0.8)
        ax.tick_params(labelsize=7)
    for ax in axs[len(order):]:
        ax.axis("off")
    handles = [Line2D([], [], color="#c0182a", lw=2.4, marker="o",
                      label=f"best data score ({chosen[0]['fit']}, {best_data:.4f})"),
               Line2D([], [], color="#f2a7b0", lw=0.9,
                      label=f"other solutions within {100 * args.tolerance:.0f} % ({len(chosen) - 1})")]
    for label, s in outside.items():
        handles.append(Line2D([], [], color=colour[label], lw=1.1, ls="--",
                              label=f"{label} best, not in set ({100 * (s['data_score'] / best_data - 1):+.1f} %)"))
    if lit:
        handles.append(Line2D([], [], color="k", ls=":", lw=1, label="literature"))
    axs[len(order)].legend(handles=handles, loc="center", fontsize=9, frameon=False)
    fig.suptitle(f"Couplings (Hz) vs mole fraction. Set: all stored solutions of {len(fits)} joint fits with weighted "
                 f"data score (prior removed) within {100 * args.tolerance:.0f} % of the best ({len(chosen)} solutions).\n"
                 f"Spread = largest range over x. RELIABLE: spread <= {args.reliable_hz:g} Hz; trend only: every "
                 f"solution changes the same way by > {args.trend_hz:g} Hz", fontsize=10)
    fig.supxlabel("mole fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return str(path)


if __name__ == "__main__":
    main()
