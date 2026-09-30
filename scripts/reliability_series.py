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


def build_joint(args):
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    series = json.load(open(args.series))
    lo, hi = (float(v) for v in args.range.split(","))
    obs = [ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), np.load(e["values"]).astype(float),
                                          [(lo, hi)], record=None, real_only=True, label=e["id"]) for e in series]
    spec = json.loads(args.structure)
    spec.setdefault("compound", "series")
    model = build_model(override_couplings(reg.structure_for(spec), json.loads(args.couplings)), ranges=[(lo, hi)])
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    settings = _settings_for(model, args.variant, base)
    joint = JointSeries(model, settings, obs, [e["x"] for e in series])
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = joint.params[0].ties.get(n, n)
            if leader in joint.coupling and leader not in key_of:
                key_of[leader] = key
    return series, joint, settings, key_of


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


def direct_data_score(joint, series, key_of, fit):
    total = 0.0
    for i, (e, f) in enumerate(zip(series, joint.forwards)):
        x = joint.params[i].vector().copy()
        for n in joint.coupling:
            x[joint.col[n]] = fit["couplings"][key_of[n]]["J_at_x"][joint.node_of[i]]
        for n, v in fit["spectrum_parameters"][e["id"]].items():
            x[joint.col[n]] = v
        r = f.predict(x).residual
        total += float(r @ r)
    return total


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
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--signal-threshold", type=float, default=0.0)
    ap.add_argument("--signal-taper", type=float, default=0.0)
    ap.add_argument("--fit", action="append", required=True, help="LABEL=path/fit.json (repeat)")
    ap.add_argument("--tolerance", type=float, default=0.03)
    ap.add_argument("--reliable-hz", type=float, default=1.6)
    ap.add_argument("--trend-hz", type=float, default=1.0)
    ap.add_argument("--literature", default="")
    ap.add_argument("--out", default="runs/processed/reliability")
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out), "reliability_series")
    series, joint, settings, key_of = build_joint(args)
    keys = [key_of[n] for n in joint.coupling]
    norm = float(np.mean([f.norm for f in joint.forwards]))
    fits, pool = {}, []
    for item in args.fit:
        label, path = item.split("=", 1)
        fit = json.load(open(path))
        if ([s["id"] for s in fit["spectra"]] != [e["id"] for e in series]
                or fit["signal_threshold"] != settings.signal_threshold
                or fit["signal_taper_hz"] != settings.signal_taper_hz):
            raise ValueError(f"{label}: other spectra or signal weighting than this comparison")
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
    result = {"criterion": {"score": "weighted data residual sum of squares (prior term removed)",
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
    nc = 5
    nr = int(np.ceil((len(order) + 1) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(17, 3.2 * nr))
    axs = axs.ravel()
    for ax, k in zip(axs, order):
        st, spread = status[k]
        for s in outside.values():
            ax.plot(x, s["J_at_x"][k], color="#bbbbbb", lw=0.8, ls="--", zorder=1)
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
        handles.append(Line2D([], [], color="#bbbbbb", lw=0.8, ls="--",
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
