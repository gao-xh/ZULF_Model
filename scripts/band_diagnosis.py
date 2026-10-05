"""Which part of a fit to adjust: misfit per isotopologue band, and the parameters that would lower it.

    python scripts/band_diagnosis.py <the fit_joint_series options of the fit> --fit RUN/fit.json \\
        [--band 186,206] [--bands 3] [--top 10] [--local-fit 3] [--figure OUT.png] [--json OUT.json]

1. Bands. The lines of every isotopologue (relative amplitude >= --min-line of that component's strongest line)
   are clustered at gaps > --gap Hz. A band is one cluster +- 1 Hz. Its cost is the fitter's weighted residual
   sum of squares on the band's points, reported with its share of the spectrum cost and its relative residual
   |data - model| / |data| (weighted). Bands are ranked by cost. --band lo,hi diagnoses one window instead.
2. Levers. At a converged fit no single parameter lowers the whole objective, so a band that stays bad is held
   there by the rest of the spectrum or lacks a parameter of its own. For the worst --bands bands, every free
   parameter p (each coupling at this node, each decay-rate family, the delay) is examined with its Jacobian column
   c split into band rows c_b and the other rows c_o (other bands, missing-peak rows, priors):
     wanted step: the one-parameter Gauss-Newton step of the band alone, dp = -(c_b . r_b) / (c_b . c_b), limited
       to a plausible size (--max-step-hz for couplings, 2x that for one-bond couplings, 0.7 for a log decay
       rate, 0.2 ms for the delay);
     band gain: the band cost that step removes (first order in the model);
     cost elsewhere: the change of the other rows for the same step;
     selectivity: c_b . c_b / c . c, the share of the parameter's effect that falls on the band.
   Each lever gets a verdict: "conflict" (the band gains, the rest loses about as much: the band and the rest want
   different values, so the model is inconsistent in this parameter; look at what else it controls), "free knob"
   (the band gains, the rest barely changes: an unconverged or new direction, refit it), or "weak" (gain below
   --min-gain of the band cost). When every lever is weak, no free parameter acts on this misfit: the model lacks
   something there (a line, a coupling held fixed or tied, another structure).
3. Test (--local-fit N > 0). The N strongest levers by band gain are refined together on the whole objective
   (all rows, everything else held), so the result shows what they can do for the band without hurting the rest.

Nothing here changes the fit; the local-fit vector is only a suggestion (write it with --save-start as a
--from-joint file to refit globally).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import fit_joint_series as fj                                            # noqa: E402
from j_tuner import TuningSession, load_fit                              # noqa: E402


def band_rows(f, lo, hi, n_rows):
    sel = np.flatnonzero((f.f >= lo) & (f.f <= hi))
    nf = len(f.f)
    return np.r_[sel, sel + nf] if n_rows == 2 * nf else sel


def parameter_names(joint, prob, s):
    """Label of every entry of z that belongs to spectrum s (couplings at its node, shared and local parameters)."""
    names = {}
    for k, n in enumerate(joint.coupling):
        names[k * joint.m + (joint.node_of[s] if joint.shape == "free" else 0)] = ("coupling", prob.key_of.get(n, n))
    for i, n in enumerate(joint.shared):
        names[joint.ntheta + i] = ("spectrum", n)
    for i, n in enumerate(joint.local):
        names[joint.nt + s * joint.nl + i] = ("spectrum", n)
    return names


def rate_label(name, labels, edges):
    """'c1.log_rate5' -> '13C@C2 rate family 5 (189.5-192.6 Hz)'."""
    if ".log_rate" not in name:
        return name
    c, k = name.split(".log_rate")
    comp = labels[int(c[1:])] if c[1:].isdigit() and int(c[1:]) < len(labels) else c
    k = int(k)
    if edges:
        lo = edges[k - 1] if k > 0 else None
        hi = edges[k] if k < len(edges) else None
        span = f"{'' if lo is None else f'{lo:g}'}-{'' if hi is None else f'{hi:g}'} Hz"
        return f"{comp} rate family {k} ({span})"
    return f"{comp} rate {k}"


def find_bands(session, z, lo, hi, min_line, gap):
    labels = session.joint.model.component_labels
    lines = session.lines(z, (lo, hi), min_relative=min_line)
    bands = []
    for c, label in enumerate(labels):
        fr = np.sort([l["frequency_hz"] for l in lines if l["component"] == c])
        if not len(fr):
            continue
        start = fr[0]
        for a, b in zip(fr[:-1], fr[1:]):
            if b - a > gap:
                bands.append({"component": label, "lo": start - 1.0, "hi": a + 1.0})
                start = b
        bands.append({"component": label, "lo": start - 1.0, "hi": fr[-1] + 1.0})
    return bands


def linearize(prob, z, s=0):
    """Residual, Jacobian and parameter bookkeeping of the whole objective at z (shared by every band)."""
    joint = prob.joint
    names = parameter_names(joint, prob, s)
    idx = np.array(sorted(names))
    for fw in joint.forwards:
        fw.jacobian_only = None
    r = joint.residual(z)
    jac = joint.jacobian(z)
    offset = 0                           # rows of spectrum s start after the rows of the spectra before it
    for t in range(s):
        rt = joint.forwards[t].predict(joint.spectrum_vector(z, t)).residual
        offset += len(rt)
        if joint.peaks is not None:
            offset += len(joint._peak_rows(t, rt)[0])
    nrow_s = len(joint.forwards[s].predict(joint.spectrum_vector(z, s)).residual)
    keys = [prob.key_of.get(n, n) for n in joint.coupling]
    hz_per_z, value = {}, {}
    for p in idx:
        kind, name = names[p]
        if kind == "coupling":
            k = keys.index(name)
            dz = np.zeros_like(z)
            dz[p] = 1e-6
            node = joint.node_of[s]
            hz_per_z[p] = float((joint.coupling_values(z + dz, k) - joint.coupling_values(z, k))[node] / 1e-6)
            value[p] = float(joint.coupling_values(z, k)[node])
    return {"r": r, "jac": jac, "names": names, "idx": idx, "offset": offset, "nrow_s": nrow_s,
            "hz_per_z": hz_per_z, "value": value}


def band_levers(prob, z, s, lo, hi, max_step_hz=1.0, min_gain=0.03, edges=(), lin=None):
    """(band cost, levers sorted by band gain) for the window lo-hi of spectrum s; see the module docstring."""
    joint = prob.joint
    lin = lin or linearize(prob, z, s)
    r, jac, names, hz_per_z = lin["r"], lin["jac"], lin["names"], lin["hz_per_z"]
    labels = joint.model.component_labels
    rows = lin["offset"] + band_rows(joint.forwards[s], lo, hi, lin["nrow_s"])
    other = np.setdiff1d(np.arange(len(r)), rows)
    rb, ro = r[rows], r[other]
    cost_b = float(rb @ rb)

    def limit(kind, name, dp, p):
        if kind == "coupling" and hz_per_z.get(p):
            cap = max_step_hz * (2.0 if abs(lin["value"][p]) > 50 else 1.0) / abs(hz_per_z[p])
            return float(np.clip(dp, -cap, cap))
        if ".log_rate" in name:
            return float(np.clip(dp, -0.7, 0.7))
        if name == "phase_delay":
            return float(np.clip(dp, -2e-4, 2e-4))
        return dp

    levers = []
    for p in lin["idx"]:
        c_b, c_o = jac[rows, p], jac[other, p]
        cc_b, cc_o = float(c_b @ c_b), float(c_o @ c_o)
        if cc_b <= 0:
            continue
        kind, name = names[p]
        dp = limit(kind, name, -float(c_b @ rb) / cc_b, p)
        gain = -(2 * dp * float(c_b @ rb) + dp * dp * cc_b)
        elsewhere = 2 * dp * float(c_o @ ro) + dp * dp * cc_o
        share = gain / cost_b if cost_b else 0.0
        if share < min_gain:
            verdict = "weak"
        elif elsewhere > 0.5 * gain:
            verdict = "conflict"
        else:
            verdict = "free knob"
        levers.append({"kind": kind, "name": name,
                       "label": name if kind == "coupling" else rate_label(name, labels, edges),
                       "band_gain": gain, "band_share": share, "elsewhere": elsewhere,
                       "elsewhere_share": elsewhere / cost_b if cost_b else 0.0,
                       "selectivity": cc_b / (cc_b + cc_o), "verdict": verdict, "step": dp,
                       "step_hz": dp * hz_per_z[p] if p in hz_per_z else None, "index": int(p)})
    levers.sort(key=lambda d: -d["band_gain"])
    return cost_b, levers


def main():
    ap = fj.make_parser()
    ap.description = __doc__
    ap.add_argument("--fit", required=True)
    ap.add_argument("--band", default="", help="lo,hi: diagnose this window only")
    ap.add_argument("--bands", type=int, default=2, help="number of worst bands to diagnose")
    ap.add_argument("--min-line", type=float, default=0.05)
    ap.add_argument("--gap", type=float, default=2.0)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--local-fit", type=int, default=3, help="levers refitted together on the whole objective (0 = no)")
    ap.add_argument("--max-step-hz", type=float, default=1.0, help="limit of a coupling's wanted step")
    ap.add_argument("--min-gain", type=float, default=0.03, help="band-gain share below which a lever is weak")
    ap.add_argument("--spectrum", type=int, default=0)
    ap.add_argument("--figure", default="")
    ap.add_argument("--json", default="")
    ap.add_argument("--save-start", default="", help="write the local-fit result as a fit.json-like start")
    args = ap.parse_args()
    args.monitor = "off"
    prob = fj.build_problem(args)
    fit = json.load(open(args.fit))
    missing = load_fit(prob, fit)
    if missing:
        print("not in the fit (start values kept):", ", ".join(missing))
    joint, s = prob.joint, args.spectrum
    joint.set_smoothing(0.0)
    joint.peak_smooth = 0.0
    session = TuningSession(prob, s)
    z = prob.z0.copy()
    f = joint.forwards[s]
    pred = f.predict(joint.spectrum_vector(z, s))
    r_s = np.asarray(pred.residual)
    spectrum_cost = float(r_s @ r_s)
    weighted_y = np.asarray(f.y) * f.weight / f.norm
    if len(weighted_y) != len(r_s):                # complex data, residual rows stacked as real then imaginary
        weighted_y = np.r_[weighted_y.real, weighted_y.imag]
    lo_all, hi_all = float(f.f.min()), float(f.f.max())
    edges = [float(v) for v in args.family_edges.split(",")] if getattr(args, "family_edges", "") else []
    labels = joint.model.component_labels

    if args.band:
        lo, hi = (float(v) for v in args.band.split(","))
        bands = [{"component": "window", "lo": lo, "hi": hi}]
    else:
        bands = find_bands(session, z, lo_all, hi_all, args.min_line, args.gap)
    for b in bands:
        rows = band_rows(f, b["lo"], b["hi"], len(r_s))
        b["points"] = int(len(rows) // (2 if len(r_s) == 2 * len(f.f) else 1))
        b["cost"] = float(r_s[rows] @ r_s[rows]) if len(rows) else 0.0
        b["share"] = b["cost"] / spectrum_cost if spectrum_cost > 0 else 0.0
        b["relative"] = (float(np.linalg.norm(r_s[rows]) / max(np.linalg.norm(weighted_y[rows]), 1e-30))
                         if len(rows) else 0.0)
    bands = [b for b in bands if b["points"] > 0]
    bands.sort(key=lambda b: -b["cost"])
    objective = float(np.sum(joint.residual(z) ** 2))
    print(f"objective {objective:.5f}; spectrum cost {spectrum_cost:.5f}")
    print("bands (worst first):")
    for b in bands:
        print(f"  {b['component']:14s} {b['lo']:7.2f}-{b['hi']:7.2f} Hz  cost {b['cost']:.5f} "
              f"({100 * b['share']:4.1f} % of the spectrum)  relative residual {b['relative']:.3f}")

    lin = linearize(prob, z, s)
    report = {"fit": args.fit, "objective": objective, "spectrum_cost": spectrum_cost, "bands": bands,
              "diagnosed": []}
    for b in bands[:args.bands]:
        cost_b, levers = band_levers(prob, z, s, b["lo"], b["hi"], args.max_step_hz, args.min_gain, edges, lin)
        top = levers[:args.top]
        print(f"\n== {b['component']} {b['lo']:.2f}-{b['hi']:.2f} Hz: cost {b['cost']:.5f} "
              f"({100 * b['share']:.1f} % of the spectrum), relative residual {b['relative']:.3f}")
        print("   lever                                     wanted step   band gain  cost elsewhere  selectivity  verdict")
        for d in top:
            step = f"{d['step_hz']:+.3f} Hz" if d["step_hz"] is not None else f"{d['step']:+.4f}"
            print(f"   {d['label']:42s} {step:>11s}   {-100 * d['band_share']:+6.1f} %    {100 * d['elsewhere_share']:+7.1f} %"
                  f"      {d['selectivity']:5.2f}     {d['verdict']}")
        # the band's own parameters: couplings of its 13C (J(site, ...)) and its component's rate families
        site = b["component"].split("@")[-1].split(" ")[0] if "@" in b["component"] else None
        comp = labels.index(b["component"]) if b["component"] in labels else None

        def is_own(d):
            if d["kind"] == "coupling":
                return site is not None and d["name"].startswith(f"J({site},")
            return comp is not None and d["name"].startswith(f"c{comp}.")
        own = [d for d in levers if is_own(d) and d["verdict"] == "weak"]
        verdicts = {d["verdict"] for d in top}
        if verdicts == {"weak"}:
            summary = "no free parameter acts on this misfit: the model lacks something in this band"
        elif "free knob" in verdicts:
            summary = "free knobs: " + ", ".join(d["label"] for d in top if d["verdict"] == "free knob")
        else:
            summary = ("conflicts: the band wants " + ", ".join(
                f"{d['label']} {d['step_hz']:+.2f} Hz" if d["step_hz"] is not None else d["label"]
                for d in top if d["verdict"] == "conflict")
                + "; the rest of the spectrum resists (the same parameter serves other lines)")
        if own:
            exhausted = ("own parameters already at the band optimum: "
                         + ", ".join(d["label"] for d in own[:6]) + (" ..." if len(own) > 6 else ""))
            print("   ->", exhausted)
        else:
            exhausted = ""
        print("   ->", summary)
        entry = {**b, "levers": top, "summary": summary, "own_exhausted": [d["label"] for d in own]}
        if args.local_fit > 0:
            from scipy.optimize import least_squares
            chosen = [d for d in levers if d["verdict"] != "weak"][:args.local_fit]
            if chosen:
                zi = np.array([d["index"] for d in chosen])
                lo_b, hi_b = prob.lower[zi], prob.upper[zi]

                def full(q):
                    zz = z.copy()
                    zz[zi] = q
                    return zz

                sol = least_squares(lambda q: joint.residual(full(q)), np.clip(z[zi], lo_b + 1e-9, hi_b - 1e-9),
                                    jac=lambda q: joint.jacobian(full(q))[:, zi], bounds=(lo_b, hi_b),
                                    x_scale="jac", max_nfev=40)
                z_new = full(sol.x)
                pr = f.predict(joint.spectrum_vector(z_new, s)).residual
                rows_s = band_rows(f, b["lo"], b["hi"], len(pr))
                band_after = float(pr[rows_s] @ pr[rows_s])
                obj_after = float(np.sum(joint.residual(z_new) ** 2))
                print(f"   whole-objective refit of {', '.join(d['label'] for d in chosen)}: band {b['cost']:.5f} -> "
                      f"{band_after:.5f} ({100 * (band_after / b['cost'] - 1):+.1f} %), whole objective "
                      f"{objective:.5f} -> {obj_after:.5f} ({100 * (obj_after / objective - 1):+.2f} %)")
                entry["local_fit"] = {"parameters": [d["label"] for d in chosen], "band_cost_after": band_after,
                                      "objective_after": obj_after}
                entry["_z_after"] = z_new
        report["diagnosed"].append(entry)

    if args.save_start and report["diagnosed"] and "_z_after" in report["diagnosed"][0]:
        z_new = report["diagnosed"][0]["_z_after"]
        start = json.loads(json.dumps(fit))
        for k, n in enumerate(joint.coupling):
            start["couplings"].setdefault(prob.key_of.get(n, n), {})["J_at_x"] = joint.coupling_values(z_new, k).tolist()
        sid = prob.series[s]["id"]
        x = joint.spectrum_vector(z_new, s)
        start.setdefault("spectrum_parameters", {}).setdefault(sid, {})
        for n in list(joint.shared) + list(joint.local):
            start["spectrum_parameters"][sid][n] = float(x[joint.col[n]])
        Path(args.save_start).write_text(json.dumps(start, indent=1))
        print(f"start written: {args.save_start}")
    if args.figure:
        figure(args.figure, f, joint, s, z, report, labels)
    if args.json:
        clean = json.loads(json.dumps(report, default=lambda o: None))
        Path(args.json).write_text(json.dumps(clean, indent=1))


def figure(path, f, joint, s, z, report, labels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figure_tools import INK, MUTED, PALETTE, SIMULATION
    diag = report["diagnosed"]
    n = len(diag)
    fig = plt.figure(figsize=(13, 1.6 + 3.4 * n))
    gs = fig.add_gridspec(n + 1, 2, height_ratios=[0.9] + [2.0] * n, width_ratios=[1.35, 1.0], hspace=0.75,
                          wspace=0.32)
    ax = fig.add_subplot(gs[0, :])
    bands = report["bands"]
    comps = [c for c in labels if c in {b["component"] for b in bands}] + \
        [c for c in dict.fromkeys(b["component"] for b in bands) if c not in labels]
    color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(comps)}
    for b in bands:
        ax.barh(0, b["hi"] - b["lo"], left=b["lo"], height=0.6, color=color[b["component"]],
                alpha=0.25 + 0.75 * b["share"] / max(bb["share"] for bb in bands), edgecolor="none")
        if b["share"] > 0.05:
            ax.text(0.5 * (b["lo"] + b["hi"]), 0.45, f"{100 * b['share']:.0f} %", ha="center", va="bottom",
                    fontsize=7, color=INK)
    for c in comps:
        ax.barh(0, 0, color=color[c], label=c)
    ax.set_yticks([])
    ax.set_ylim(-0.5, 1.0)
    ax.set_xlabel("Frequency [Hz]", fontsize=8)
    ax.legend(fontsize=7, ncol=len(comps), frameon=False, loc="upper right", bbox_to_anchor=(1, 1.45))
    ax.set_title("share of the spectrum cost per isotopologue band (darker = worse)", fontsize=9, loc="left")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=7)
    nf = len(f.f)
    cm = lambda v: np.asarray(v)[:nf] + 1j * np.asarray(v)[nf:] if len(v) == 2 * nf else np.asarray(v)
    y = cm(f.y)
    m0 = cm(f.predict(joint.spectrum_vector(z, s)).model)
    for i, d in enumerate(diag):
        a = fig.add_subplot(gs[i + 1, 0])
        w = (f.f >= d["lo"] - 2) & (f.f <= d["hi"] + 2)
        a.plot(f.f[w], y.real[w], color=INK, lw=0.8, label="data")
        a.plot(f.f[w], m0.real[w], color=SIMULATION, lw=1.0, label="fit")
        if "_z_after" in d:
            m1 = cm(f.predict(joint.spectrum_vector(d["_z_after"], s)).model)
            a.plot(f.f[w], m1.real[w], color=PALETTE[3], lw=1.0, ls="--", label="refit of the levers")
        a.axvspan(d["lo"], d["hi"], color=color.get(d["component"], MUTED), alpha=0.08, lw=0)
        a.set_title(f"{d['component']} band {d['lo']:.1f}-{d['hi']:.1f} Hz: {100 * d['share']:.0f} % of the "
                    f"spectrum cost, relative residual {d['relative']:.2f}", fontsize=9, loc="left")
        a.legend(fontsize=7, frameon=False)
        a.tick_params(labelsize=7)
        a.set_xlabel("Frequency [Hz]", fontsize=8)
        b = fig.add_subplot(gs[i + 1, 1])
        lev = d["levers"][::-1]
        yy = np.arange(len(lev))
        gain = [100 * l["band_share"] for l in lev]
        cost = [100 * l["elsewhere_share"] for l in lev]
        cap = max(max(gain), 3.0) * 3.0                 # cost elsewhere is often far larger: clip and label it
        b.barh(yy + 0.18, gain, color=PALETTE[1], height=0.34, label="band cost removed")
        b.barh(yy - 0.18, [min(c, cap) for c in cost], color="#b03030", alpha=0.7, height=0.34,
               label="cost added elsewhere")
        for k, l, c in zip(yy, lev, cost):
            if c > cap:
                b.text(cap, k - 0.18, f" {c:.0f} %", va="center", fontsize=6, color="#b03030")
            step = f"{l['step_hz']:+.2f} Hz" if l["step_hz"] is not None else f"{l['step']:+.3g}"
            b.text(cap * 1.32, k, f"{step}  {l['verdict']}", va="center", fontsize=6.5,
                   color={"conflict": "#b03030", "free knob": PALETTE[1]}.get(l["verdict"], MUTED))
        b.set_yticks(yy)
        b.set_yticklabels([l["label"] for l in lev], fontsize=7)
        b.set_xlabel("share of the band cost, for the band's wanted step [%]", fontsize=8)
        b.tick_params(labelsize=7)
        b.set_xlim(min(0, min(cost)) * 1.1, cap * 2.1)
        b.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.16), ncol=2)
        if d.get("own_exhausted"):
            b.text(0.0, -0.33, "own parameters already at the band optimum: " + ", ".join(d["own_exhausted"][:5]),
                   transform=b.transAxes, fontsize=6.5, color=MUTED)
        b.set_title(d["summary"] if len(d["summary"]) < 90 else d["summary"][:87] + "...", fontsize=7.5, loc="left")
        lf = d.get("local_fit")
        if lf:
            a.text(0.01, 0.02, f"refit of the levers on the whole objective: band "
                   f"{100 * (lf['band_cost_after'] / d['cost'] - 1):+.0f} %, whole "
                   f"{100 * (lf['objective_after'] / report['objective'] - 1):+.1f} %", transform=a.transAxes,
                   fontsize=7, color=MUTED)
        for sp in ("right", "top"):
            b.spines[sp].set_visible(False)
            a.spines[sp].set_visible(False)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"figure: {path}")


if __name__ == "__main__":
    main()
