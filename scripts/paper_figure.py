"""Publication-style figure of a fit_joint_series result: (a) the whole spectrum, experiment above simulation;
(b) the line bands, experiment and simulation overlaid, with isotopologue bars and structure insets.

    python scripts/paper_figure.py <the fit_joint_series options of the fit> --fit RUN/fit.json \\
        --fid DATA/average_fid.npy [--display-window 0.1] [--wide 5,300] [--segments "104,158;222,276"] \\
        [--gains 1,2.5] [--insets insets.json] [--colors '{"13C@C1": "#2f8f5b"}'] [--title NAME] [--figure OUT.png]

Display processing (not the fit's): the FID from --display-crop s to its end, exponential window
--display-window 1/s, zero fill --display-zero-fill. Sampling rate: --sampling-rate, else the series entry's
record (sampling_rate_hz), else 4000 Hz. Zero-order phase: --phase0-deg, else the series entry's phasing, else the
instrument calibration (configs/confirmed_samples.json); the first-order phase from this FID's switching edge. The fitted model (couplings, decay rates, delay of --fit) is rendered through
this same processing; its gains are solved on the fit's ranges and kept for the whole range; the residual
zero-order phase of the gains (modulo 180 deg) is removed from data and model alike.

Display baseline (zulf_processing.anchor_spline_baseline, the same steps for data and model): (1) a cubic spline
through anchor points more than --protect Hz from every model line above 20 Hz and away from the power-line
harmonics and from narrow data peaks near the lines, knots every --knots Hz, protected runs longer than two knot
spacings bridged by a straight line (the noise-free model without outlier rejection); (2) within --lift Hz of the
model lines only, an AsLS 1.5 Hz baseline, which lifts the negative crop wings between close lines. The rolling
baseline (period about 1 / (crop + delay)) and those wings are processing effects of the record start, not
signal. Power-line harmonics are kept and marked. Nothing here changes the fit.

--insets: JSON list of {"component": label, "smiles": "...", "atom": heavy-atom index of the 13C,
"segment": 0, "rect": [f0_hz, y0, width_hz, height]} (RDKit drawings, scripts/figure_tools.py; skipped without
RDKit; "bond_note": a coupling key, e.g. "J(C2,HC2)", writes that J, with sigma from --uncertainty, on the 13C-H
bond). Writes OUT.png (and --formats pdf,svg) and OUT.caption.txt. --manual labels the model "manual parameters
(not a fit)" (a parameter file written by ZULF Studio from moved sliders).
"""
import json
import sys
from io import BytesIO
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import fit_joint_series as fj                                            # noqa: E402
from figure_tools import GRID, INK, MUTED, PALETTE, SIMULATION, isotopologue_png  # noqa: E402
from j_tuner import load_fit                                             # noqa: E402

MAINS_HZ = 60.06


def display_spectrum(fid_path, crop_s, window, zero_fill, fs=4000.0, phase0_rad=None):
    from zulf_core.render.phasing import correction_phasor, reference_delay_s
    from zulf_processing import plan_for_dataset, process_dataset
    from zulf_processing.diagnostics import switching_edge
    cal = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]["phase_calibration"]
    fid = np.load(fid_path).astype(float)
    edge = switching_edge(fid, fs)["edge_time_s"] + cal["delay_offset_s"]
    start = int(round(crop_s * fs))
    plan = plan_for_dataset(len(fid), fs, None, {"start_sample": start, "stop_sample": len(fid), "zero_fill": zero_fill,
                                                "apodization_rate_per_s": window, "ranges": [[1.0, fs / 2 - 1.0]]})
    ds = process_dataset(fid, fs, plan=plan, phase_criterion=None)
    acq = plan.acquisition()
    phi0 = float(np.radians(cal["phase0_deg"])) if phase0_rad is None else float(phase0_rad)
    delay = float(-(edge + acq.time_origin_s))
    spec = ds.spectrum * correction_phasor(ds.frequencies_hz, phi0, delay + reference_delay_s(acq))
    return ds.frequencies_hz, spec, acq, {"phase0_rad": phi0, "delay_s": delay}, (len(fid) - start) / fs


def render_model(prob, z, f, spec, acq, phasing, wide, fit_ranges):
    """Model on the display grid: gains solved on fit_ranges, then kept for `wide`."""
    from zulf_core.solver import ObservedSpectrum
    from zulf_core.solver.forward import MixtureForward
    from zulf_core.solver.refine import _signal_kwargs
    j = prob.joint
    f0 = j.forwards[0]

    def forward(ranges):
        obs = ObservedSpectrum.from_spectrum(f, spec, ranges, record=acq.to_dict(), phasing=phasing, real_only=False)
        return obs, MixtureForward(j.params[0], obs, f0.protocol, j.settings.gain_model, j.settings.background_order,
                                   j.settings.band_weighting, **_signal_kwargs(j.settings))
    x = j.spectrum_vector(z, 0)
    _, fw_fit = forward(fit_ranges)
    gains = fw_fit.predict(x).gains
    obs_w, fw_w = forward([wide])
    pred = fw_w.predict(x, fixed_gains=gains)
    values = j.params[0].values(x)
    lines = []
    for c, system in enumerate(j.params[0].systems(values)):
        tl = fw_w.transitions(values, c, system)
        a = np.abs(np.asarray(tl.amplitudes))
        if len(a):
            lines.append((c, np.asarray(tl.frequencies_hz), a / a.max()))
    return np.asarray(obs_w.frequencies_hz), np.asarray(obs_w.values), np.asarray(pred.model), gains, lines


def display_baseline(f, y, m, lines, protect_hz, knots_hz, lift_hz):
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    from zulf_processing import anchor_spline_baseline, asls_baseline, line_mask
    model_lines = [v for _, fr, a in lines for v in fr[a >= 0.03] if v > 20.0]
    step = float(np.median(np.diff(f)))
    quiet = (f > f.max() - 10) & (f <= f.max())
    noise = 1.4826 * np.median(np.abs(y[quiet] - np.median(y[quiet])))
    peaks, _ = find_peaks(np.abs(y), prominence=6 * noise, width=(None, 1.5 / step))
    mains = [MAINS_HZ * k for k in range(1, int(f.max() / MAINS_HZ) + 1)]
    protect = (line_mask(f, model_lines, protect_hz) | line_mask(f, mains, 0.6)
               | (line_mask(f, f[peaks], 0.8) & line_mask(f, model_lines, 5.0)))
    yc = y - anchor_spline_baseline(y, f, protect, knot_spacing_hz=knots_hz)
    mc = m - anchor_spline_baseline(m, f, protect, knot_spacing_hz=knots_hz, k_sigma=np.inf)
    taper = np.clip(gaussian_filter1d(line_mask(f, model_lines, lift_hz).astype(float), 1.0 / step), 0, 1)
    yc = yc - taper * asls_baseline(yc, f, smooth_hz=1.5, p=0.01)
    mc = mc - taper * asls_baseline(mc, f, smooth_hz=1.5, p=0.01)
    return yc, mc, mains


def clusters(fr, amp, lo, hi, threshold=0.12, gap=2.6, pad=0.6):
    v = np.sort(fr[(amp >= threshold) & (fr >= lo) & (fr <= hi)])
    out, start = [], None
    for i, x in enumerate(v):
        start = x if start is None else start
        if i == len(v) - 1 or v[i + 1] - x > gap:
            out.append((start - pad, x + pad))
            start = None
    return out


def main():
    ap = fj.make_parser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--fid", required=True)
    ap.add_argument("--display-crop", type=float, default=0.1)
    ap.add_argument("--display-window", type=float, default=0.1)
    ap.add_argument("--display-zero-fill", type=int, default=4)
    ap.add_argument("--wide", default="5,300")
    ap.add_argument("--segments", default="", help="lo,hi;lo,hi detail panels (default: the fit's ranges hull)")
    ap.add_argument("--gains", default="", help="per detail panel display factor, e.g. 1,2.5")
    ap.add_argument("--protect", type=float, default=1.0)
    ap.add_argument("--knots", type=float, default=1.5)
    ap.add_argument("--lift", type=float, default=2.5)
    ap.add_argument("--stacked-detail", action="store_true", help="detail panels stacked instead of overlaid")
    ap.add_argument("--insets", default="")
    ap.add_argument("--colors", default="{}", help="JSON {component label: colour}")
    ap.add_argument("--uncertainty", default="", help="coupling_diagram.py OUT.json: sigma for the bond notes")
    ap.add_argument("--title", default="")
    ap.add_argument("--figure", default="")
    ap.add_argument("--sampling-rate", type=float, default=0.0, help="Hz (default: the series entry's record, else 4000)")
    ap.add_argument("--phase0-deg", type=float, default=None, help="display zero-order phase (default: series phasing)")
    ap.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--manual", action="store_true", help="the parameters are not a fit (label and caption say so)")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    prob = fj.build_problem(args)
    fit = json.load(open(args.fit))
    load_fit(prob, fit)
    z = prob.z0
    j = prob.joint
    labels = list(j.model.component_labels)
    colors = {lab: PALETTE[i % len(PALETTE)] for i, lab in enumerate(labels)}
    colors.update(json.loads(args.colors))
    wide = tuple(float(v) for v in args.wide.split(","))
    fit_ranges = [tuple(r) for r in prob.series[0].get("ranges", [(prob.lo, prob.hi)])]
    entry = prob.series[0]
    fs = args.sampling_rate or float((entry.get("record") or {}).get("sampling_rate_hz", 4000.0))
    phase0 = (np.radians(args.phase0_deg) if args.phase0_deg is not None
              else (entry.get("phasing") or {}).get("phase0_rad"))
    f, spec, acq, phasing, record_s = display_spectrum(args.fid, args.display_crop, args.display_window,
                                                       args.display_zero_fill, fs=fs, phase0_rad=phase0)
    sim_label = "Simulation (manual parameters, not a fit)" if args.manual else "Simulation (fit)"
    f, y, m, gains, lines = render_model(prob, z, f, spec, acq, phasing, wide, fit_ranges)
    phi = float(np.angle(gains[0]))
    phi -= np.pi * round(phi / np.pi)
    y, m = (y * np.exp(-1j * phi)).real, (m * np.exp(-1j * phi)).real
    yc, mc, mains = display_baseline(f, y, m, lines, args.protect, args.knots, args.lift)
    segs = ([tuple(float(v) for v in s.split(",")) for s in args.segments.split(";")] if args.segments
            else [(min(r[0] for r in fit_ranges), max(r[1] for r in fit_ranges))])
    seg_gain = [float(v) for v in args.gains.split(",")] if args.gains else [1.0] * len(segs)

    fig = plt.figure(figsize=(15, 8.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25], hspace=0.32, left=0.055, right=0.99, top=0.96, bottom=0.08)
    # (a) whole range, stacked
    ax0 = fig.add_subplot(outer[0])
    sel = (f >= wide[0]) & (f <= wide[1])
    inband = np.zeros(len(f), bool)
    for lo, hi in segs:
        inband |= (f >= lo) & (f <= hi)
    sc0 = max(float(np.max(yc[sel & inband])), 1e-30)
    off0 = -1.15
    ax0.plot(f[sel], yc[sel] / sc0, color=INK, lw=0.6)
    ax0.plot(f[sel], off0 + mc[sel] / sc0, color=SIMULATION, lw=0.6)
    for v in mains:
        if wide[0] < v < wide[1]:
            ax0.plot([v], [1.1], marker="v", ms=5, color=MUTED, ls="none")
    ax0.set_xlim(*wide)
    ax0.set_ylim(off0 - 0.15, 1.2)
    for yy in (0, off0):
        ax0.axhline(yy, color=GRID, lw=0.6, zorder=0)
    ax0.spines[["top", "right"]].set_visible(False)
    ax0.set_yticks([0, 0.5, 1])
    ax0.spines["left"].set_bounds(0, 1)
    ax0.tick_params(colors=MUTED, labelcolor=INK, labelsize=10)
    ax0.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
    ax0.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))
    ax0.set_ylabel("Signal [a.u.]", color=INK, fontsize=11)
    ax0.set_xlabel("Frequency [Hz]", color=INK, fontsize=10.5)
    for txt, yy, fc in (("Experimental", 0.55, INK), (sim_label, off0 + 0.55, SIMULATION)):
        ax0.text(wide[0] + 1, yy, txt, color="white", fontsize=10, fontweight="bold", va="center",
                 bbox=dict(boxstyle="round,pad=0.3", fc=fc, ec="none"))
    ax0.text(wide[1], 1.1, f"v  {MAINS_HZ:g} Hz power-line harmonics", ha="right", va="center", fontsize=8.5, color=MUTED)
    for a, b in segs:
        ax0.add_patch(plt.Rectangle((a, off0 - 0.12), b - a, 1.2 - off0, fc="none", ec=GRID, lw=0.8, ls="--"))
    sp0 = dict(zip(j.local, z[j.nt:j.nt + j.nl]))
    sp0.update(dict(zip(j.shared, z[j.ntheta:j.ntheta + len(j.shared)])))
    bt0, bz0 = sp0.get("field_transverse_ut"), sp0.get("field_z_ut")
    field_label = ("zero field" if bt0 is None and bz0 is None else
                   f"B transverse {1e3 * (bt0 or 0):.0f} nT, B z {1e3 * (bz0 or 0):.0f} nT, "
                   f"|B| {1e3 * np.hypot(bt0 or 0, bz0 or 0):.0f} nT")
    ax0.text(0.0, 1.03, f"a  {args.title + ', ' if args.title else ''}{field_label}; whole spectrum, "
             f"{wide[0]:g}-{wide[1]:g} Hz",
             transform=ax0.transAxes, fontsize=11, fontweight="bold", color=INK)
    # (b) detail panels, overlaid
    inner = outer[1].subgridspec(1, len(segs), width_ratios=[b - a for a, b in segs], wspace=0.035)
    axes = [fig.add_subplot(inner[0, i]) for i in range(len(segs))]
    scale = sc0                       # one scale for every detail panel: the largest data value in any of them
    off = -1.18 if args.stacked_detail else 0.0
    for i, (ax, (a, b)) in enumerate(zip(axes, segs)):
        sl = (f >= a) & (f <= b)
        g = seg_gain[i] / scale
        ax.plot(f[sl], off + mc[sl] * g, color=SIMULATION, lw=0.9 if args.stacked_detail else 1.6,
                alpha=1 if args.stacked_detail else 0.8)
        ax.plot(f[sl], yc[sl] * g, color=INK, lw=0.8)
        for v in mains:
            if a < v < b:
                ax.plot([v], [1.02], marker="v", ms=5, color=MUTED, ls="none")
        ax.set_xlim(a, b)
        ax.set_ylim(off - 0.08 if args.stacked_detail else -0.14, 1.24)
        ax.axhline(0, color=GRID, lw=0.6, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=10)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(2))
        if i:
            ax.spines["left"].set_visible(False)
            ax.tick_params(left=False, labelleft=False)
        else:
            ax.set_yticks([0, 0.5, 1])
            ax.spines["left"].set_bounds(0, 1)
            ax.set_ylabel("Signal [a.u.]", color=INK, fontsize=11)
            ax.text(0.0, 1.03, "b  the line bands", transform=ax.transAxes, fontsize=11, fontweight="bold", color=INK)
        if seg_gain[i] != 1.0:
            ax.text(b - 0.5, 0.97, f"x{seg_gain[i]:g}", ha="right", fontsize=10, color=MUTED)
        for c, fr, amp in lines:
            for lo, hi in clusters(fr, amp, a, b):
                ax.add_patch(FancyBboxPatch((lo, 1.13), hi - lo, 0.05, boxstyle="round,pad=0,rounding_size=0.25",
                                            mutation_aspect=0.12, fc=colors[labels[c]], ec="none", alpha=0.45))
    for k in range(len(axes) - 1):                      # broken-axis marks
        for ax, side in ((axes[k], 1.0), (axes[k + 1], 0.0)):
            for dx in (-0.006, 0.006):
                ax.plot((side + dx - 0.006, side + dx + 0.006), (-0.025, 0.025), transform=ax.transAxes, color=INK,
                        clip_on=False, lw=1.0)
    fig.text(0.52, 0.035, "Frequency [Hz]", ha="center", color=INK, fontsize=10.5)
    if not args.stacked_detail:
        axes[-1].legend([Line2D([], [], color=INK, lw=1.2), Line2D([], [], color=SIMULATION, lw=2, alpha=0.8)],
                        ["Experimental", sim_label], loc="upper right", bbox_to_anchor=(1.0, 0.74),
                        frameon=False, fontsize=10)
    if args.insets:
        import matplotlib.image as mpimg
        sig = {}
        if args.uncertainty:
            sig = {r["key"]: r["sigma"] for r in json.load(open(args.uncertainty))}
        for ins in json.load(open(args.insets)):
            note = None
            if ins.get("bond_note") in fit["couplings"]:          # e.g. "J(C2,HC2)": its 1J on the C-H bond
                key = ins["bond_note"]
                jv = fit["couplings"][key]["J_at_x"][0]
                note = (f"{jv:.2f} \u00b1 {sig[key]:.2f} Hz" if key in sig else f"{jv:.2f} Hz")
            out_png = isotopologue_png(ins["smiles"], ins["atom"], colors[ins["component"]], return_coords=True,
                                       mark_bond=bool(note))
            if out_png is None:
                print("RDKit not installed: structure insets skipped")
                break
            png, pts = out_png
            seg = axes[ins.get("segment", 0)]
            a_ = seg.inset_axes(ins["rect"], transform=seg.transData)
            a_.imshow(mpimg.imread(BytesIO(png), format="png"))
            if note:
                # 1J under the drawing, joined to the marked 13C-H bond
                mid = 0.5 * (np.array(pts["c"]) + np.array(pts["h"]))
                a_.annotate("$^1J$ = " + note, xy=mid, xytext=(0.5 * pts["size"][0], pts["size"][1] * 1.04),
                            textcoords="data", ha="center", va="top", fontsize=9.5, color=colors[ins["component"]],
                            fontweight="bold", annotation_clip=False,
                            arrowprops=dict(arrowstyle="-", color=colors[ins["component"]], lw=0.9, alpha=0.8,
                                            shrinkA=2, shrinkB=4))
            a_.axis("off")
    J = {prob.key_of.get(n, n): j.coupling_values(z, k)[0] for k, n in enumerate(j.coupling)}
    sp = dict(zip(j.local, z[j.nt:j.nt + j.nl]))
    sp.update(dict(zip(j.shared, z[j.ntheta:j.ntheta + len(j.shared)])))
    bt, bz = sp.get("field_transverse_ut"), sp.get("field_z_ut")
    field = ("at zero field" if bt is None and bz is None else
             f"in a static field (B transverse {1e3 * (bt or 0):.1f} nT, B z {1e3 * (bz or 0):.1f} nT, "
             f"{'set by hand' if args.manual else 'fitted'}, D47)")
    caption = (f"{args.title or prob.series[0]['id']} {field}, sampling {fs:g} Hz. "
               + ("Manual parameters (not a fit): " if args.manual else "Fit: ")
               + ", ".join(f"{k} {v:.2f}" for k, v in J.items()) + " Hz. "
               f"Display processing: record {args.display_crop:g}-{args.display_crop + record_s:.1f} s, window "
               f"{args.display_window:g} 1/s; the model rendered through the same processing (gains of the fit "
               f"ranges). Display baseline corrected, the same for data and model: cubic spline through anchor "
               f"points more than {args.protect:g} Hz from the model lines (knots every {args.knots:g} Hz), then "
               f"AsLS 1.5 Hz within {args.lift:g} Hz of the lines. Power-line harmonics kept and marked. Bars: lines "
               f"of each 13C isotopologue.")
    out = Path(args.figure or Path(args.fit).parent / "paper_figure.png")
    for fmt in [v.strip() for v in args.formats.split(",") if v.strip()]:
        fig.savefig(out.with_suffix("." + fmt), dpi=args.dpi)
    out.with_suffix(".caption.txt").write_text(caption + "\n")
    print(out)
    print(caption)


if __name__ == "__main__":
    main()
