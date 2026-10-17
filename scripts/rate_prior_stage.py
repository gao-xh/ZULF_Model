"""Stage 4 (D69): decay-rate and line-position priors from the line-shape probe, then a full refit.

    python scripts/rate_prior_stage.py RUN_FOLDER OUT_FOLDER [--windows auto | lo,hi;lo,hi] [--peaks 3]
        [--half-width 2.0] [--groups 8] [--max-shift 0.6] [--rate-up 1.5] [--sigma 0.3] [--weight 0.1]
        [--line-sigma 0.1] [--line-weight 0.1] [--max-nfev 40] [--probe-only]

From a finished fit (normally stage 3; stages 1-3 are not changed):
  1. windows: `auto` takes the `--peaks` largest peaks of the plain residual inside the data's line regions
     with a model line nearby (fit_joint_series plain_residual_peaks / at_a_line), +-`--half-width` Hz around
     each, overlapping windows merged; or the windows given;
  2. in every window the line-shape probe (scripts/line_shape_probe.py) frees the positions of the strongest line
     groups and the decay-rate families there, so a rate is not widened to cover a misplaced line;
  3. priors in a refit from RUN_FOLDER with its recorded command (same family edges; one start; `--max-nfev`
     evaluations; no residual refine; one BLAS thread), written to OUT_FOLDER/rate_prior.json:
     - rates: the probed rates (geometric mean where windows share a family; a rate below 0.2 1/s or at the 1.5 x
       ceiling is a bound, not a measurement, and is left out) as Gaussian priors on the log rates (`--sigma`
       log units; `--weight`: all rate priors 1 sigma off together cost that fraction of the start objective);
     - line positions (owner, 2026-10-17): every probed group moved by more than 0.02 Hz (not at the shift
       bound) pulls its line toward the probed position, linearised in the couplings by the group's df/dJ
       (`--line-sigma` Hz, `--line-weight` likewise); the couplings move the lines as far as the structure allows;
  4. acceptance: the refit is kept only when the data objective of RUN_FOLDER's problem (no prior rows) falls;
     OUT_FOLDER/stage4.json gets objective_before / objective_after / accepted (run_batch then reports stage 4 only
     when accepted; otherwise the input stays the result).
The refit is a conditional numerical result of the given structure. OUT_FOLDER/stage4.json records the windows,
the probe results, the prior table and the acceptance.
"""
import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def auto_windows(joint, z, peaks=3, half_width=2.0, sigma=4.0):
    """+-half_width Hz around the largest plain-residual peaks inside the line regions of spectrum 0, merged."""
    import fit_joint_series as fj
    # a model line (>= 10 % of its isotopologue's strongest) within 0.8 Hz: the probe moves model lines, so a
    # data feature without one (41.6 Hz in sec-butylamine, 86 Hz in DMEDA) gives it nothing to do
    found = [pk for pk in fj.plain_residual_peaks(joint, z, sigma) if pk[1] == 0
             and fj.at_a_line(joint, z, 0, pk[2], 0.1)][:peaks]
    spans = sorted((p[2] - half_width, p[2] + half_width) for p in found)
    merged = []
    for lo, hi in spans:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [tuple(w) for w in merged], [(round(p[2], 3), round(p[0], 1)) for p in found]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--windows", default="auto")
    ap.add_argument("--peaks", type=int, default=3)
    ap.add_argument("--half-width", type=float, default=2.0)
    ap.add_argument("--groups", type=int, default=8)
    ap.add_argument("--max-shift", type=float, default=0.6)
    ap.add_argument("--rate-up", type=float, default=1.5)
    ap.add_argument("--sigma", type=float, default=0.3)
    ap.add_argument("--weight", type=float, default=0.1)
    ap.add_argument("--line-sigma", type=float, default=0.1, help="width (Hz) of the line-position priors")
    ap.add_argument("--line-weight", type=float, default=0.1, help="line priors: fraction of the start objective "
                    "at 1 sigma; 0 = rate priors only")
    ap.add_argument("--max-nfev", type=int, default=40)
    ap.add_argument("--probe-only", action="store_true", help="write the priors, do not refit")
    args = ap.parse_args()
    import line_shape_probe as lsp
    from fit_monitor import SpectrumCache
    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    joint = SpectrumCache().get(str(run)).joint
    z = np.array(json.load(open(run / "fit.json"))["z_final"], float)
    if args.windows == "auto":
        windows, peaks = auto_windows(joint, z, args.peaks, args.half_width)
    else:
        windows = [tuple(float(v) for v in w.split(",")) for w in args.windows.split(";")]
        peaks = []
    print(f"stage 4: windows {[(round(a, 2), round(b, 2)) for a, b in windows]} (peaks {peaks})", flush=True)
    probes, logs = [], {}
    for w in windows:
        res = lsp.probe(joint, z, w, groups=args.groups, max_shift=args.max_shift, rate_up=args.rate_up)
        print(f"  probe {w[0]:.2f}-{w[1]:.2f} Hz: rms {res['rms_window'][0]:.4f} -> {res['rms_window'][1]:.4f}; "
              + ", ".join(f"{n} {a:.2f}->{b:.2f}" for n, (a, b) in res["rates"].items()), flush=True)
        lsp.figure(res, out / f"probe_{w[0]:.1f}-{w[1]:.1f}.png", title=run.name)
        probes.append({k: v for k, v in res.items() if not k.startswith("_")})
        for n, (a, b) in res["rates"].items():
            if b < 0.2 or b >= 0.98 * args.rate_up * a:      # at or near the probe's bounds: not a measurement
                continue
            logs.setdefault(n, []).append(np.log(b))
    rates = {n: float(np.exp(np.mean(v))) for n, v in logs.items()}
    # line positions (owner, 2026-10-17): every probed group moved by more than 0.02 Hz and not at the shift bound
    # becomes a line-position prior, linearised in the couplings by its df/dJ (fit_joint_series --line-prior-*)
    lines = [{"component": g["component"], "frequency_hz": g["frequency_hz"], "target_hz": g["new_frequency_hz"],
              "relative_amplitude": g["relative_amplitude"], "df_dJ": g["df_dJ"]}
             for pr in probes for g in pr["groups"]
             if 0.02 < abs(g["shift_hz"]) < 0.98 * args.max_shift and g["df_dJ"]]
    json.dump({"rates": rates, "lines": lines}, open(out / "rate_prior.json", "w"), indent=1)
    print(f"stage 4: {len(rates)} rate priors, {len(lines)} line-position priors", flush=True)
    json.dump({"run": str(run), "windows": windows, "peaks": peaks, "probes": probes, "rates": rates,
               "settings": vars(args)}, open(out / "stage4.json", "w"), indent=1)
    if args.probe_only or not rates:
        print("stage 4: no refit" + ("" if rates else " (no rate families in the windows)"))
        return
    argv = list(json.load(open(run / "monitor" / "status.json"))["argv"])

    def setopt(name, value):
        if name in argv:
            argv[argv.index(name) + 1] = value
        else:
            argv.extend([name, value])
    for name, value in (("--from-joint", str(run / "fit.json")), ("--out", str(out)), ("--starts", "1"),
                        ("--max-nfev", str(args.max_nfev)), ("--residual-refine", "0"), ("--workers", "1"),
                        ("--rate-prior", str(out / "rate_prior.json")), ("--rate-prior-sigma", str(args.sigma)),
                        ("--rate-prior-weight", str(args.weight)), ("--line-prior-sigma", str(args.line_sigma)),
                        ("--line-prior-weight", str(args.line_weight))):
        setopt(name, value)
    for flag in ("--sign-scan", "--shape-refine", "--peak-refine", "--one-bond-from-data"):
        while flag in argv:
            argv.remove(flag)
    setopt("--component-search", "off")
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    cmd = [sys.executable] + argv
    print("stage 4 refit: " + " ".join(shlex.quote(c) for c in cmd[:3]) + " ...", flush=True)
    code = subprocess.call(cmd, cwd=ROOT, env=env)
    if code or not (out / "fit.json").exists():
        sys.exit(code or 1)
    # acceptance (owner, 2026-10-17): the data objective of the input's problem (no prior rows) must fall
    before = float(np.sum(joint.residual(z) ** 2))
    z4 = np.array(json.load(open(out / "fit.json"))["z_final"], float)
    after = float(np.sum(joint.residual(z4) ** 2))
    accepted = after < before
    record = json.load(open(out / "stage4.json"))
    record.update(objective_before=before, objective_after=after, accepted=accepted)
    json.dump(record, open(out / "stage4.json", "w"), indent=1)
    print(f"stage 4: data objective {before:.5f} -> {after:.5f}: {'kept' if accepted else 'not kept (stage 3 stays)'}",
          flush=True)


if __name__ == "__main__":
    main()
