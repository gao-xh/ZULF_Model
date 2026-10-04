"""Line table of a structure: every line in a band, its isotopologue, amplitude, df/dJ for every coupling and the
exact split of its frequency into coupling contributions (zulf_core.physics.lines).

    python scripts/line_table.py --structure '{"motif": "(CH3)2CH-NH2", "one_bond": {"C1": 133, "C2": 125}}' \\
        [--exchange fast] [--couplings '{"J(HC2,HC3)": 0.2}'] [--fit runs/.../fit.json] [--x 1.0] \\
        [--band 115,140] [--min-amplitude 0.05] [--shift '{"J(C2,HC1)": -0.5}'] [--second-order] [--json out.json]

--fit takes the couplings of a fit_joint_series result (at the node nearest --x); --couplings overrides on top.
--shift prints the first-order positions after changing the named couplings by the given amounts (Hz).
--second-order adds d2f/dJ2 and the cross terms d2f/dJk dJl, the distance to the nearest other line and a trust
step (the coupling change at which the largest second-order term reaches 0.02 Hz); --shift then also gives the
second-order position. Near-degenerate lines (nearest line much closer than the trust step suggests) need an exact
recomputation, not a Taylor step.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from zulf_core.physics.lines import line_table              # noqa: E402
from zulf_hypothesis import build_model, exchange_variants  # noqa: E402


def lines_for(spec, exchange="slow", couplings=None, band=None, min_amplitude=0.0, second_order=False):
    """[(component label, [Line])] with couplings named by their structure keys (e.g. J(C2,HC1))."""
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    spec = dict(spec)
    spec.setdefault("compound", "structure")
    fragment = reg.structure_for(spec)
    if couplings:
        fragment = override_couplings(fragment, couplings)
    fragment = exchange_variants(fragment, exchange)[0]
    model = build_model(fragment)
    from zulf_core.solver import RefineSettings
    params = model.settings(RefineSettings(), fix_unspecified=False).parameterize(model.interpretation)
    key_of = {n: key for key, names in model.coupling_names.items() for n in names}
    values = params.values(params.vector())
    out = []
    for c, sys_ in enumerate(params.systems(values)):
        names = {}
        for n in params.coupling_names(c):
            names[tuple(params.parameters[n].detail)] = key_of.get(n, n)
        out.append((model.component_labels[c], line_table(sys_, names, band_hz=band,
                                                          min_relative_amplitude=min_amplitude,
                                                          second_order=second_order)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--exchange", default="slow", choices=["slow", "fast"])
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--fit", default="")
    ap.add_argument("--x", type=float, default=None)
    ap.add_argument("--band", default="")
    ap.add_argument("--min-amplitude", type=float, default=0.05)
    ap.add_argument("--shift", default="{}")
    ap.add_argument("--second-order", action="store_true")
    ap.add_argument("--json", default="")
    args = ap.parse_args()
    couplings = {}
    if args.fit:
        fit = json.load(open(args.fit))
        node = int(np.argmin(np.abs(np.asarray(fit["x"]) - (args.x if args.x is not None else fit["x"][0]))))
        couplings = {k: c["J_at_x"][node] for k, c in fit["couplings"].items()}
    couplings.update(json.loads(args.couplings))
    band = tuple(float(v) for v in args.band.split(",")) if args.band else None
    shift = json.loads(args.shift)
    table = lines_for(json.loads(args.structure), args.exchange, couplings, band, args.min_amplitude,
                      args.second_order)
    record = []
    for label, lines in table:
        print(f"== {label}: {len(lines)} lines")
        for line in lines:
            moved = f"  -> {line.shifted(shift):.3f} Hz (1st order)" if shift else ""
            if shift and args.second_order:
                moved += f", {line.shifted(shift, second_order=True):.3f} Hz (2nd order)"
            top = sorted(line.contribution_hz.items(), key=lambda kv: -abs(kv[1]))
            print(f"  {line.frequency_hz:8.3f} Hz  amp {line.relative_amplitude:.2f}{moved}")
            print("      df/dJ: " + ", ".join(f"{n} {v:+.3f}" for n, v in sorted(line.df_dj.items(), key=lambda kv: -abs(kv[1]))))
            print("      Hz from: " + ", ".join(f"{n} {v:+.3f}" for n, v in top))
            if line.hessian is not None:
                keys = list(line.df_dj)
                diag = sorted(((line.hessian[(k, k)], k) for k in keys), key=lambda t: -abs(t[0]))
                cross = sorted(((line.hessian[(a, b)], a, b) for i, a in enumerate(keys) for b in keys[i + 1:]),
                               key=lambda t: -abs(t[0]))[:3]
                print("      d2f/dJ2: " + ", ".join(f"{k} {v:+.3f}" for v, k in diag if abs(v) >= 1e-3))
                print("      cross: " + ", ".join(f"{a} x {b} {v:+.3f}" for v, a, b in cross if abs(v) >= 1e-3))
                warn = ("  NEAR-DEGENERATE: Taylor steps unreliable, recompute exactly"
                        if line.near_degenerate() else "")
                print(f"      nearest other line {line.nearest_line_hz:.3f} Hz, trust step {line.trust_step_hz:.2f} Hz{warn}")
            record.append({"component": label, "frequency_hz": line.frequency_hz,
                           "relative_amplitude": line.relative_amplitude, "df_dJ": line.df_dj,
                           "contribution_hz": line.contribution_hz,
                           "nearest_line_hz": line.nearest_line_hz,
                           **({"hessian": {f"{a}|{b}": v for (a, b), v in line.hessian.items()},
                               "trust_step_hz": line.trust_step_hz} if line.hessian is not None else {}),
                           **({"shifted_hz": line.shifted(shift)} if shift else {})})
    if args.json:
        json.dump(record, open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
