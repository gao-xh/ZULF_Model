"""Line table of a structure: every line in a band, its isotopologue, amplitude, df/dJ for every coupling and the
exact split of its frequency into coupling contributions (zulf_core.physics.lines).

    python scripts/line_table.py --structure '{"motif": "(CH3)2CH-NH2", "one_bond": {"C1": 133, "C2": 125}}' \\
        [--exchange fast] [--couplings '{"J(HC2,HC3)": 0.2}'] [--fit runs/.../fit.json] [--x 1.0] \\
        [--band 115,140] [--min-amplitude 0.05] [--shift '{"J(C2,HC1)": -0.5}'] [--json out.json]

--fit takes the couplings of a fit_joint_series result (at the node nearest --x); --couplings overrides on top.
--shift prints the first-order positions after changing the named couplings by the given amounts (Hz).
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


def lines_for(spec, exchange="slow", couplings=None, band=None, min_amplitude=0.0):
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
                                                          min_relative_amplitude=min_amplitude)))
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
    table = lines_for(json.loads(args.structure), args.exchange, couplings, band, args.min_amplitude)
    record = []
    for label, lines in table:
        print(f"== {label}: {len(lines)} lines")
        for line in lines:
            moved = f"  -> {line.shifted(shift):.3f} Hz" if shift else ""
            top = sorted(line.contribution_hz.items(), key=lambda kv: -abs(kv[1]))
            print(f"  {line.frequency_hz:8.3f} Hz  amp {line.relative_amplitude:.2f}{moved}")
            print("      df/dJ: " + ", ".join(f"{n} {v:+.3f}" for n, v in sorted(line.df_dj.items(), key=lambda kv: -abs(kv[1]))))
            print("      Hz from: " + ", ".join(f"{n} {v:+.3f}" for n, v in top))
            record.append({"component": label, "frequency_hz": line.frequency_hz,
                           "relative_amplitude": line.relative_amplitude, "df_dJ": line.df_dj,
                           "contribution_hz": line.contribution_hz,
                           **({"shifted_hz": line.shifted(shift)} if shift else {})})
    if args.json:
        json.dump(record, open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
