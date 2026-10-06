"""Best vector of a running (or killed) fit_joint_series run as a fit.json-like start file.

    python scripts/snapshot_fit.py runs/processed/NAME OUT.json

Reads NAME/monitor (fit_monitor records): every start's best vector ("z" with its own cost "best"), rebuilds the
problem from the recorded command line and rescores every candidate on the plain objective (no smoothing, hard
missing-peak rows): the recorded costs of different stages (smoothing, model-line weights) are not comparable.
Writes couplings, spectrum parameters and the rate-family edges, readable with --from-joint and
j_tuner.load_fit. Use it before a cloud container is reclaimed, or to restart a run from where it stopped.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_joint_series import make_parser, build_problem  # noqa: E402
from run_problem import problem_argv  # noqa: E402


def problem_from_status(status):
    return build_problem(make_parser().parse_args(problem_argv(list(status["argv"][1:])) + ["--monitor", "off"]))


def candidates(monitor_dir):
    """(file name, recorded cost, z): the last best vector of every record file."""
    out = []
    for path in sorted(Path(monitor_dir).glob("*.jsonl")):
        last = None
        for line in open(path):
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:            # a partly written last line of a killed run
                continue
            if "z" in rec and rec.get("best") is not None:
                last = (float(rec["best"]), rec["z"])
        if last is not None:
            out.append((path.stem, last[0], np.asarray(last[1], float)))
    return out


def snapshot(run_dir, out_path):
    run_dir = Path(run_dir)
    with open(run_dir / "monitor" / "status.json") as fh:
        status = json.load(fh)
    prob = problem_from_status(status)
    j = prob.joint
    j.set_smoothing(0.0)
    j.peak_smooth = 0.0
    found = [c for c in candidates(run_dir / "monitor") if len(c[2]) == len(prob.z0)]
    if not found:
        raise SystemExit(f"no best vectors of this problem's size in {run_dir / 'monitor'}")
    scored = sorted(((float(np.sum(j.residual(z) ** 2)), name, rec, z) for name, rec, z in found),
                    key=lambda t: t[0])
    score, name, recorded, z = scored[0]
    fit = {"x": j.nodes.tolist(), "scores": [score], "couplings": {}, "spectrum_parameters": {},
           "family_edges_hz": [float(v) for v in j.params[0].policy.family_edges_hz],
           "snapshot": {"run": str(run_dir), "record": name, "recorded_cost": recorded,
                        "candidates": [{"record": n, "objective": s, "recorded_cost": r} for s, n, r, _ in scored]}}
    for k, n in enumerate(j.coupling):
        fit["couplings"][prob.key_of.get(n, n)] = {"J_at_x": j.coupling_values(z, k).tolist()}
    for s, e in enumerate(prob.series):
        x = j.spectrum_vector(z, s)
        fit["spectrum_parameters"][e["id"]] = {n: float(x[j.col[n]]) for n in list(j.shared) + list(j.local)}
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(fit, fh, indent=1)
    return fit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="output directory of a fit_joint_series run (with monitor/)")
    ap.add_argument("out", help="fit.json-like file to write")
    args = ap.parse_args()
    fit = snapshot(args.run, args.out)
    snap = fit["snapshot"]
    print(f"snapshot {fit['scores'][0]:.6g} (record {snap['record']}, {len(snap['candidates'])} candidates) -> "
          f"{args.out}")


if __name__ == "__main__":
    main()
