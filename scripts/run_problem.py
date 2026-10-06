"""Rebuild the fit problem of a fit_joint_series run from its output directory.

    from run_problem import load_run
    run = load_run("runs/processed/NAME")      # run.prob (build_problem result at the fit), run.fit, run.argv

The command comes from RUN_LOG.md ("## Command" block), else from monitor/status.json. Options that point to
start files or outputs are dropped, and the run's own fit.json (or a given fit file) becomes --from-joint, so
prob.z0 is the fitted vector. Used by plot_runs.py, plot_components.py and snapshot_fit.py.
"""
import json
import shlex
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_joint_series import make_parser, build_problem  # noqa: E402

# options that point to start files or outputs, or only steer the multi-start; not needed to rebuild the problem
DROP = ("--seeds", "--from-joint", "--out", "--monitor", "--trace", "--starts", "--workers", "--seed", "--spread",
        "--max-nfev")


def run_argv(run_dir):
    """The fit_joint_series options of a run (script name removed)."""
    run_dir = Path(run_dir)
    log = run_dir / "RUN_LOG.md"
    if log.exists() and "## Command" in log.read_text():
        cmd = log.read_text().split("## Command")[1].split("```")[1]
        argv = shlex.split(cmd)
    else:
        with open(run_dir / "monitor" / "status.json") as fh:
            argv = json.load(fh)["argv"]
    if argv and argv[0].endswith(".py"):
        argv = argv[1:]
    elif argv and argv[0] == "python":
        argv = argv[2:]
    return argv


def problem_argv(argv):
    keep, skip = [], False
    for a in argv:
        if skip:
            skip = False
            continue
        if a in DROP:
            skip = True
            continue
        keep.append(a)
    return keep


def load_run(run_dir, fit_path=None):
    run_dir = Path(run_dir)
    argv = run_argv(run_dir)
    fit_path = Path(fit_path) if fit_path else run_dir / "fit.json"
    extra = ["--from-joint", str(fit_path)] if fit_path.exists() else []
    prob = build_problem(make_parser().parse_args(problem_argv(argv) + extra + ["--monitor", "off"]))
    fit = None
    if fit_path.exists():
        with open(fit_path) as fh:
            fit = json.load(fh)
    return SimpleNamespace(run_dir=run_dir, argv=argv, prob=prob, fit=fit, name=run_dir.name)


def spectrum_arrays(prob, s=0, z=None):
    """Frequencies, data, model, component spectra (complex, on the fitted points) and the prediction of
    spectrum s at z (default prob.z0)."""
    j = prob.joint
    f = j.forwards[s]
    z = prob.z0 if z is None else z
    pred = f.predict(j.spectrum_vector(z, s))
    nf = len(f.f)

    def cm(v):
        v = np.asarray(v)
        return v[:nf] + 1j * v[nf:] if len(v) == 2 * nf else v.astype(complex)
    return (np.asarray(f.f), cm(f.y), cm(pred.model), [cm(c) for c in pred.component_spectra], pred)


def segments(x, gap_hz=0.2):
    """Index runs of x without gaps wider than gap_hz (excluded ranges are not bridged in plots)."""
    return np.split(np.arange(len(x)), np.flatnonzero(np.diff(x) > gap_hz) + 1)
