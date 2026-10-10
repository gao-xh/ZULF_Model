"""Unattended batch of known-structure fits: scans on the drive -> averages -> whole-grid spectrum -> two fits.

    python scripts/run_batch.py configs/batch_2026-10-10.json [--parallel 3] [--workers 3] [--only ID,...]

Per entry (resumable: a step whose output exists is skipped):
  1. averages of the scans (scripts/average_scans.py, read from the drive): `all-scans/` and `z5/`
     (--exclude-z 5) in ~/research/zulf/data/processed/<measurement>/[<scans folder>/] (the folder name is added
     when one measurement holds several runs, e.g. a dilution series);
  2. the whole-grid spectrum of the z5 average (scripts/make_series_entry.py, `--exclude 81.5,86`; entries may set
     `grid` and `sampling_rate_factor`, e.g. 0.5 for double-acquisition data): runs/series/batch_<id>/;
  3. stage 1 fit (fit_joint_series.py): one decay rate per isotopologue, component search at the start, 6 starts,
     field fitted only when the entry says `field` (else zero field): runs/processed/batch_<id>_s1/;
  4. stage 2 fit from stage 1 (--from-joint): one decay rate per model line (--family-edges lines), rate bounds
     0.1-40 1/s, 2 starts: runs/processed/batch_<id>_s2/.
Waits (up to --wait hours) for an entry's scans folder to appear on the drive. Writes a results table
(docs/analysis/2026-10-10_batch_results.md by default) after every finished entry; a failed step is recorded and
the entry skipped. Fit results are conditional numerical results of the given structure, not assignments.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
DRIVE = Path("/Volumes/Extreme/research/zulf/data")
PROCESSED = Path.home() / "research/zulf/data/processed"
LOCK = threading.Lock()
STATE: dict = {}


def measurement_of(scans: str):
    parts = Path(scans).parts                      # original|raw / <measurement> / <folder> [/ data]
    return parts[1], parts[2]


def run(cmd, log: Path) -> int:
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        return subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode


def note(eid, **kw):
    with LOCK:
        STATE.setdefault(eid, {}).update(kw)


def summary(fit_json: Path) -> dict:
    d = json.loads(fit_json.read_text())
    sp = next(iter(d.get("spectrum_parameters", {}).values()), {})
    out = {"score": d["scores"][0], "residual": next(iter(d.get("data_region_residuals", {}).values()), None),
           "J": {k: (v["J_at_x"][0], v.get("J_std_at_x", [float("nan")])[0]) for k, v in d["couplings"].items()}}
    if "field_z_ut" in sp or "field_transverse_ut" in sp:
        out["field_nt"] = (1e3 * sp.get("field_transverse_ut", 0.0), 1e3 * sp.get("field_z_ut", 0.0))
    rates = [math.exp(v) for k, v in sp.items() if ".log_rate" in k]
    out["rates_at_bound"] = sum(1 for r in rates if r > 39.5 or r < 0.105)
    return out


def process(entry, args):
    eid = entry["id"]
    scans = DRIVE / entry["scans"]
    _, meas = measurement_of(entry["scans"])
    log = ROOT / "runs/processed" / f"batch_{eid}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    while not scans.is_dir():
        if time.time() - t0 > args.wait * 3600:
            note(eid, status="no data on the drive")
            return
        note(eid, status="waiting for data")
        time.sleep(60)
    folder = Path(entry["scans"]).parts[3] if "series" in meas or "dilution" in meas else ""
    proc = PROCESSED / meas / folder
    note(eid, status="averaging", measurement=meas)
    for name, extra in (("z5", ["--exclude-z", "5"]), ("all-scans", [])):
        if not (proc / name / "average_fid.npy").exists():
            if run([PY, "scripts/average_scans.py", str(scans), str(proc / name), *extra], log):
                note(eid, status=f"failed: average {name}")
                return
    series = ROOT / "runs/series" / f"batch_{eid}"
    if not (series / "series.json").exists():
        note(eid, status="spectrum")
        cmd = [PY, "scripts/make_series_entry.py", "--fid", str(proc / "z5/average_fid.npy"), "--id", eid,
               "--out", str(series.relative_to(ROOT)), "--exclude", "81.5,86"]
        if entry.get("grid"):
            cmd += ["--grid", entry["grid"]]
        if entry.get("sampling_rate_factor"):
            scans_json = json.loads((proc / "z5/scans.json").read_text())
            fs = float(scans_json.get("sampling_rate_hz") or scans_json.get("settings", {}).get("sampling_rate_hz"))
            cmd += ["--sampling-rate", f"{fs * entry['sampling_rate_factor']:.6f}"]
        if run(cmd, log):
            note(eid, status="failed: spectrum")
            return
    lo, hi = (entry.get("grid") or "20,380").split(",")
    common = [PY, "scripts/fit_joint_series.py", "--series", str((series / "series.json").relative_to(ROOT)),
              "--real-only", "false", "--shape", "free", "--exchange", "fast", "--range", f"{lo},{hi}",
              "--structure", json.dumps(entry["structure"]), "--signal-threshold", "2.5", "--signal-taper", "4.0",
              "--peak-penalty", "5", "--peak-smooth", "0.03", "--peak-min-sigma", "3",
              "--workers", str(args.workers), "--max-nfev", "300"]
    if entry.get("field"):
        common += ["--fit-field", "--field-start", entry.get("field_start", "0.037,0.045")]
    s1 = ROOT / "runs/processed" / f"batch_{eid}_s1"
    if not (s1 / "fit.json").exists():
        note(eid, status="stage 1 fit")
        if run(common + ["--starts", "6", "--rate-bounds", "0.2,15", "--component-search", "start",
                         "--out", str(s1.relative_to(ROOT))], log):
            note(eid, status="failed: stage 1")
            return
    note(eid, s1=summary(s1 / "fit.json"))
    s2 = ROOT / "runs/processed" / f"batch_{eid}_s2"
    if not (s2 / "fit.json").exists():
        note(eid, status="stage 2 fit")
        if run(common + ["--starts", "2", "--rate-bounds", "0.1,40", "--component-search", "off",
                         "--family-edges", "lines", "--from-joint", str((s1 / "fit.json").relative_to(ROOT)),
                         "--out", str(s2.relative_to(ROOT))], log):
            note(eid, status="failed: stage 2")
            return
    note(eid, s2=summary(s2 / "fit.json"), status="done", minutes=round((time.time() - t0) / 60))


def write_table(path: Path, entries):
    rows = ["# Overnight batch 2026-10-10: results", "",
            "Conditional numerical results of the given structures (scripts/run_batch.py, configs/batch_2026-10-10.json);",
            "stage 1: one rate per isotopologue; stage 2: one rate per model line. J in Hz (J_std: linearised).", "",
            "| id | status | stage 1 score | stage 2 score / residual | field (nT, transverse / z) | couplings (stage 2) |",
            "|---|---|---|---|---|---|"]
    with LOCK:
        for e in entries:
            st = STATE.get(e["id"], {})
            s1, s2 = st.get("s1"), st.get("s2")
            f = (s2 or s1 or {}).get("field_nt")
            jtxt = ", ".join(f"{k} {v[0]:.2f}" for k, v in (s2 or s1 or {}).get("J", {}).items())
            rows.append(f"| {e['id']} | {st.get('status', 'queued')} | {s1['score']:.4f} |" if s1 else
                        f"| {e['id']} | {st.get('status', 'queued')} | - |")
            rows[-1] += (f" {s2['score']:.4f} / {s2['residual']:.3f} |" if s2 else " - |")
            rows[-1] += (f" {f[0]:.1f} / {f[1]:.1f} |" if f else " zero field |" if not e.get("field") else " - |")
            rows[-1] += f" {jtxt} |"
    path.write_text("\n".join(rows) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plan")
    ap.add_argument("--parallel", type=int, default=3)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--wait", type=float, default=6.0, help="hours to wait for an entry's data on the drive")
    ap.add_argument("--only", default="")
    ap.add_argument("--table", default="docs/analysis/2026-10-10_batch_results.md")
    args = ap.parse_args()
    entries = json.loads(Path(args.plan).read_text())["entries"]
    if args.only:
        keep = set(args.only.split(","))
        entries = [e for e in entries if e["id"] in keep]
    table = ROOT / args.table

    def job(e):
        try:
            process(e, args)
        except Exception as exc:                       # one entry failing never stops the batch
            note(e["id"], status=f"failed: {type(exc).__name__}: {exc}"[:200])
        write_table(table, entries)
        print(f"{e['id']}: {STATE.get(e['id'], {}).get('status')}", flush=True)
    write_table(table, entries)
    with ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(job, entries))
    write_table(table, entries)


if __name__ == "__main__":
    main()
