"""Fit-trace viewer: a self-contained HTML page with a slider over the recorded frames of a fit (fit_joint_series
--trace N): data and model (real part, and imaginary part for complex data), residual, the objective along the
path and the couplings of the current frame.

    python scripts/trace_view.py runs/processed/RUN [--out RUN/trace.html] [--band 115,140]
"""
import argparse
import json
from pathlib import Path

import numpy as np

PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fit trace viewer</title>
<style>
:root { --bg: #ffffff; --fg: #1d2330; --muted: #6b7280; --data: #1d2330; --model: #c0392b; --res: #9aa3ad;
        --grid: #e5e7eb; --accent: #2563eb; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg: #12151b; --fg: #e5e7eb;
        --muted: #9aa3ad; --data: #e5e7eb; --model: #f87171; --res: #6b7280; --grid: #262b35; --accent: #60a5fa; } }
:root[data-theme="dark"] { --bg: #12151b; --fg: #e5e7eb; --muted: #9aa3ad; --data: #e5e7eb; --model: #f87171;
        --res: #6b7280; --grid: #262b35; --accent: #60a5fa; }
body { background: var(--bg); color: var(--fg); font: 14px/1.4 system-ui, sans-serif; margin: 0; padding: 16px; }
h1 { font-size: 18px; margin: 0 0 4px; } .muted { color: var(--muted); }
.row { display: flex; gap: 16px; flex-wrap: wrap; align-items: center; margin: 8px 0; }
input[type=range] { flex: 1; min-width: 200px; } canvas { width: 100%; height: 320px; display: block; }
#cost { height: 140px; } table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
td, th { padding: 2px 10px; text-align: right; border-bottom: 1px solid var(--grid); } th { color: var(--muted); }
button { background: none; border: 1px solid var(--grid); color: var(--fg); padding: 4px 10px; border-radius: 6px; }
</style></head><body>
<h1>Fit trace</h1><div class="muted" id="sub"></div>
<div class="row"><button id="play">Play</button><input id="slider" type="range" min="0" value="0">
<span id="label"></span></div>
<div class="row"><label>spectrum <select id="spec"></select></label>
<label>from <input id="lo" size="6"> to <input id="hi" size="6"> Hz</label>
<label><input type="checkbox" id="imag"> imaginary part</label></div>
<canvas id="plot"></canvas><canvas id="cost"></canvas>
<table id="jt"></table>
<script>
const T = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const slider = document.getElementById("slider"), label = document.getElementById("label");
slider.max = T.frames.length - 1;
document.getElementById("sub").textContent = T.origin + ", " + T.evaluations + " evaluations, " + T.frames.length + " frames";
const spec = document.getElementById("spec");
T.spectra.forEach((s, i) => { const o = document.createElement("option"); o.value = i; o.textContent = s.id + " (x " + s.x + ")"; spec.appendChild(o); });
document.getElementById("lo").value = T.band[0]; document.getElementById("hi").value = T.band[1];
function canvas(id) { const c = document.getElementById(id), r = window.devicePixelRatio || 1;
  c.width = c.clientWidth * r; c.height = c.clientHeight * r; const g = c.getContext("2d"); g.scale(r, r);
  return [g, c.clientWidth, c.clientHeight]; }
function draw() {
  const k = +slider.value, s = +spec.value, part = document.getElementById("imag").checked ? 1 : 0;
  const lo = +document.getElementById("lo").value, hi = +document.getElementById("hi").value;
  const fr = T.frames[k];
  label.textContent = "frame " + (k + 1) + "/" + T.frames.length + ", evaluation " + fr.evaluation + ", " + fr.stage + ", objective " + fr.objective.toPrecision(5);
  const f = T.f[s], y = T.y[s][part], m = T.model[s][k][part];
  const idx = []; for (let i = 0; i < f.length; i++) if (f[i] >= lo && f[i] <= hi) idx.push(i);
  let ymin = Infinity, ymax = -Infinity;
  idx.forEach(i => { ymin = Math.min(ymin, y[i], m[i], y[i] - m[i] - 0.0); ymax = Math.max(ymax, y[i], m[i]); });
  const span = ymax - ymin || 1; ymin -= 0.35 * span;
  const [g, W, H] = canvas("plot"); g.clearRect(0, 0, W, H);
  const X = v => 40 + (v - lo) / (hi - lo) * (W - 50), Y = v => 10 + (ymax - v) / (ymax - ymin) * (H - 30);
  g.strokeStyle = css("--grid"); g.beginPath(); g.moveTo(40, H - 20); g.lineTo(W - 10, H - 20); g.stroke();
  g.fillStyle = css("--muted"); g.font = "11px system-ui";
  for (let t = Math.ceil(lo); t <= hi; t += Math.max(1, Math.round((hi - lo) / 10))) { g.fillText(t, X(t) - 8, H - 6); }
  const line = (arr, color, off) => { g.strokeStyle = color; g.lineWidth = 1.2; g.beginPath();
    idx.forEach((i, j) => { const yy = Y(arr(i) + off); j ? g.lineTo(X(f[i]), yy) : g.moveTo(X(f[i]), yy); }); g.stroke(); };
  line(i => y[i], css("--data"), 0); line(i => m[i], css("--model"), 0);
  const roff = ymin + 0.15 * span; line(i => y[i] - m[i], css("--res"), roff);
  g.fillStyle = css("--data"); g.fillText("data", W - 120, 20); g.fillStyle = css("--model"); g.fillText("model", W - 80, 20);
  g.fillStyle = css("--res"); g.fillText("residual (offset)", W - 120, 34);
  const [c, CW, CH] = canvas("cost"); c.clearRect(0, 0, CW, CH);
  const obj = T.frames.map(q => Math.log10(q.objective)); const omin = Math.min(...obj), omax = Math.max(...obj);
  const CX = j => 40 + j / Math.max(T.frames.length - 1, 1) * (CW - 50), CY = v => 8 + (omax - v) / ((omax - omin) || 1) * (CH - 26);
  c.strokeStyle = css("--accent"); c.beginPath(); obj.forEach((v, j) => j ? c.lineTo(CX(j), CY(v)) : c.moveTo(CX(j), CY(v))); c.stroke();
  c.fillStyle = css("--model"); c.beginPath(); c.arc(CX(k), CY(obj[k]), 4, 0, 7); c.fill();
  c.fillStyle = css("--muted"); c.fillText("log10 objective along the path", 44, CH - 6);
  const jt = document.getElementById("jt"); const first = T.frames[0].J;
  jt.innerHTML = "<tr><th>coupling</th><th>this frame (Hz)</th><th>start</th><th>final</th></tr>" +
    Object.keys(fr.J).map(n => "<tr><td>" + n + "</td><td>" + fr.J[n].map(v => v.toFixed(3)).join(", ") + "</td><td>" +
      first[n].map(v => v.toFixed(2)).join(", ") + "</td><td>" + T.frames[T.frames.length - 1].J[n].map(v => v.toFixed(2)).join(", ") + "</td></tr>").join("");
}
let timer = null; document.getElementById("play").onclick = () => {
  if (timer) { clearInterval(timer); timer = null; return; }
  timer = setInterval(() => { slider.value = (+slider.value + 1) % T.frames.length; draw(); }, 150); };
["slider", "spec", "lo", "hi", "imag"].forEach(id => document.getElementById(id).addEventListener("input", draw));
window.addEventListener("resize", draw); draw();
</script></body></html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--out", default="")
    ap.add_argument("--band", default="")
    ap.add_argument("--digits", type=int, default=5, help="significant digits kept for the spectra")
    args = ap.parse_args()
    run = Path(args.run)
    meta = json.loads((run / "trace.json").read_text())
    arrays = np.load(run / "trace.npz")
    n = len(meta["spectra"])
    scale = [float(np.abs(arrays[f"y{s}"]).max()) or 1.0 for s in range(n)]

    def two(v, s):
        v = np.asarray(v) / scale[s]
        return [np.round(v.real, args.digits).tolist(), np.round(v.imag, args.digits).tolist()]
    f0 = arrays["f0"]
    band = [float(v) for v in args.band.split(",")] if args.band else [float(f0.min()), float(f0.max())]
    data = {"origin": meta["origin"], "evaluations": meta["evaluations"], "frames": meta["frames"],
            "spectra": meta["spectra"], "band": band,
            "f": [np.round(arrays[f"f{s}"], 4).tolist() for s in range(n)],
            "y": [two(arrays[f"y{s}"], s) for s in range(n)],
            "model": [[two(m, s) for m in arrays[f"model{s}"]] for s in range(n)]}
    out = Path(args.out) if args.out else run / "trace.html"
    out.write_text(PAGE.replace("__DATA__", json.dumps(data)))
    print(out)


if __name__ == "__main__":
    main()
