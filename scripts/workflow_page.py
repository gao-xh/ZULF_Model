"""HTML page of the workflow map: docs/WORKFLOW.md (map and W1-W7) rendered as one self-contained page.

    python scripts/workflow_page.py [--source docs/WORKFLOW.md] [--out runs/workflow/index.html]

The page is generated from the document, so it never disagrees with it: regenerate after WORKFLOW.md changes.
Mermaid blocks are emitted as <pre class="mermaid"> (rendered by the artifact viewer natively, or by mermaid.js
when --mermaid-cdn is given for a local browser). Only the part before "W3 in detail" is included.
"""
import argparse
import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MERMAID_CDN = "https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/mermaid.min.js"
MATHJAX_CDN = "https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"      # SVG output: no fonts or CSS to load
STATUS = {"W1": "in use", "W2": "in use", "W3": "in use", "W4": "interactive", "W5": "in use", "W6": "in use",
          "W7": "not trained"}


def inline(text: str) -> str:
    """Markdown inline subset: `code`, **bold**."""
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    return re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", out)


def blocks(lines):
    """Yield (kind, payload) for a markdown section: mermaid, command, list, table, paragraph."""
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```mermaid") or line.startswith("```math"):
            j = i + 1
            while not lines[j].startswith("```"):
                j += 1
            yield line[3:].strip(), "\n".join(lines[i + 1:j])
            i = j + 1
        elif line.startswith("    "):
            j = i
            while j < len(lines) and (lines[j].startswith("    ") or not lines[j].strip()):
                j += 1
            yield "command", "\n".join(l[4:] for l in lines[i:j]).strip("\n")
            i = j
        elif line.startswith("- "):
            items = []
            while i < len(lines) and (lines[i].startswith("- ") or lines[i].startswith("  ")):
                if lines[i].startswith("- "):
                    items.append(lines[i][2:])
                else:
                    items[-1] += " " + lines[i].strip()
                i += 1
            yield "list", items
        elif line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(set(c) <= set("-: ") for c in cells):
                    rows.append(cells)
                i += 1
            yield "table", rows
        elif line.strip():
            para = [line]
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("- ", "|", "    ", "```")):
                para.append(lines[i])
                i += 1
            yield "paragraph", " ".join(para)
        else:
            i += 1


def render_blocks(lines) -> str:
    out = []
    algo_title = None                      # a "Key algorithms ..." paragraph: its list becomes one boxed block
    for kind, payload in blocks(lines):
        if kind == "paragraph" and payload.startswith("Key algorithms"):
            algo_title = payload.rstrip(":")
            continue
        if algo_title and kind == "list":
            head = inline(algo_title).replace("Key algorithms", "<span>Key algorithms</span>", 1)
            out.append(f'<div class="algo"><h3>{head}</h3><ul>' + "".join(f"<li>{inline(t)}</li>" for t in payload)
                       + "</ul></div>")
            algo_title = None
            continue
        algo_title = None
        if kind == "math":
            out.append(f'<div class="math">\\[{html.escape(payload, quote=False)}\\]</div>')
        elif kind == "mermaid":
            out.append(f'<div class="diagram"><pre class="mermaid">{html.escape(payload, quote=False)}</pre></div>')
        elif kind == "command":
            out.append(f'<pre class="cmd"><code>{html.escape(payload, quote=False)}</code></pre>')
        elif kind == "list":
            out.append("<ul>" + "".join(f"<li>{inline(t)}</li>" for t in payload) + "</ul>")
        elif kind == "table":
            head, *body = payload
            rows = "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in body)
            out.append('<div class="table-wrap"><table><thead><tr>' + "".join(f"<th>{inline(c)}</th>" for c in head)
                       + f"</tr></thead><tbody>{rows}</tbody></table></div>")
        else:
            label = re.match(r"^([A-Z][a-z]+(?: \([^)]*\))?):", payload)
            text = inline(payload)
            if label and len(label.group(1)) < 80:
                text = f'<span class="label">{inline(label.group(1))}</span>' + inline(payload[len(label.group(1)) + 1:])
            out.append(f"<p>{text}</p>")
    return "\n".join(out)


def build(source: Path, mermaid_cdn: bool) -> str:
    text = source.read_text()
    text = text.split("## W3 in detail")[0].rstrip().rstrip("-").rstrip()
    text = re.sub(r"\b(sections? [\d, -]+?) below", r"WORKFLOW.md \1", text)    # the method stays in the document
    text = text.replace("in the second part of this file", "in the second part of docs/WORKFLOW.md")
    text = text.replace("a full command at the end of this file", "a full command at the end of docs/WORKFLOW.md")
    parts = re.split(r"^## ", text, flags=re.M)
    intro_lines = parts[0].splitlines()[1:]
    reviewed = next((l for l in intro_lines if l.startswith("Last reviewed")), "")
    intro = [l for l in intro_lines if not l.startswith("Last reviewed")]
    sections = []
    for part in parts[1:]:
        title, *body = part.splitlines()
        m = re.match(r"(W\d)\. (.*)", title)
        sections.append((m.group(1) if m else "map", m.group(2) if m else title, body))
    nav = "".join(f'<a href="#{key.lower()}"><b>{key}</b> {html.escape(title.split(":")[0])}</a>'
                  for key, title, _ in sections if key != "map")
    body_html = []
    for key, title, body in sections:
        sid = key.lower()
        if key == "map":
            body_html.append(f'<section id="map" class="map"><h2>Map</h2>{render_blocks(body)}</section>')
            continue
        status = STATUS.get(key, "")
        cls = {"in use": "ok", "interactive": "alt", "not trained": "todo"}.get(status, "ok")
        body_html.append(
            f'<section id="{sid}" class="flow"><header><span class="key">{key}</span>'
            f'<h2>{inline(title)}</h2><span class="chip {cls}">{status}</span></header>{render_blocks(body)}</section>')
    mathjax = ('<script>window.MathJax = {tex: {inlineMath: [["$", "$"]], displayMath: [["\\\\[", "\\\\]"]]}, '
               'svg: {fontCache: "global"}};</script>'
               f'<script src="{MATHJAX_CDN}"></script>') if "```math" in text or "$" in text else ""
    script = mathjax + (f'<script src="{MERMAID_CDN}"></script><script>mermaid.initialize({{startOnLoad: true, '
              f'theme: matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "neutral"}});</script>'
              if mermaid_cdn else "")
    return PAGE.format(intro=render_blocks(intro), reviewed=inline(reviewed), nav=nav, body="\n".join(body_html),
                       script=script)


PAGE = """<title>ZULF Workflow Map</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,500;6..72,600&family=Public+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: a lab-notebook column; one map up top, then one block per workflow (diagram, command, notes). */
:root {{
  --bg: #f6f7f5; --panel: #ffffff; --ink: #1d2422; --muted: #5d6965; --line: #d9dfdc;
  --accent: #19657a; --accent-soft: #e2eef1; --ok: #2f7a4d; --ok-soft: #e3f1e8; --alt: #7a5a19;
  --alt-soft: #f4ecd9; --todo: #8a3b3b; --todo-soft: #f5e3e3; --code-bg: #eef1ef;
  --display: "Newsreader", Georgia, serif; --body: "Public Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #121716; --panel: #1a201f; --ink: #e3e9e6; --muted: #9aa7a2; --line: #2c3533;
  --accent: #6cc0d4; --accent-soft: #19333a; --ok: #7fcb9a; --ok-soft: #1b3326; --alt: #e0bd72;
  --alt-soft: #3a2f17; --todo: #e59a9a; --todo-soft: #3d2222; --code-bg: #202826; color-scheme: dark; }} }}
:root[data-theme="dark"] {{
  --bg: #121716; --panel: #1a201f; --ink: #e3e9e6; --muted: #9aa7a2; --line: #2c3533;
  --accent: #6cc0d4; --accent-soft: #19333a; --ok: #7fcb9a; --ok-soft: #1b3326; --alt: #e0bd72;
  --alt-soft: #3a2f17; --todo: #e59a9a; --todo-soft: #3d2222; --code-bg: #202826; color-scheme: dark; }}
* {{ box-sizing: border-box; }}
body {{ background: var(--bg); color: var(--ink); font: 15px/1.6 var(--body); margin: 0; }}
.wrap {{ max-width: 1040px; margin: 0 auto; padding-inline: 16px; padding-block: 32px 64px;
  display: grid; gap: 28px; }}
.top h1 {{ font: 600 clamp(28px, 4vw, 40px)/1.15 var(--display); margin: 0 0 8px; text-wrap: balance; }}
.top .eyebrow {{ font: 500 12px var(--mono); letter-spacing: .08em; text-transform: uppercase; color: var(--accent); }}
.top p {{ max-width: 70ch; margin: 8px 0; color: var(--muted); }}
.top .reviewed {{ font: 12px var(--mono); color: var(--muted); }}
nav {{ position: sticky; top: env(safe-area-inset-top, 0px); z-index: 2; background: var(--bg);
  display: flex; flex-wrap: wrap; gap: 6px; padding-block: 10px; border-bottom: 1px solid var(--line); }}
nav a {{ text-decoration: none; color: var(--ink); font-size: 13px; padding: 4px 10px; border: 1px solid var(--line);
  border-radius: 999px; background: var(--panel); }}
nav a b {{ font-family: var(--mono); color: var(--accent); font-weight: 500; }}
nav a:hover, nav a:focus-visible {{ border-color: var(--accent); outline: none; }}
section {{ display: grid; gap: 14px; min-width: 0; scroll-margin-top: 64px; }}
section.flow {{ background: var(--panel); border: 1px solid var(--line); border-radius: 10px; padding: 22px; }}
section h2 {{ font: 600 22px/1.25 var(--display); margin: 0; text-wrap: balance; }}
section header {{ display: flex; align-items: baseline; flex-wrap: wrap; gap: 10px; }}
.key {{ font: 500 14px var(--mono); color: var(--accent); background: var(--accent-soft); padding: 2px 8px;
  border-radius: 6px; }}
.chip {{ font: 500 11px var(--mono); letter-spacing: .06em; text-transform: uppercase; padding: 2px 8px;
  border-radius: 999px; margin-left: auto; }}
.chip.ok {{ color: var(--ok); background: var(--ok-soft); }}
.chip.alt {{ color: var(--alt); background: var(--alt-soft); }}
.chip.todo {{ color: var(--todo); background: var(--todo-soft); }}
.diagram {{ overflow-x: auto; background: var(--bg); border: 1px solid var(--line); border-radius: 8px; padding: 12px; }}
.map .diagram {{ background: var(--panel); }}
pre.mermaid {{ margin: 0; display: flex; justify-content: center; font-family: var(--body); }}
pre.cmd {{ margin: 0; overflow-x: auto; background: var(--code-bg); border-radius: 8px; padding: 12px 14px;
  font: 13px/1.55 var(--mono); }}
code {{ font-family: var(--mono); font-size: .9em; background: var(--code-bg); padding: 0 4px; border-radius: 4px; }}
pre.cmd code {{ background: none; padding: 0; font-size: inherit; }}
ul {{ margin: 0; padding-left: 20px; display: grid; gap: 4px; }}
p {{ margin: 0; max-width: 78ch; }}
.label {{ font-weight: 600; color: var(--accent); margin-right: 4px; }}
.math {{ overflow-x: auto; padding: 4px 0; font-size: 15px; }}
.math mjx-container {{ margin: 6px 0 !important; }}
.algo {{ border-top: 1px solid var(--line); padding-top: 14px; display: grid; gap: 8px; }}
.algo h3 {{ margin: 0; font: 500 12px var(--mono); letter-spacing: .06em; color: var(--muted); }}
.algo h3 span {{ text-transform: uppercase; color: var(--accent); }}
.algo li {{ font-size: 14px; }}
.table-wrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13.5px; }}
th, td {{ text-align: left; padding: 8px 10px; border-bottom: 1px solid var(--line); vertical-align: top; }}
th {{ font: 500 11px var(--mono); letter-spacing: .06em; text-transform: uppercase; color: var(--muted); }}
td:first-child {{ white-space: nowrap; font-weight: 600; }}
footer {{ font-size: 13px; color: var(--muted); }}
@media (max-width: 560px) {{ section.flow {{ padding: 16px; }} .chip {{ margin-left: 0; }} }}
</style>
<div class="wrap">
  <div class="top">
    <div class="eyebrow">ZULF_Model &middot; docs/WORKFLOW.md</div>
    <h1>ZULF workflows, from scans to structures</h1>
    {intro}
    <p class="reviewed">{reviewed}</p>
  </div>
  <nav aria-label="Workflows"><a href="#map"><b>Map</b></a>{nav}</nav>
  {body}
  <footer>Generated from docs/WORKFLOW.md by scripts/workflow_page.py. The step-by-step fitting method (sections
  1-10) stays in the document.</footer>
</div>
{script}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=str(ROOT / "docs" / "WORKFLOW.md"))
    ap.add_argument("--out", default=str(ROOT / "runs" / "workflow" / "index.html"))
    ap.add_argument("--mermaid-cdn", action="store_true", help="load mermaid.js (for opening the file locally)")
    args = ap.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(Path(args.source), args.mermaid_cdn))
    print(out)


if __name__ == "__main__":
    main()
