"""Command line interface.

    zulf-model tools list
    zulf-model tools export --format anthropic|openai [--output FILE]
    zulf-model tool NAME --json '{...}'        (or --request FILE, or JSON on stdin)
    zulf-model submit NAME --json '{...}'      (background job; poll with `tool get_job`)
    zulf-model serve-mcp
    zulf-model train RUN_CONFIG [--max-steps N]
    zulf-model diagnose AVERAGE.npy 0.ini
    zulf-model throughput [--samples N]

Every command prints JSON on stdout; errors print {"status": "error", ...} and
exit with code 1.
"""
from __future__ import annotations

import argparse
import json
import sys


def _arguments(args) -> dict:
    if getattr(args, "json", None):
        return json.loads(args.json)
    if getattr(args, "request", None):
        with open(args.request, encoding="utf-8") as handle:
            return json.load(handle)
    if not sys.stdin.isatty():
        text = sys.stdin.read().strip()
        return json.loads(text) if text else {}
    return {}


def main(argv=None) -> int:
    from .agent import registry
    reg = registry()
    parser = argparse.ArgumentParser(prog="zulf-model", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    tools = sub.add_parser("tools")
    tools_sub = tools.add_subparsers(dest="action", required=True)
    tools_sub.add_parser("list")
    export = tools_sub.add_parser("export")
    export.add_argument("--format", choices=["anthropic", "openai"], default="anthropic")
    export.add_argument("--output")
    for name in ("tool", "submit"):
        p = sub.add_parser(name)
        p.add_argument("name")
        p.add_argument("--json")
        p.add_argument("--request")
    sub.add_parser("serve-mcp")
    train = sub.add_parser("train")
    train.add_argument("run_config")
    train.add_argument("--max-steps", type=int)
    diag = sub.add_parser("diagnose")
    diag.add_argument("average_npy")
    diag.add_argument("ini")
    thr = sub.add_parser("throughput")
    thr.add_argument("--samples", type=int, default=50)
    thr.add_argument("--points", type=int, default=16384)
    args = parser.parse_args(argv)
    try:
        if args.command == "tools" and args.action == "list":
            result = {"tools": [{"name": t.name, "long_running": t.long_running, "description": t.description}
                                for t in reg.tools.values()]}
        elif args.command == "tools":
            result = reg.anthropic_tools() if args.format == "anthropic" else reg.openai_tools()
            if args.output:
                with open(args.output, "w", encoding="utf-8") as handle:
                    json.dump(result, handle, indent=2)
                result = {"written": args.output, "count": len(result)}
        elif args.command == "tool":
            result = reg.call(args.name, _arguments(args))
        elif args.command == "submit":
            result = reg.call("submit_job", {"tool": args.name, "arguments": _arguments(args)})
        elif args.command == "serve-mcp":
            from .agent.mcp_server import main as serve
            serve()
            return 0
        elif args.command == "train":
            payload = {"run_config": args.run_config}
            if args.max_steps:
                payload["max_steps"] = args.max_steps
            result = reg.call("train_model", payload)
        elif args.command == "diagnose":
            result = reg.call("diagnose_fid", {"source": {"average_npy": args.average_npy, "ini": args.ini}})
        elif args.command == "throughput":
            from .throughput import measure
            result = measure(args.samples, args.points)
        print(json.dumps(result, indent=2, default=float))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "error", "error": f"{type(exc).__name__}: {exc}"}))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
