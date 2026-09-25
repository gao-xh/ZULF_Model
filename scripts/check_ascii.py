"""Fail when tracked source, config or documentation files contain non-ASCII text.

Usage: python scripts/check_ascii.py [root]
"""
import sys
from pathlib import Path

SUFFIXES = {".py", ".md", ".json", ".toml", ".yaml", ".yml", ".txt", ".cfg"}
SKIP_DIRS = {".git", "runs", "__pycache__", ".venv", "build", "dist"}


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1])
    bad = []
    for path in root.rglob("*"):
        if any(part in SKIP_DIRS or part.endswith(".egg-info") for part in path.parts) or path.suffix not in SUFFIXES:
            continue
        if not path.is_file():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if any(ord(ch) > 127 for ch in line):
                bad.append(f"{path.relative_to(root)}:{number}: {line.strip()[:80]!r}")
    for entry in bad:
        print(entry)
    print(f"{len(bad)} non-ASCII line(s).")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
