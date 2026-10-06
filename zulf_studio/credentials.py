"""API keys of the Studio AI assistant in the macOS Keychain (service "zulf-studio", one item per variable).

Keys are entered once in the AI assistant tab ("Remember in the macOS Keychain") and put into this process's
environment when Studio starts; a key already set in the environment wins. Nothing is written to a file or to
the log. Uses the system `security` tool; on other systems `available()` is False and nothing is stored.
Remove a key with Forget in the tab, or in the Keychain Access app (search "zulf-studio").
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from typing import Callable, List, Optional

SERVICE = "zulf-studio"
VARIABLES = ("ANTHROPIC_API_KEY", "OPENAI_API_KEY")


def _run(args: List[str], runner: Callable = subprocess.run):
    return runner(args, capture_output=True, text=True)


def available() -> bool:
    return sys.platform == "darwin" and shutil.which("security") is not None


def save(variable: str, key: str, runner: Callable = subprocess.run) -> bool:
    """Store (or replace) the key of `variable`; True on success."""
    if variable not in VARIABLES or not key:
        raise ValueError("unknown variable or empty key")
    r = _run(["security", "add-generic-password", "-U", "-s", SERVICE, "-a", variable, "-w", key], runner)
    return r.returncode == 0


def load(variable: str, runner: Callable = subprocess.run) -> Optional[str]:
    r = _run(["security", "find-generic-password", "-s", SERVICE, "-a", variable, "-w"], runner)
    key = r.stdout.strip() if r.returncode == 0 else ""
    return key or None


def delete(variable: str, runner: Callable = subprocess.run) -> bool:
    r = _run(["security", "delete-generic-password", "-s", SERVICE, "-a", variable], runner)
    return r.returncode == 0


def load_into_environment(runner: Callable = subprocess.run, environ=os.environ) -> List[str]:
    """Put remembered keys into the environment where the variable is not set yet; returns the variables set."""
    if runner is subprocess.run and not available():
        return []
    loaded = []
    for var in VARIABLES:
        if environ.get(var):
            continue
        key = load(var, runner)
        if key:
            environ[var] = key
            loaded.append(var)
    return loaded
