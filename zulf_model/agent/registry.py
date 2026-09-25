"""Typed tool registry for AI agents (Claude, GPT) and the command line.

Every tool is a pure function taking a JSON-compatible dict and returning a
JSON-compatible dict. Large arrays are written to files inside the workspace
and returned as absolute paths. Input data paths are read-only. The same
registry backs the JSON CLI, the optional MCP server and the exported
Anthropic / OpenAI tool definitions, so behaviour is identical everywhere.
"""
from __future__ import annotations

import json
import os
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

WORKSPACE_ENV = "ZULF_MODEL_WORKSPACE"


def workspace() -> Path:
    root = Path(os.environ.get(WORKSPACE_ENV, "runs")).resolve()
    root.mkdir(parents=True, exist_ok=True)
    return root


def new_run_dir(tool: str) -> Path:
    path = workspace() / "tool_runs" / f"{tool}_{uuid.uuid4().hex[:12]}"
    path.mkdir(parents=True)
    return path


@dataclass
class Tool:
    name: str
    description: str
    input_schema: dict
    handler: Callable[[dict], dict]
    long_running: bool = False
    tags: List[str] = field(default_factory=list)

    def validate(self, arguments: dict) -> dict:
        schema = self.input_schema
        props = schema.get("properties", {})
        missing = [k for k in schema.get("required", []) if k not in arguments]
        if missing:
            raise ValueError(f"{self.name}: missing required arguments {missing}")
        unknown = {k for k in arguments if not k.startswith("_")} - set(props)
        if unknown and not schema.get("additionalProperties", False):
            raise ValueError(f"{self.name}: unknown arguments {sorted(unknown)}")
        out = dict(arguments)
        for key, spec in props.items():
            if key not in out and "default" in spec:
                out[key] = spec["default"]
        return out


class ToolRegistry:
    def __init__(self):
        self.tools: Dict[str, Tool] = {}

    def register(self, tool: Tool) -> Tool:
        if tool.name in self.tools:
            raise ValueError(f"Tool {tool.name} already registered.")
        self.tools[tool.name] = tool
        return tool

    def tool(self, name: str, description: str, input_schema: dict, long_running: bool = False, tags=()):
        def decorator(func):
            self.register(Tool(name, description, input_schema, func, long_running, list(tags)))
            return func
        return decorator

    def call(self, name: str, arguments: Optional[dict] = None) -> dict:
        if name not in self.tools:
            raise KeyError(f"Unknown tool '{name}'. Available: {sorted(self.tools)}")
        tool = self.tools[name]
        result = tool.handler(tool.validate(arguments or {}))
        json.dumps(result)  # results must be JSON-serializable
        return result

    def names(self) -> List[str]:
        return sorted(self.tools)

    # -- exports -----------------------------------------------------------------------
    def anthropic_tools(self) -> List[dict]:
        """Definitions for the Anthropic Messages API `tools` parameter."""
        return [{"name": t.name, "description": t.description, "input_schema": t.input_schema}
                for t in self.tools.values()]

    def openai_tools(self) -> List[dict]:
        """Definitions for OpenAI function calling (`tools` with type function)."""
        return [{"type": "function", "function": {"name": t.name, "description": t.description,
                                                  "parameters": t.input_schema}}
                for t in self.tools.values()]


REGISTRY = ToolRegistry()
