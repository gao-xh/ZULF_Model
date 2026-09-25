"""AI-agent interface: one tool registry for the JSON CLI, MCP and exported API tool schemas."""
from .registry import REGISTRY, Tool, ToolRegistry, workspace


def registry() -> ToolRegistry:
    from . import tools  # noqa: F401  (registers tools on import)
    return REGISTRY


__all__ = ["REGISTRY", "Tool", "ToolRegistry", "workspace", "registry"]
