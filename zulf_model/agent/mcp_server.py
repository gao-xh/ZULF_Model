"""Optional MCP (Model Context Protocol) server over STDIO.

Exposes every registry tool with its JSON schema, so Claude Code, Claude
Desktop, Codex or any MCP client can call them. Requires the `mcp` package
(`pip install "mcp>=1.12,<2"`). Register with, for example:

    claude mcp add zulf-model -- python -m zulf_model.agent.mcp_server
"""
from __future__ import annotations

import asyncio
import json


def build_server():
    import mcp.types as types
    from mcp.server import Server

    from .tools import REGISTRY

    server = Server("zulf-model")

    @server.list_tools()
    async def list_tools():
        return [types.Tool(name=t.name, description=t.description, inputSchema=t.input_schema)
                for t in REGISTRY.tools.values()]

    @server.call_tool()
    async def call_tool(name: str, arguments: dict):
        try:
            result = await asyncio.to_thread(REGISTRY.call, name, arguments or {})
            text = json.dumps(result, default=float)
        except Exception as exc:
            text = json.dumps({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
        return [types.TextContent(type="text", text=text)]

    return server


async def _serve():
    from mcp.server.stdio import stdio_server
    server = build_server()
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


def main():
    asyncio.run(_serve())


if __name__ == "__main__":
    main()
