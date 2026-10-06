"""AI assistant of ZULF Studio: a language model drives the session through the studio's own tools.

Providers:
- "anthropic": Claude through the Anthropic SDK (Messages API, client tools, manual loop). Credentials as the
  SDK resolves them: ANTHROPIC_API_KEY, ANTHROPIC_AUTH_TOKEN or an `ant auth login` profile. Default model
  claude-opus-5-5; the server-side refusal fallback is on (fallbacks="default").
- "openai": OpenAI models, e.g. Codex models, through the OpenAI SDK (Responses API, function calling).
  Credentials: OPENAI_API_KEY. The model name is required (argument, or OPENAI_MODEL).

The tools are the studio API tools (zulf_studio.api.TOOLS) executed in-process on the same session the window
shows; every call and its outcome goes to the session log (source "ai"). Large array results (simulate) are
summarised before they go back to the model. The model is told that fitted values are conditional numerical
results (AGENTS.md). Keys are never stored by the studio.
"""
from __future__ import annotations

import json
import os
import threading
from typing import Callable, List, Optional

from .api import StudioAPI

ANTHROPIC_MODEL = "claude-opus-5-5"
MAX_RESULT_CHARS = 12000

SYSTEM = """You operate ZULF Studio, a zero- to ultralow-field NMR simulator and fitter, through its tools.
The session holds a molecular structure (natural-abundance 13C isotopologues), couplings J (Hz), a static field
(B transverse and B z, nT; zero field is the default), one decay rate (1/s) for the quick-look simulation, an
optional processed experimental spectrum, fit jobs (fit_joint_series, the model rendered through the data
processing) with traces, and publication figures.
Work economically: read get_state first; prefer lines, fit_status and the scale / rms residual of simulate
over long arrays. Change one thing at a time and check its effect. A fit or figure runs in the background:
poll fit_status / figure_status. A fitted result is a conditional numerical result: report it with its objective,
residual, settings and remaining misfit, never as a determined molecular assignment. Say what you changed in
the session. Answer in the language of the user's request."""


def _compact(tool: str, result) -> str:
    """Tool result as text for the model; simulate arrays are replaced by a summary."""
    if tool == "simulate" and isinstance(result, dict):
        keep = {k: v for k, v in result.items() if k not in ("f", "sim_re", "sim_im", "data_re", "data_im")}
        keep["points"] = len(result.get("f", []))
        keep["lines"] = [{"component": r["component"], "frequency_hz": round(r["frequency_hz"], 4),
                          "relative": round(r["relative"], 4)} for r in result.get("lines", []) if r["relative"] >= 0.05]
        result = keep
    text = json.dumps(result)
    if len(text) > MAX_RESULT_CHARS:
        text = text[:MAX_RESULT_CHARS] + f'... [truncated, {len(text)} characters]'
    return text


class StudioAssistant:
    """One conversation with a model about the session. `ask` runs the tool loop (in a thread when
    background=True) and calls on_message(role, text) for the transcript: "user", "assistant", "tool", "error"."""

    def __init__(self, session, provider: str = "anthropic", model: str = "", max_steps: int = 30,
                 on_message: Optional[Callable[[str, str], None]] = None, client=None):
        if provider not in ("anthropic", "openai"):
            raise ValueError("provider must be 'anthropic' or 'openai'")
        self.session, self.provider, self.max_steps = session, provider, int(max_steps)
        self.model = model or (ANTHROPIC_MODEL if provider == "anthropic" else os.environ.get("OPENAI_MODEL", ""))
        if not self.model:
            raise ValueError("OpenAI needs a model name (field 'model' or OPENAI_MODEL), e.g. a Codex model "
                             "available to your account")
        self.api = StudioAPI(session)
        self.on_message = on_message or (lambda role, text: None)
        self._client = client
        self.history: List = []            # provider-specific conversation
        self.stop_requested = False
        self.running = False

    # ---- plumbing ------------------------------------------------------------------------------
    def client(self):
        if self._client is None:
            if self.provider == "anthropic":
                import anthropic
                self._client = anthropic.Anthropic()
            else:
                import openai
                self._client = openai.OpenAI()
        return self._client

    def _run_tool(self, name: str, args) -> tuple:
        try:
            if isinstance(args, str):
                args = json.loads(args or "{}")
            result = self.api.call(name, args or {})
            text, error = _compact(name, result), False
        except Exception as exc:                       # the model sees the error and can correct itself
            text, error = f"{type(exc).__name__}: {exc}", True
        self.session.log(f"{name} {json.dumps(args)[:200]} -> {'error: ' if error else ''}{text[:200]}", "ai")
        self.on_message("tool", f"{name}({json.dumps(args)[:200]}) -> {text[:300]}")
        return text, error

    def stop(self):
        self.stop_requested = True

    # ---- conversation --------------------------------------------------------------------------
    def ask(self, prompt: str, background: bool = False):
        if self.running:
            raise RuntimeError("the assistant is still working on the previous request")
        self.stop_requested = False
        self.on_message("user", prompt)
        self.session.log(f"[{self.provider}:{self.model}] {prompt}", "ai")
        if background:
            threading.Thread(target=self._safe_turn, args=(prompt,), daemon=True).start()
            return None
        return self._safe_turn(prompt)

    def _safe_turn(self, prompt):
        self.running = True
        try:
            text = self._turn_anthropic(prompt) if self.provider == "anthropic" else self._turn_openai(prompt)
            self.on_message("assistant", text)
            self.session.log(f"answer: {text[:300]}", "ai")
            return text
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            self.on_message("error", msg)
            self.session.log(f"error: {msg}", "ai")
            return msg
        finally:
            self.running = False

    def _turn_anthropic(self, prompt: str) -> str:
        client = self.client()
        tools = self.api.tools("anthropic")
        self.history.append({"role": "user", "content": prompt})
        for _ in range(self.max_steps):
            if self.stop_requested:
                return "(stopped)"
            response = client.beta.messages.create(
                model=self.model, max_tokens=16000, system=SYSTEM, tools=tools, messages=self.history,
                betas=["server-side-fallback-2026-07-01"], fallbacks="default")
            self.history.append({"role": "assistant", "content": response.content})
            if response.stop_reason == "refusal":
                return "(the model declined this request)"
            calls = [b for b in response.content if getattr(b, "type", None) == "tool_use"]
            if not calls:
                if response.stop_reason == "max_tokens":
                    return "".join(b.text for b in response.content if b.type == "text") + " [output cut at max_tokens]"
                return "".join(b.text for b in response.content if b.type == "text")
            results = []
            for b in calls:
                text, error = self._run_tool(b.name, b.input)
                results.append({"type": "tool_result", "tool_use_id": b.id, "content": text, "is_error": error})
            self.history.append({"role": "user", "content": results})
        return f"(stopped after {self.max_steps} steps)"

    def _turn_openai(self, prompt: str) -> str:
        client = self.client()
        tools = [{"type": "function", **t} for t in self.api.tools()]      # name, description, parameters
        self.history.append({"role": "user", "content": prompt})
        for _ in range(self.max_steps):
            if self.stop_requested:
                return "(stopped)"
            response = client.responses.create(model=self.model, instructions=SYSTEM, input=self.history,
                                               tools=tools)
            self.history += list(response.output)
            calls = [o for o in response.output if getattr(o, "type", None) == "function_call"]
            if not calls:
                return getattr(response, "output_text", "") or ""
            for c in calls:
                text, _ = self._run_tool(c.name, c.arguments)
                self.history.append({"type": "function_call_output", "call_id": c.call_id, "output": text})
        return f"(stopped after {self.max_steps} steps)"

    def reset(self):
        self.history = []
