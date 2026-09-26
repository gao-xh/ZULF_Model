import asyncio
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


class AgentToolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        os.environ["ZULF_MODEL_WORKSPACE"] = cls.tmp.name
        from zulf_model.agent import registry
        cls.reg = registry()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_exports_are_consistent(self):
        anthropic = self.reg.anthropic_tools()
        openai = self.reg.openai_tools()
        self.assertEqual([t["name"] for t in anthropic], [t["function"]["name"] for t in openai])
        for tool in anthropic:
            self.assertEqual(tool["input_schema"]["type"], "object")
            json.dumps(tool)

    def test_validation_errors(self):
        with self.assertRaises(ValueError):
            self.reg.call("simulate_transitions", {})
        with self.assertRaises(ValueError):
            self.reg.call("simulate_transitions", {"system": {}, "unexpected": 1})
        with self.assertRaises(KeyError):
            self.reg.call("no_such_tool", {})

    def test_simulate_and_render(self):
        system = {"isotopes": ["13C", "1H", "1H", "1H"],
                  "couplings_hz": [[0, 140, 140, 140], [140, 0, 0, 0], [140, 0, 0, 0], [140, 0, 0, 0]]}
        out = self.reg.call("simulate_transitions", {"system": system})
        freqs = sorted(line["frequency_hz"] for line in out["strongest"])
        np.testing.assert_allclose(freqs, [140.0, 280.0], atol=1e-9)
        for mode in ("pure", "continuous"):
            rendered = self.reg.call("render_spectrum", {"interpretation": {"components": [{"system": system}]},
                                                          "mode": mode, "points": 2048})
            self.assertTrue(Path(rendered["arrays"]).exists())

    def test_refine_on_processed_spectrum(self):
        from zulf_core.physics import compute_transitions
        from zulf_core.render.acquisition import evaluate_spectrum, process_record
        from zulf_core.render.phasing import phase_correct
        from zulf_core.spinsystem import SpinSystem
        from zulf_model.render import Acquisition, Renderer
        j = np.zeros((4, 4)); j[0, 1:] = j[1:, 0] = 140.0
        acq = Acquisition(1000.0, 4096, start_sample=20)
        fid = Renderer(acq).synthesize(compute_transitions(SpinSystem(("13C", "1H", "1H", "1H"), j)), 1.0,
                                       gain=np.exp(0.7j), phase_delay_s=0.0004)
        f = np.arange(100.0, 320.0, 0.1)
        phased = phase_correct(evaluate_spectrum(process_record(fid, acq), acq, f), f, 0.7, 0.0004, acq)
        path = Path(self.tmp.name) / "phased.csv"
        np.savetxt(path, np.column_stack([f, phased.real]), delimiter=",", header="frequency_hz,real")
        start = [[0, 138.5, 138.5, 138.5], [138.5, 0, 0, 0], [138.5, 0, 0, 0], [138.5, 0, 0, 0]]
        spectrum = {"path": str(path), "record": acq.to_dict(), "phasing": {"phase0_rad": 0.7, "delay_s": 0.0004}}
        args = {"candidates": [{"components": [{"system": {"isotopes": ["13C", "1H", "1H", "1H"],
                                                            "couplings_hz": start}}]}],
                "spectrum": spectrum, "ranges": [[120.0, 160.0], [260.0, 300.0]],
                "settings": {"starts": 1}}
        out = self.reg.call("refine_candidates", args)
        self.assertTrue(out["observation"]["real_only"])
        self.assertEqual(out["observation"]["lineshape"], "finite_record")
        couplings = np.asarray(out["results"][0]["interpretation"]["components"][0]["system"]["couplings_hz"])
        self.assertLess(abs(couplings[0, 1] - 140.0), 1e-2)
        with self.assertRaises(ValueError):
            self.reg.call("refine_candidates", {"candidates": args["candidates"], "ranges": args["ranges"]})

    def test_fid_tools_and_background_job(self):
        from zulf_core.physics import compute_transitions
        from zulf_model.render import Acquisition, Renderer
        from zulf_core.spinsystem import SpinSystem
        j = np.zeros((4, 4)); j[0, 1:] = j[1:, 0] = 140.0
        acq = Acquisition(1000.0, 2048)
        fid = Renderer(acq).synthesize(compute_transitions(SpinSystem(("13C", "1H", "1H", "1H"), j)), 2.0)
        path = Path(self.tmp.name) / "fid.npy"
        np.save(path, fid)
        source = {"fid_npy": str(path), "sampling_rate_hz": 1000.0}
        diag = self.reg.call("diagnose_fid", {"source": source})
        self.assertIn("candidate_recipes", diag)
        processed = self.reg.call("process_fid", {"source": source, "ranges": [[100, 300]]})
        self.assertGreater(processed["points"], 10)
        job = self.reg.call("submit_job", {"tool": "simulate_transitions", "arguments": {
            "system": {"isotopes": ["13C", "1H"], "couplings_hz": [[0, 140], [140, 0]]}}})
        for _ in range(150):
            status = self.reg.call("get_job", {"job_id": job["job_id"]})
            if status["status"] in ("complete", "failed"):
                break
            time.sleep(0.2)
        self.assertEqual(status["status"], "complete", msg=json.dumps(status))
        self.assertEqual(status["result"]["n_transitions"], 1)
        self.assertTrue(self.reg.call("list_jobs", {})["jobs"])


class CliTests(unittest.TestCase):
    def run_cli(self, *args, stdin=None):
        env = dict(os.environ, ZULF_MODEL_WORKSPACE=tempfile.mkdtemp())
        proc = subprocess.run([sys.executable, "-m", "zulf_model.cli", *args], input=stdin, capture_output=True,
                              text=True, env=env, cwd=REPO, timeout=300)
        return proc.returncode, json.loads(proc.stdout)

    def test_cli_tool_and_errors(self):
        code, out = self.run_cli("tool", "simulate_transitions", "--json",
                                 json.dumps({"system": {"isotopes": ["13C", "1H"], "couplings_hz": [[0, 140], [140, 0]]}}))
        self.assertEqual(code, 0)
        self.assertEqual(out["n_transitions"], 1)
        code, out = self.run_cli("tool", "simulate_transitions", "--json", "{}")
        self.assertEqual(code, 1)
        self.assertEqual(out["status"], "error")
        code, out = self.run_cli("tools", "export", "--format", "openai")
        self.assertEqual(code, 0)
        self.assertEqual(out[0]["type"], "function")


class McpTests(unittest.TestCase):
    def test_stdio_list_and_call(self):
        try:
            from mcp import ClientSession, StdioServerParameters
            from mcp.client.stdio import stdio_client
        except ImportError:
            self.skipTest("mcp not installed")

        async def session():
            params = StdioServerParameters(command=sys.executable, args=["-m", "zulf_model.agent.mcp_server"],
                                           env=dict(os.environ, ZULF_MODEL_WORKSPACE=tempfile.mkdtemp()), cwd=str(REPO))
            async with stdio_client(params) as (read, write):
                async with ClientSession(read, write) as client:
                    await client.initialize()
                    tools = await client.list_tools()
                    result = await client.call_tool("describe_project", {})
                    return [t.name for t in tools.tools], json.loads(result.content[0].text)

        names, described = asyncio.run(asyncio.wait_for(session(), 120))
        self.assertIn("refine_candidates", names)
        self.assertIn("rules", described)


if __name__ == "__main__":
    unittest.main()
