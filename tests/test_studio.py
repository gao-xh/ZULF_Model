import json
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

from zulf_studio.api import TOOLS, StudioAPI, serve
from zulf_studio.session import StudioSession

GAMMA_H, GAMMA_C = 42.577478, 10.708395          # Hz / uT


def session(spec=None, tmp=None):
    return StudioSession(spec, workspace=tmp or tempfile.mkdtemp(), log_file=False)


METHYL = {"compound": "methyl", "chain": {"groups": [["C1", "C", 3]], "bonds": []}, "one_bond": {"C1": 136.0}}
METHINE = {"compound": "methine", "chain": {"groups": [["C1", "C", 1]], "bonds": []}, "one_bond": {"C1": 140.0}}


class SessionPhysicsTests(unittest.TestCase):
    def test_methyl_at_zero_field_has_lines_at_j_and_2j(self):
        # Closed form: an XA3 group (13C with three equivalent protons) at zero field has lines at J and 2J only.
        lines = session(METHYL).lines(min_relative=0.01)
        np.testing.assert_allclose(sorted(l["frequency_hz"] for l in lines), [136.0, 272.0], atol=1e-9)

    def test_transverse_field_splits_a_ch_line_at_first_order(self):
        # Closed form (docs/CONVENTIONS.md): a transverse field B splits the 13C-1H line by +-(gH + gC) B / 2.
        s = session(METHINE)
        s.set_field(transverse_nt=100.0)
        near = sorted(l["frequency_hz"] for l in s.lines(min_relative=0.05) if l["frequency_hz"] > 100)
        self.assertEqual(len(near), 2)
        self.assertAlmostEqual((near[1] - near[0]) / 2, (GAMMA_H + GAMMA_C) * 0.1 / 2, delta=0.005)
        self.assertAlmostEqual(np.mean(near), 140.0, delta=0.05)

    def test_couplings_override_and_new_isotopologue(self):
        spec = {"compound": "acetonitrile", "chain": {"groups": [["C1", "C", 3], ["C2", "C", 0]], "bonds": [["C1", "C2"]]},
                "one_bond": {"C1": 136.3}}
        s = session(spec)
        self.assertEqual([c["label"] for c in s.components()], ["13C@C1"])
        s.set_couplings({"J(C2,HC1)": -10.0})                  # a coupling makes the nitrile 13C isotopologue visible
        self.assertEqual([c["label"] for c in s.components()], ["13C@C1", "13C@C2"])
        c2 = [l["frequency_hz"] for l in s.lines(0.01) if l["component"] == "13C@C2"]
        np.testing.assert_allclose(sorted(c2), [10.0, 20.0], atol=1e-9)      # XA3 with J = -10 Hz: |J| and 2|J|
        with self.assertRaises(ValueError):
            s.set_couplings({"not a key": 1.0})

    def test_simulation_peaks_at_the_lines(self):
        s = session(METHYL)
        s.set_view(130.0, 142.0)
        s.set_linewidth(1.0)
        sim = s.simulate(points=4001)
        f, re = np.asarray(sim["f"]), np.asarray(sim["sim_re"])
        self.assertAlmostEqual(f[np.argmax(re)], 136.0, delta=0.01)
        half = re > re.max() / 2                                     # FWHM of a Lorentzian = rate / pi
        self.assertAlmostEqual(f[half].max() - f[half].min(), 1.0 / np.pi, delta=0.01)


class DisplayScaleTests(unittest.TestCase):
    def test_scale_is_continuous_while_a_coupling_moves_and_can_be_locked(self):
        # Data with dispersive lines (phase 90 deg): a least-squares scale on the real part crosses zero as the
        # simulated line moves through the dispersion; the magnitude scale must change smoothly and stay > 0.
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            f = np.arange(120.0, 290.0, 0.01)
            g = 0.3
            y = sum(1j * a * g / (g + 1j * (f - f0)) for a, f0 in ((1.0, 136.3), (0.8, 272.6)))
            np.save(tmp / "f.npy", f)
            np.save(tmp / "y.npy", y)
            json.dump([{"id": "x", "x": 1.0, "freq": str(tmp / "f.npy"), "values": str(tmp / "y.npy"),
                        "ranges": [[125, 150], [255, 290]]}], open(tmp / "series.json", "w"))
            s = session(METHYL, d)
            s.load_spectrum(series=str(tmp / "series.json"))
            scales = []
            for j in np.arange(135.5, 137.5, 0.02):
                s.set_couplings({"J(C1,HC1)": float(j)})
                scales.append(s.simulate(points=2000)["scale"])
            scales = np.asarray(scales)
            self.assertTrue(np.all(scales > 0))
            self.assertLess(np.max(np.abs(np.diff(np.log(scales)))), np.log(1.2))
            locked = s.lock_scale(True)["scale"]
            s.set_couplings({"J(C1,HC1)": 140.0})
            self.assertEqual(s.simulate()["scale"], locked)
            s.lock_scale(False)
            self.assertNotEqual(s.simulate()["scale"], locked)


class AutoPhaseTests(unittest.TestCase):
    def test_known_phase_and_delay_are_recovered(self):
        # Data = simulation x exp(-i (phi + 2 pi f tau)) (+ noise): the display must apply +phi, +tau.
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            s = session(METHYL, d)
            s.set_linewidth(1.0)
            f = np.arange(120.0, 290.0, 0.01)
            g = 1.0 / (2 * np.pi)
            sim = sum(a * g / (g + 1j * (f - f0)) for a, f0 in ((1.0, 136.0), (0.8, 272.0)))
            phi, tau = np.radians(40.0), 1.2e-3
            rng = np.random.default_rng(1)
            data = sim * np.exp(-1j * (phi + 2 * np.pi * f * tau)) + 0.002 * (rng.normal(size=len(f)) +
                                                                            1j * rng.normal(size=len(f)))
            np.save(tmp / "f.npy", f)
            np.save(tmp / "y.npy", data)
            json.dump([{"id": "x", "x": 1.0, "freq": str(tmp / "f.npy"), "values": str(tmp / "y.npy"),
                        "ranges": [[125, 150], [255, 290]]}], open(tmp / "s.json", "w"))
            s.load_spectrum(series=str(tmp / "s.json"))
            r = s.auto_phase("model")
            self.assertAlmostEqual(r["phase_deg"], 40.0, delta=0.5)
            self.assertAlmostEqual(r["delay_ms"], 1.2, delta=0.01)
            self.assertGreater(r["match"], 0.9)
            s.set_display(phase_deg=0.0, delay_ms=1.2)
            r = s.auto_phase("data", fit_delay=False)                      # model-free, delay held at the truth
            self.assertAlmostEqual((r["phase_deg"] - 40.0 + 90) % 180 - 90, 0.0, delta=2.0)   # modulo pi
            self.assertEqual(r["delay_ms"], 1.2)


class FitPlumbingTests(unittest.TestCase):
    def test_fit_command_needs_a_series_and_starts_away_from_zero_field(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            s = session(METHYL, d)
            with self.assertRaises(ValueError):
                s.fit_command()
            f = np.linspace(120, 290, 200)
            np.save(tmp / "f.npy", f)
            np.save(tmp / "y.npy", np.zeros(len(f), complex))
            json.dump([{"id": "x", "x": 1.0, "freq": str(tmp / "f.npy"), "values": str(tmp / "y.npy"),
                        "ranges": [[125, 150], [255, 290]]}], open(tmp / "series.json", "w"))
            s.load_spectrum(series=str(tmp / "series.json"))
            argv = s.fit_command(starts=3)
            self.assertEqual(argv[argv.index("--range") + 1], "125.0,290.0")
            self.assertEqual(argv[argv.index("--field-start") + 1], "0.02,0.02")      # zero field: nonzero start
            self.assertEqual(json.loads(argv[argv.index("--couplings") + 1]), {"J(C1,HC1)": 136.0})
            s.set_field(37.0, 45.0)
            argv = s.fit_command(fit_field=False)
            self.assertNotIn("--fit-field", argv)

    def test_apply_fit_and_trace_frames(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"
            run.mkdir()
            json.dump({"couplings": {"J(C1,HC1)": {"J_at_x": [135.5]}}, "scores": [0.1],
                       "spectrum_parameters": {"x": {"field_perp_ut": 0.03, "field_z_ut": 0.05,
                                                     "c0.log_rate0": np.log(2.0), "c0.log_rate1": np.log(4.0),
                                                     "phase_delay": -0.004}}}, open(run / "fit.json", "w"))
            f = np.linspace(120, 150, 50)
            np.savez(run / "trace.npz", f0=f, y0=np.zeros(50, complex), model0=np.zeros((2, 50), complex))
            json.dump({"evaluations": 7, "frames": [{"evaluation": 0, "stage": "fit", "objective": 1.0,
                                                     "J": {"J(C1,HC1)": [136.0]}},
                                                    {"evaluation": 7, "stage": "final", "objective": 0.1,
                                                     "J": {"J(C1,HC1)": [135.5]}}]}, open(run / "trace.json", "w"))
            s = session(METHYL, d)
            r = s.apply_fit(str(run))
            self.assertEqual(r["couplings"], {"J(C1,HC1)": 135.5})
            np.testing.assert_allclose(s.field_nt, [30.0, 50.0])        # the pre-rename field name is read too
            self.assertAlmostEqual(s.rate_per_s, 3.0)                     # median of the family rates
            self.assertEqual(s.trace_index, 1)
            fr = s.trace_frame(0, apply=True)
            self.assertEqual(fr["objective"], 1.0)
            self.assertAlmostEqual(s.couplings()[0]["value"], 136.0)


class ApiTests(unittest.TestCase):
    def test_http_tools_state_and_errors(self):
        s = session(METHYL)
        events = []
        s.listeners.append(events.append)
        server = serve(s, port=0)
        base = f"http://127.0.0.1:{server.server_address[1]}"

        def post(path, body):
            req = urllib.request.Request(base + path, data=json.dumps(body).encode(), method="POST")
            return json.loads(urllib.request.urlopen(req).read())
        try:
            tools = json.loads(urllib.request.urlopen(base + "/api/tools?format=anthropic").read())
            self.assertEqual({t["name"] for t in tools}, set(TOOLS))
            self.assertTrue(all("input_schema" in t for t in tools))
            post("/api/set_field", {"transverse_nt": 12.5, "z_nt": 40})
            post("/api/call", {"tool": "set_couplings", "args": {"values": {"J(C1,HC1)": 130.0}}})
            state = json.loads(urllib.request.urlopen(base + "/api/state").read())
            self.assertEqual(state["field_nt"]["transverse"], 12.5)
            self.assertEqual(state["couplings"][0]["value"], 130.0)
            self.assertIn("field", events)
            sim = post("/api/simulate", {"points": 200})["result"]
            self.assertEqual(len(sim["f"]), 200)
            with self.assertRaises(urllib.error.HTTPError):
                post("/api/no_such_tool", {})
            self.assertTrue(any(e["source"] == "api" for e in s.read_log()))
        finally:
            server.shutdown()

    def test_every_tool_maps_to_a_session_method(self):
        api = StudioAPI(session(METHYL))
        for name, (method, _, schema) in TOOLS.items():
            self.assertTrue(callable(getattr(api.session, method)), name)
            self.assertEqual(schema["type"], "object")


class WindowSmokeTests(unittest.TestCase):
    def test_window_and_session_follow_each_other(self):
        try:
            import os
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            from PySide6.QtWidgets import QApplication
            from zulf_studio.app import StudioWindow
        except ImportError:
            self.skipTest("PySide6 not installed")
        app = QApplication.instance() or QApplication([])
        s = session(METHYL)
        w = StudioWindow(s)
        app.processEvents()
        w.b_z.spin.setValue(25.0)                                     # GUI -> session
        self.assertEqual(s.field_nt[1], 25.0)
        done = threading.Event()
        threading.Thread(target=lambda: (s.set_couplings({"J(C1,HC1)": 131.0}), done.set())).start()
        done.wait(5)
        for _ in range(20):                                           # session (another thread) -> GUI
            app.processEvents()
        self.assertAlmostEqual(w.coupling_rows["J(C1,HC1)"].value(), 131.0)
        w.redraw()
        self.assertGreater(w.lines_table.rowCount(), 0)
        w.close()


if __name__ == "__main__":
    unittest.main()


def _synthetic_fid_2khz(n=16000, j=136.0, fs=2000.0):
    """A methyl-like FID at 2 kHz: saturated plateau to 3.5 ms, then lines at J and 2J (decaying cosines)."""
    t = np.arange(n) / fs
    fid = 40.0 * np.exp(-1.0 * t) * (np.cos(2 * np.pi * j * t) + 0.8 * np.cos(2 * np.pi * 2 * j * t))
    fid += np.random.default_rng(0).normal(0, 0.05, n)
    fid[:7] = 3000.0
    return fid


class FigureTests(unittest.TestCase):
    def test_display_spectrum_uses_the_sampling_rate(self):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
        import paper_figure
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "fid.npy"
            np.save(path, _synthetic_fid_2khz())
            f, spec, _, _, _ = paper_figure.display_spectrum(path, 0.1, 0.1, 2, fs=2000.0, phase0_rad=0.0)
            band = (f > 100) & (f < 200)
            self.assertAlmostEqual(f[band][np.argmax(np.abs(spec[band]))], 136.0, delta=0.05)

    def test_studio_figure_of_manual_parameters_on_a_2khz_series(self):
        import subprocess
        import sys
        import time
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            np.save(tmp / "fid.npy", _synthetic_fid_2khz())
            subprocess.run([sys.executable, str(root / "scripts" / "make_series_entry.py"), "--fid", str(tmp / "fid.npy"),
                            "--id", "synthetic", "--out", str(tmp / "series"), "--sampling-rate", "2000",
                            "--ranges", "125,150;260,285"], check=True, capture_output=True, cwd=root)
            entry = json.loads((tmp / "series" / "series.json").read_text())[0]
            self.assertEqual(Path(entry["source_fid"]), (tmp / "fid.npy").resolve())
            s = session(METHYL, d)
            s.load_spectrum(series=str(tmp / "series" / "series.json"))
            self.assertTrue(s.figure_command()["manual"])               # no applied fit: manual parameters
            s.make_figure(formats="png,svg", dpi=80, title="synthetic")
            t0 = time.time()
            while s.figure_status()["running"] and time.time() - t0 < 120:
                time.sleep(0.2)
            st = s.figure_status()
            self.assertEqual(st["returncode"], 0, s.read_log(30))
            self.assertEqual(st["files"], ["figure.caption.txt", "figure.png", "figure.svg"])
            caption = (Path(st["directory"]) / "figure.caption.txt").read_text()
            self.assertIn("sampling 2000 Hz", caption)
            self.assertIn("Manual parameters (not a fit)", caption)
            out = s.export_figure(str(tmp / "export"))
            self.assertIn("parameters_fit.json", out["files"])


class AssistantTests(unittest.TestCase):
    """The tool loop with fake provider clients (no network): the model's tool calls reach the session."""

    @staticmethod
    def _anthropic_client(script):
        from types import SimpleNamespace as NS
        calls = iter(script)

        class Messages:
            def __init__(self):
                self.requests = []

            def create(self, **kw):
                self.requests.append({**kw, "messages": list(kw["messages"])})      # snapshot of the history
                kind, payload = next(calls)
                if kind == "tool":
                    name, args = payload
                    return NS(stop_reason="tool_use", content=[NS(type="tool_use", id=f"t{len(self.requests)}",
                                                                  name=name, input=args)])
                return NS(stop_reason="end_turn", content=[NS(type="text", text=payload)])
        messages = Messages()
        return NS(beta=NS(messages=messages)), messages

    def test_anthropic_loop_runs_tools_on_the_session(self):
        from zulf_studio.assistant import StudioAssistant
        s = session(METHYL)
        client, messages = self._anthropic_client([("tool", ("set_field", {"transverse_nt": 20, "z_nt": 30})),
                                                   ("tool", ("no_such_tool", {})),
                                                   ("text", "Field set to 20 / 30 nT.")])
        seen = []
        a = StudioAssistant(s, "anthropic", client=client, on_message=lambda r, t: seen.append(r))
        self.assertEqual(a.ask("set the field"), "Field set to 20 / 30 nT.")
        self.assertEqual(s.field_nt, [20.0, 30.0])
        req = messages.requests[0]
        self.assertEqual(req["model"], "claude-opus-5-5")
        self.assertEqual({t["name"] for t in req["tools"]}, set(TOOLS))
        bad = messages.requests[2]["messages"][-1]["content"][0]            # the unknown tool came back as an error
        self.assertTrue(bad["is_error"])
        self.assertEqual(seen, ["user", "tool", "tool", "assistant"])
        self.assertTrue(any(e["source"] == "ai" for e in s.read_log()))

    def test_openai_loop_and_simulate_summary(self):
        from types import SimpleNamespace as NS
        from zulf_studio.assistant import StudioAssistant
        s = session(METHYL)
        outputs = iter([
            [NS(type="function_call", name="set_couplings", arguments=json.dumps({"values": {"J(C1,HC1)": 131.0}}),
                call_id="c1")],
            [NS(type="function_call", name="simulate", arguments="{}", call_id="c2")],
            [NS(type="message")],
        ])
        requests = []

        def create(**kw):
            requests.append({**kw, "input": list(kw["input"])})
            out = next(outputs)
            return NS(output=out, output_text="done" if out[0].type == "message" else "")
        client = NS(responses=NS(create=create))
        with self.assertRaises(ValueError):
            StudioAssistant(s, "openai", model="", client=client)            # OpenAI needs a model name
        a = StudioAssistant(s, "openai", model="test-model", client=client)
        self.assertEqual(a.ask("set J"), "done")
        self.assertEqual(s.couplings()[0]["value"], 131.0)
        self.assertEqual(requests[0]["tools"][0]["type"], "function")
        sim_out = next(i for i in requests[2]["input"] if isinstance(i, dict) and i.get("call_id") == "c2")["output"]
        self.assertNotIn("sim_re", sim_out)                                 # arrays summarised for the model
        self.assertIn("points", json.loads(sim_out))


class AssistantSetupTests(unittest.TestCase):
    def test_session_key_is_used_and_never_logged(self):
        try:
            import os
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            from PySide6.QtWidgets import QApplication
            from zulf_studio.app import StudioWindow
        except ImportError:
            self.skipTest("PySide6 not installed")
        import os
        from unittest import mock
        from zulf_studio import credentials
        app = QApplication.instance() or QApplication([])
        saved = {k: os.environ.pop(k, None) for k in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY")}
        stored = {}
        patches = [mock.patch.object(credentials, "available", lambda: True),               # never the real Keychain
                   mock.patch.object(credentials, "save", lambda v, k: stored.update({v: k}) or True),
                   mock.patch.object(credentials, "delete", lambda v: stored.pop(v, None) is not None)]
        for p in patches:
            p.start()
        try:
            s = session(METHYL)
            w = StudioWindow(s)
            app.processEvents()
            secret = "sk-ant-test-0123456789"
            w.ai_key.setText(secret)
            w._ai_use_key()
            self.assertEqual(os.environ.get("ANTHROPIC_API_KEY"), secret)
            self.assertEqual(w.ai_key.text(), "")
            self.assertIn("credentials found", w.ai_creds.text())
            self.assertFalse(any(secret in e["message"] for e in s.read_log(1000)))
            w._ai_forget_key()
            self.assertNotIn("ANTHROPIC_API_KEY", os.environ)
            self.assertEqual(stored, {})                                   # not remembered: Keychain untouched
            w.ai_key.setText(secret)
            w.ai_remember.setChecked(True)
            w._ai_use_key()
            self.assertEqual(stored, {"ANTHROPIC_API_KEY": secret})
            self.assertFalse(any(secret in e["message"] for e in s.read_log(1000)))
            w._ai_forget_key()
            self.assertEqual(stored, {})
            w.close()
        finally:
            for p in patches:
                p.stop()
            for k, v in saved.items():
                if v is not None:
                    os.environ[k] = v


class CredentialsTests(unittest.TestCase):
    def test_keychain_commands_and_environment(self):
        from types import SimpleNamespace as NS
        from zulf_studio import credentials
        calls, store = [], {"OPENAI_API_KEY": "sk-openai-x"}

        def runner(args, **kw):
            calls.append(args)
            account = args[args.index("-a") + 1]
            if args[1] == "find-generic-password":
                return NS(returncode=0 if account in store else 44, stdout=store.get(account, "") + "\n")
            if args[1] == "add-generic-password":
                store[account] = args[args.index("-w") + 1]
            if args[1] == "delete-generic-password":
                return NS(returncode=0 if store.pop(account, None) else 44, stdout="")
            return NS(returncode=0, stdout="")
        self.assertTrue(credentials.save("ANTHROPIC_API_KEY", "sk-ant-y", runner))
        self.assertEqual(calls[-1][:2], ["security", "add-generic-password"])
        self.assertIn("-U", calls[-1])                                     # replaces an existing item
        self.assertEqual(calls[-1][calls[-1].index("-s") + 1], "zulf-studio")
        env = {"OPENAI_API_KEY": "from-shell"}
        loaded = credentials.load_into_environment(runner, env)
        self.assertEqual(loaded, ["ANTHROPIC_API_KEY"])                    # the shell's key wins
        self.assertEqual(env, {"OPENAI_API_KEY": "from-shell", "ANTHROPIC_API_KEY": "sk-ant-y"})
        self.assertTrue(credentials.delete("ANTHROPIC_API_KEY", runner))
        self.assertIsNone(credentials.load("ANTHROPIC_API_KEY", runner))
        with self.assertRaises(ValueError):
            credentials.save("HOME", "x", runner)


class FilesAndJobsTests(unittest.TestCase):
    """Session files, opening by content, imports and the job list (2026-10-08 UI work)."""

    def test_session_file_round_trip_ignores_unknown_keys_and_reports_missing_data(self):
        a = session(METHYL)
        a.set_couplings({"J(C1,HC1)": 131.5})
        a.set_field(transverse_nt=7.0, z_nt=52.0)
        a.set_linewidth(2.5)
        path = a.save_session(Path(a.workspace) / "methyl")["path"]
        self.assertTrue(path.endswith(".zulfstudio"))
        d = json.loads(Path(path).read_text())
        d["future_key"] = 1                                           # written by a newer Studio
        d["data"] = {"series": str(Path(a.workspace) / "gone" / "series.json"), "index": 0}
        Path(path).write_text(json.dumps(d))
        b = session()
        r = b.open_path(path)                                         # dispatched by content
        self.assertEqual(r["opened"], "session")
        self.assertEqual(len(r["missing"]), 1)
        self.assertAlmostEqual({c["key"]: c["value"] for c in b.couplings()}["J(C1,HC1)"], 131.5)
        self.assertEqual(b.field_nt, [7.0, 52.0])
        self.assertAlmostEqual(b.rate_per_s, 2.5)
        with self.assertRaises(ValueError):
            b.open_path(str(Path(a.workspace) / "methyl.zulfstudio.txt"))

    def test_import_fid_makes_and_loads_a_whole_grid_series(self):
        import time
        tmp = Path(tempfile.mkdtemp())
        np.save(tmp / "average_fid.npy", _synthetic_fid_2khz())
        (tmp / "scans.json").write_text(json.dumps({"sampling_rate_hz": 2000.0}))     # rate found next to the FID
        s = session(METHYL)
        r = s.open_path(str(tmp / "average_fid.npy"))
        self.assertEqual(r["opened"], "FID")
        t0 = time.time()
        while any(job.running for job in s.jobs) and time.time() - t0 < 120:
            time.sleep(0.2)
        time.sleep(0.5)
        self.assertIsNotNone(s.data)
        ranges = s.data["ranges"]
        self.assertAlmostEqual(ranges[0][0], 20.0)
        self.assertAlmostEqual(ranges[-1][1], 380.0)
        self.assertTrue(all(b[0] - a[1] < 1.0 for a, b in zip(ranges, ranges[1:])))   # only small gaps (mains)
        jobs = s.job_list()
        self.assertEqual(jobs[0]["kind"], "import")
        self.assertEqual(jobs[0]["returncode"], 0)
        self.assertEqual(s.blind_command(workers=2)[2], str((tmp / "average_fid.npy").resolve()))  # source_fid

    def test_blind_command_needs_a_fid(self):
        s = session(METHYL)
        with self.assertRaises(ValueError):
            s.blind_command()
        fid = Path(tempfile.mkdtemp()) / "x.npy"
        np.save(fid, np.zeros(10))
        argv = s.blind_command(fid=str(fid), workers=3, structure=METHYL, labeling="15N")
        self.assertIn("analyze_sample.py", argv[1])
        self.assertEqual(argv[argv.index("--workers") + 1], "3")
        self.assertEqual(json.loads(argv[argv.index("--structure") + 1]), METHYL)
        self.assertEqual(argv[argv.index("--labeling") + 1], "15N")

    def test_analysis_processes_counts_runs_not_their_forked_workers(self):
        from zulf_studio.session import analysis_processes
        py = "/opt/env/bin/python"
        ps = "\n".join([
            f"100 1 {py} scripts/fit_joint_series.py --series s.json --starts 8 --workers 3 --out runs/a",
            f"101 100 {py} scripts/fit_joint_series.py --series s.json --starts 8 --workers 3 --out runs/a",
            f"102 100 {py} scripts/fit_joint_series.py --series s.json --starts 8 --workers 3 --out runs/a",
            f"200 1 {py} scripts/fit_joint_series.py --series s.json --starts 2 --workers 6 --out runs/b",
            f"300 1 {py} scripts/analyze_sample.py fid.npy --id x",
            "400 1 /bin/zsh -c python scripts/fit_joint_series.py --workers 9",          # a shell, not a run
            f"500 1 {py} -c import x"])
        procs = analysis_processes(ps)
        self.assertEqual([p["pid"] for p in procs], [100, 200, 300])
        self.assertEqual([p["workers"] for p in procs], [3, 2, 4])     # min(workers, starts); script default
        self.assertEqual(procs[0]["out"], "runs/a")

    def test_machine_status_and_window_job_widgets(self):
        try:
            import os
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            from PySide6.QtWidgets import QApplication
            from zulf_studio.app import StudioWindow
        except ImportError:
            self.skipTest("PySide6 not installed")
        s = session(METHYL)
        m = s.machine_status()
        self.assertGreaterEqual(m["cores"]["logical"], 1)
        self.assertGreaterEqual(m["suggested_workers"], 1)
        app = QApplication.instance() or QApplication([])
        w = StudioWindow(s)
        w.update_jobs()
        self.assertIn("idle", w.activity.text.text())
        self.assertIn("cores", w.machine.text())
        self.assertIn("Jobs", [w.tabs.tabText(i) for i in range(w.tabs.count())])
        s.set_field(z_nt=10.0)
        app.processEvents()
        self.assertIn("\u2022", w.windowTitle())                      # unsaved change
        s.save_session(Path(s.workspace) / "m")
        app.processEvents()
        self.assertNotIn("\u2022", w.windowTitle())
        w.close()
