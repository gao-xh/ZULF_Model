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
            y = sum(1j * a * g / (g - 1j * (f - f0)) for a, f0 in ((1.0, 136.3), (0.8, 272.6)))
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
