import json
import sys
import tempfile
import threading
import time
import unittest
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fit_joint_series as fj                                          # noqa: E402
import j_tuner as jt                                                   # noqa: E402

from zulf_core.render.acquisition import Acquisition

STRUCTURE = {"compound": "ethyl test", "chain": {"groups": [["C1", "C", 2], ["C2", "C", 3]], "bonds": [["C1", "C2"]]},
             "one_bond": {"C1": 130.0, "C2": 125.0}}
TRUTH = {"J(C1,HC2)": -4.2, "J(HC1,HC2)": 7.1}


def _args(tmp, extra=()):
    return fj.make_parser().parse_args(
        ["--series", str(tmp / "series.json"), "--real-only", "false", "--shape", "free", "--exchange", "fast",
         "--range", "110,260", "--structure", json.dumps(STRUCTURE), "--signal-threshold", "2.5", *extra])


def _problem(tmp):
    """A synthetic ethyl spectrum: the model at TRUTH (unit gains) plus noise, as a one-entry series file."""
    f = np.arange(110.0, 260.0, 0.05)
    acq = Acquisition.pure(1000.0, 4000)
    np.save(tmp / "f.npy", f)
    np.save(tmp / "y.npy", np.zeros(len(f), complex))
    entry = {"id": "synthetic", "x": 1.0, "freq": str(tmp / "f.npy"), "values": str(tmp / "y.npy"),
             "record": acq.to_dict(), "ranges": [[110.0, 260.0]]}
    json.dump([entry], open(tmp / "series.json", "w"))
    prob = fj.build_problem(_args(tmp))
    s = jt.TuningSession(prob)
    s.set_couplings(TRUTH)
    x = prob.joint.spectrum_vector(s.z, 0)
    clean = prob.joint.forwards[0].predict(x, fixed_gains=np.ones(len(prob.model.component_labels), complex)).model
    rng = np.random.default_rng(5)
    noisy = clean + 0.002 * np.abs(clean).max() * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
    np.save(tmp / "y.npy", noisy)
    return fj.build_problem(_args(tmp)), s.z.copy()


class JTunerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.tmp = Path(cls.tmpdir.name)
        cls.prob, cls.truth_z = _problem(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def session(self):
        prob = fj.build_problem(_args(self.tmp))
        return jt.TuningSession(prob)

    def test_set_and_evaluate_match_the_fitter_objective(self):
        s = self.session()
        before = {c["key"]: c["value"] for c in s.couplings()}
        s.set_couplings({"J(HC1,HC2)": 6.5})
        after = {c["key"]: c["value"] for c in s.couplings()}
        self.assertAlmostEqual(after["J(HC1,HC2)"], 6.5, places=9)
        self.assertEqual({k: v for k, v in after.items() if k != "J(HC1,HC2)"},
                         {k: v for k, v in before.items() if k != "J(HC1,HC2)"})
        view = [180.0, 210.0]
        out = s.evaluate(view)
        j = s.joint
        self.assertAlmostEqual(out["scores"]["objective"], float(np.sum(j.residual(s.z) ** 2)), places=12)
        r = j.forwards[0].predict(j.spectrum_vector(s.z, 0)).residual
        f = j.forwards[0].f
        sel = np.flatnonzero((f >= view[0]) & (f <= view[1]))
        self.assertAlmostEqual(out["scores"]["window"], float(np.sum(r[np.r_[sel, sel + len(f)]] ** 2)), places=12)
        self.assertTrue(all(view[0] <= v <= view[1] for v in out["spectrum"]["f"]))
        self.assertTrue(all(view[0] <= l["frequency_hz"] <= view[1] for l in out["lines"]))
        self.assertIn("residual_peaks", out)

    def test_refinement_recovers_the_truth_and_moves_only_the_chosen_couplings(self):
        s = self.session()
        s.set_couplings({"J(C1,HC2)": TRUTH["J(C1,HC2)"] + 0.4, "J(HC1,HC2)": TRUTH["J(HC1,HC2)"] - 0.4})
        start = s.scores(s.z)[0]["objective"]
        held = {c["key"]: c["value"] for c in s.couplings() if c["key"] not in TRUTH}
        s.refine(list(TRUTH), "none", None, 40)
        while s.progress()["running"]:
            time.sleep(0.05)
        got = {c["key"]: c["value"] for c in s.couplings()}
        for k, v in TRUTH.items():
            self.assertAlmostEqual(got[k], v, delta=0.05, msg=k)
        self.assertEqual({k: got[k] for k in held}, held)
        self.assertLess(s.scores(s.z)[0]["objective"], start)
        s.undo()                                        # back to the state before the refinement
        self.assertAlmostEqual({c["key"]: c["value"] for c in s.couplings()}["J(HC1,HC2)"], TRUTH["J(HC1,HC2)"] - 0.4)

    def test_stop_keeps_the_best_point(self):
        s = self.session()
        s.set_couplings({"J(HC1,HC2)": TRUTH["J(HC1,HC2)"] - 0.5})
        start = s.scores(s.z)[0]["objective"]
        original = s.joint.residual
        calls = [0]

        def counting(z):
            calls[0] += 1
            if calls[0] == 4:
                s.stop()
            return original(z)
        s.joint.residual = counting
        try:
            s.refine(["J(HC1,HC2)"], "none", None, 200)
            while s.progress()["running"]:
                time.sleep(0.05)
        finally:
            del s.joint.residual
        self.assertEqual(s.progress()["message"], "stopped")
        self.assertLessEqual(s.scores(s.z)[0]["objective"], start)

    def test_gradient_against_finite_differences(self):
        s = self.session()
        s.set_couplings({"J(HC1,HC2)": 6.8})
        g = {c["key"]: c["gradient"] for c in s.gradient()}
        h = 1e-4
        for key in ("J(HC1,HC2)", "J(C1,HC2)"):
            v = {c["key"]: c["value"] for c in s.couplings()}[key]
            zp = s.set_couplings({key: v + h}, s.z.copy())
            zm = s.set_couplings({key: v - h}, s.z.copy())
            numeric = (np.sum(s.joint.residual(zp) ** 2) - np.sum(s.joint.residual(zm) ** 2)) / (2 * h)
            self.assertAlmostEqual(g[key], numeric, delta=1e-3 * max(abs(numeric), 1e-3), msg=key)

    def test_saved_result_loads_back_and_reads_as_from_joint(self):
        s = self.session()
        s.set_couplings({"J(HC1,HC2)": 6.9})
        path = s.save(self.tmp / "out" / "tuned.json")
        with open(path) as fh:
            saved = json.load(fh)
        prob = fj.build_problem(_args(self.tmp, ["--from-joint", path]))   # the fitter reads it
        np.testing.assert_allclose(prob.z0, s.z, atol=1e-9)
        prob2 = fj.build_problem(_args(self.tmp))
        self.assertEqual(jt.load_fit(prob2, saved), [])
        np.testing.assert_allclose(prob2.z0, s.z, atol=1e-9)

    def test_http_api(self):
        s = self.session()
        server = jt.serve(s, "127.0.0.1", 0, str(self.tmp / "http"))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        base = f"http://127.0.0.1:{server.server_address[1]}"
        try:
            page = urllib.request.urlopen(base + "/").read().decode()
            self.assertIn("J Tuner", page)
            state = json.loads(urllib.request.urlopen(base + "/api/state").read())
            self.assertEqual([c["key"] for c in state["couplings"]], s.keys)
            req = urllib.request.Request(base + "/api/evaluate", json.dumps(
                {"values": {"J(HC1,HC2)": 6.6}, "view": [180, 210], "commit": True}).encode(),
                {"Content-Type": "application/json"})
            out = json.loads(urllib.request.urlopen(req).read())
            self.assertAlmostEqual({c["key"]: c["value"] for c in out["couplings"]}["J(HC1,HC2)"], 6.6, places=9)
            self.assertGreater(out["scores"]["objective"], 0)
            req = urllib.request.Request(base + "/api/undo", json.dumps({"view": [180, 210]}).encode(),
                                         {"Content-Type": "application/json"})
            out = json.loads(urllib.request.urlopen(req).read())
            self.assertNotAlmostEqual({c["key"]: c["value"] for c in out["couplings"]}["J(HC1,HC2)"], 6.6)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()


class FitMonitorTests(unittest.TestCase):
    """The monitor observes the fit; the optimizer's path and result must be identical with and without it."""

    def test_monitoring_does_not_change_a_start_and_records_it(self):
        import fit_monitor
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, _ = _problem(tmp)
            joint = prob.joint
            s = jt.TuningSession(prob)
            z = s.set_couplings({"J(C1,HC2)": TRUTH["J(C1,HC2)"] + 0.5, "J(HC1,HC2)": TRUTH["J(HC1,HC2)"] - 0.5})
            keys = [prob.key_of.get(n, n) for n in joint.coupling]
            fj._JOINT_TASK = (joint, prob.lower, prob.upper, 30, [], None, False)
            try:
                fj._MONITOR = None
                plain = fj._solve_indexed((0, z.copy()))
                fj._MONITOR = (str(tmp / "monitor"), keys)
                status = fit_monitor.RunStatus(tmp / "monitor", ["fit_joint_series.py"], 1, keys, joint.nodes.tolist(),
                                               {k: [0.0] for k in keys})
                watched = fj._solve_indexed((0, z.copy()))
                status.finished(0, watched[0])
                status.set(phase="finished")
            finally:
                fj._MONITOR = None
                fj._JOINT_TASK = None
            self.assertEqual(plain[0], watched[0])
            np.testing.assert_array_equal(plain[1], watched[1])
            self.assertIsNone(joint.monitor)
            run = fit_monitor.read_run(tmp)
            rec = run["starts"]["start_000"]
            self.assertFalse(rec["running"])
            self.assertGreater(rec["evaluations"], 2)
            self.assertAlmostEqual(rec["best"], float(np.sum(joint.residual(np.asarray(rec["z"])) ** 2)), places=12)
            self.assertLessEqual(rec["best"], rec["history"][0][1])
            self.assertAlmostEqual(rec["J0"]["J(HC1,HC2)"][0], TRUTH["J(HC1,HC2)"] - 0.5, places=9)
            self.assertEqual(run["status"]["phase"], "finished")
            self.assertEqual(fit_monitor.find_runs(tmp), [tmp])


class BandDiagnosisTests(unittest.TestCase):
    """scripts/band_diagnosis.py on the synthetic ethyl spectrum with one coupling detuned."""

    @classmethod
    def setUpClass(cls):
        import band_diagnosis as bd
        cls.bd = bd
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.tmp = Path(cls.tmpdir.name)
        cls.prob, cls.truth_z = _problem(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_detuned_coupling_is_the_free_knob(self):
        prob, bd = self.prob, self.bd
        s = jt.TuningSession(prob)
        s.z = self.truth_z.copy()
        s.set_couplings({"J(C1,HC2)": TRUTH["J(C1,HC2)"] + 0.1})
        z = s.z.copy()
        lines = s.lines(z, (110.0, 260.0), min_relative=0.05)
        c1 = [l["frequency_hz"] for l in lines if l["label"].startswith("13C@C1")]
        lo, hi = min(c1) - 1.0, max(c1) + 1.0
        cost_b, levers = bd.band_levers(prob, z, 0, lo, hi, max_step_hz=1.0, min_gain=0.03)
        top = levers[0]
        self.assertEqual(top["name"], "J(C1,HC2)")
        self.assertEqual(top["verdict"], "free knob")
        self.assertAlmostEqual(top["step_hz"], -0.1, delta=0.025)          # Gauss-Newton: near-linear at 0.1 Hz
        # independent check: the band cost actually evaluated at the stepped vector
        z2 = z.copy()
        z2[top["index"]] += top["step"]
        f = prob.joint.forwards[0]
        r2 = np.asarray(f.predict(prob.joint.spectrum_vector(z2, 0)).residual)
        rows = bd.band_rows(f, lo, hi, len(r2))
        actual_gain = cost_b - float(r2[rows] @ r2[rows])
        self.assertGreater(actual_gain, 0.0)
        self.assertAlmostEqual(top["band_gain"] / actual_gain, 1.0, delta=0.25)
