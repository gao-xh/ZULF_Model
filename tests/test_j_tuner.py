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


class ComponentSearchTests(unittest.TestCase):
    """fit_joint_series.component_search: a coupling started in a wrong basin is found again on its window."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.tmp = Path(cls.tmpdir.name)
        cls.prob, cls.truth_z = _problem(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_wrong_basin_is_found(self):
        prob = self.prob
        s = jt.TuningSession(prob)
        s.z = self.truth_z.copy()
        s.set_couplings({"J(C1,HC2)": 6.0})                      # truth -4.2 Hz: far outside a local step
        z_bad = s.z.copy()
        windows = fj._component_windows(prob.joint, z_bad, 0)
        self.assertTrue(any(prob.joint.model.component_labels[c] == "13C@C1" for c, _, _ in windows))
        args = _args(self.tmp, ["--component-search-starts", "12", "--seed", "3"])
        args.workers = 1
        z_new, rec = fj.component_search(prob.joint, prob, z_bad, args, "start")
        self.assertLess(rec["objective_after"], 0.5 * rec["objective_before"])
        k = prob.joint.coupling.index(next(n for n in prob.joint.coupling if prob.key_of.get(n, n) == "J(C1,HC2)"))
        self.assertAlmostEqual(float(prob.joint.coupling_values(z_new, k)[0]), TRUTH["J(C1,HC2)"], delta=0.1)
        self.assertTrue(any(w["accepted"] for w in rec["windows"]))


class GroupIndexTests(unittest.TestCase):
    def test_exchanging_groups_found_in_a_symmetric_component(self):
        # AA'BB' ethylenediamine with the N-H protons kept: J(HC1a,HN1) also maps to the mirror pairs (HC2a-HN2),
        # which made the old intersection empty. Reference: give J(C1,HN1) and J(N1,HN1) unique values and read
        # which group carries them in each component's coupling matrix.
        import regression_confirmed as reg
        from fit_processed_spectrum import override_couplings
        from fit_staged import group_index
        from zulf_hypothesis import exchange_variants
        from zulf_hypothesis.builder import build_model
        spec = {"compound": "eda", "motif": "H2N-CH2-CH2-NH2 (AA'BB')", "one_bond": {"C1": 131.0, "N1": 65.0}}
        frag = override_couplings(reg.structure_for(spec), {"J(C1,HN1)": 3.217, "J(N1,HN1)": -71.93})
        model = build_model(exchange_variants(frag, "slow")[0], ranges=[(85, 240)])
        for c, (heavy, value) in enumerate(((0, 3.217), (0, -71.93))):
            label = model.component_labels[c]
            j = np.asarray(model.interpretation.components[c].system.group_couplings())
            expected = int(np.flatnonzero(np.isclose(j[heavy], value))[0])
            self.assertEqual(group_index(model, "HN1", c), expected, label)
            self.assertNotEqual(group_index(model, "HN2", c), expected, label)


class RateRemapTests(unittest.TestCase):
    def test_new_families_take_the_rate_of_the_old_family_at_their_centre(self):
        # old edges 180: families (-180) (180-); new edges 150, 180, 220: centres 149.5, 165, 200, 220.5
        old = {"c0.log_rate0": 0.1, "c0.log_rate1": 0.2, "c1.log_rate0": 1.1, "c1.log_rate1": 1.2, "phase_delay": 3e-3}
        new = fj.remap_family_rates(old, [180.0], [150.0, 180.0, 220.0])
        self.assertEqual(new, {"c0.log_rate0": 0.1, "c0.log_rate1": 0.1, "c0.log_rate2": 0.2, "c0.log_rate3": 0.2,
                               "c1.log_rate0": 1.1, "c1.log_rate1": 1.1, "c1.log_rate2": 1.2, "c1.log_rate3": 1.2,
                               "phase_delay": 3e-3})
        # fewer families: edges 150, 180, 220 -> 200 (centres 199.5, 200.5): old families 2 and 3
        back = fj.remap_family_rates(new, [150.0, 180.0, 220.0], [200.0])
        self.assertEqual({k: v for k, v in back.items() if k.startswith("c0")}, {"c0.log_rate0": 0.2, "c0.log_rate1": 0.2})
        self.assertEqual(fj.remap_family_rates(old, None, [150.0]), old)        # edges not recorded: unchanged
        self.assertEqual(fj.remap_family_rates(old, [180.0], [180.0]), old)     # same edges: unchanged

    def test_from_joint_loads_rates_across_different_edges(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, _ = _problem(tmp)
            labels = prob.model.component_labels
            fit = {"x": [1.0], "couplings": {}, "family_edges_hz": [180.0],
                   "spectrum_parameters": {"synthetic": {f"c{c}.log_rate{k}": 0.3 + c + 0.1 * k
                                                         for c in range(len(labels)) for k in range(2)}}}
            json.dump(fit, open(tmp / "old.json", "w"))
            p = fj.build_problem(_args(tmp, ["--family-edges", "150,180,220", "--from-joint", str(tmp / "old.json")]))
            j = p.joint
            for c in range(len(labels)):
                got = [float(p.z0[j.nt + j.local.index(f"c{c}.log_rate{k}")]) for k in range(4)]
                np.testing.assert_allclose(got, [0.3 + c, 0.3 + c, 0.4 + c, 0.4 + c])
            # and j_tuner.load_fit does the same
            p2 = fj.build_problem(_args(tmp, ["--family-edges", "150,180,220"]))
            jt.load_fit(p2, fit)
            np.testing.assert_allclose(p2.z0[p2.joint.nt:], p.z0[j.nt:])


class TiedRatesTests(unittest.TestCase):
    def test_tied_families_equal_one_family(self):
        # Reference: the same structure without family edges has exactly one rate per component. With edges and
        # every component tied, the residual and the Jacobian (rate column = sum of the family columns) must match.
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, z_truth = _problem(tmp)
            one = fj.build_problem(_args(tmp))
            tied = fj.build_problem(_args(tmp, ["--family-edges", "150,200", "--tie-rates", "."]))
            j1, j2 = one.joint, tied.joint
            self.assertEqual(j1.local, j2.local)                  # followers are not free parameters
            z = z_truth.copy()
            for i, n in enumerate(j1.local):
                if ".log_rate" in n:
                    z[j1.nt + i] = 0.7 + 0.2 * i
            np.testing.assert_allclose(j2.residual(z), j1.residual(z), rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(j2.jacobian(z), j1.jacobian(z), rtol=1e-7, atol=1e-9)
            partly = fj.build_problem(_args(tmp, ["--family-edges", "150,200", "--tie-rates", "C2"]))
            free = [n for n in partly.joint.local if ".log_rate" in n]
            labels = partly.model.component_labels
            c2 = labels.index(next(l for l in labels if "C2" in l))
            self.assertEqual([n for n in free if n.startswith(f"c{c2}.")], [f"c{c2}.log_rate0"])
            self.assertEqual(len([n for n in free if not n.startswith(f"c{c2}.")]), 3)
            with self.assertRaises(ValueError):
                fj.build_problem(_args(tmp, ["--family-edges", "150", "--tie-rates", "no such component"]))
