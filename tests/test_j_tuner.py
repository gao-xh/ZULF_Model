import copy
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


def _argv(tmp, extra=()):
    return ["--series", str(tmp / "series.json"), "--real-only", "false", "--shape", "free", "--exchange", "fast",
            "--range", "110,260", "--structure", json.dumps(STRUCTURE), "--signal-threshold", "2.5", *extra]


def _args(tmp, extra=()):
    return fj.make_parser().parse_args(_argv(tmp, extra))


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

    def test_wrong_basin_is_found_with_held_gains(self):
        prob = self.prob
        s = jt.TuningSession(prob)
        s.z = self.truth_z.copy()
        s.set_couplings({"J(C1,HC2)": 6.0})
        args = _args(self.tmp, ["--component-search-starts", "12", "--seed", "3", "--component-search-hold-gains"])
        args.workers = 1
        z_new, rec = fj.component_search(prob.joint, prob, s.z.copy(), args, "start")
        self.assertLess(rec["objective_after"], 0.5 * rec["objective_before"])
        k = prob.joint.coupling.index(next(n for n in prob.joint.coupling if prob.key_of.get(n, n) == "J(C1,HC2)"))
        self.assertAlmostEqual(float(prob.joint.coupling_values(z_new, k)[0]), TRUTH["J(C1,HC2)"], delta=0.15)


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


class HeldGainLocalFitTests(unittest.TestCase):
    def test_held_gain_cost_is_the_window_part_of_the_full_residual(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, z_truth = _problem(tmp)
            j = prob.joint
            f = j.forwards[0]
            s = jt.TuningSession(prob)
            s.z = z_truth.copy()
            s.set_couplings({"J(C1,HC2)": -3.6})                 # truth -4.2 Hz
            z_bad = s.z.copy()
            labels = prob.model.component_labels
            c1 = labels.index("13C@C1")
            window = next((lo, hi) for c, lo, hi in fj._component_windows(j, z_bad, 0) if c == c1)
            name = next(n for n in j.coupling if prob.key_of.get(n, n) == "J(C1,HC2)")
            out = j.local_fit(z_bad, 0, window, [name], lower=prob.lower, upper=prob.upper, max_nfev=60,
                              hold_gains=True)
            cost, z_new = out[0]
            k = j.coupling.index(name)
            # the gains stay at the start's global solution (fitted with the wrong coupling), which biases the window
            # optimum slightly; the global refit that follows a local candidate removes it
            self.assertAlmostEqual(float(j.coupling_values(z_new, k)[0]), TRUTH["J(C1,HC2)"], delta=0.1)
            # reference: the full forward (every line rendered) with the gains of z_bad held, on the window rows
            x_bad = j.spectrum_vector(z_bad, 0)
            ref = f.predict(x_bad)
            r = np.asarray(f.predict(j.spectrum_vector(z_new, 0), fixed_gains=ref.gains,
                                     fixed_background=ref.background).residual)
            idx = np.flatnonzero((f.f >= window[0]) & (f.f <= window[1]))
            rows = np.concatenate([idx, idx + len(f.f)]) if len(r) == 2 * len(f.f) else idx
            self.assertAlmostEqual(cost, float(r[rows] @ r[rows]), delta=1e-3 * max(cost, 1e-12) + 1e-12)
            # at the start point the held-gain residual is the plain residual (the gains are the global solution)
            np.testing.assert_allclose(np.asarray(f.predict(x_bad, fixed_gains=ref.gains,
                                                            fixed_background=ref.background).residual),
                                       np.asarray(ref.residual), rtol=1e-9, atol=1e-12)


class SnapshotTests(unittest.TestCase):
    def test_snapshot_picks_the_best_vector_by_the_plain_objective(self):
        # Two record files: one claims a tiny cost for a worse vector (as a smoothing stage can), the other a
        # large cost for the truth. The snapshot must rescore and keep the truth, and --from-joint must load it.
        import snapshot_fit
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, z_truth = _problem(tmp)
            s = jt.TuningSession(prob)
            s.z = z_truth.copy()
            s.set_couplings({"J(C1,HC2)": -2.0})
            z_bad = s.z.copy()
            run = tmp / "run"
            (run / "monitor").mkdir(parents=True)
            extra = ["--family-edges", "150,200"]
            p_edges = fj.build_problem(_args(tmp, extra))
            zt, zb = p_edges.z0.copy(), p_edges.z0.copy()     # the couplings of the two vectors, default rates
            zt[:p_edges.joint.nt] = z_truth[:prob.joint.nt]
            zb[:p_edges.joint.nt] = z_bad[:prob.joint.nt]
            json.dump({"argv": ["scripts/fit_joint_series.py"] + _argv(tmp, extra + ["--out", str(run)])},
                      open(run / "monitor" / "status.json", "w"))
            with open(run / "monitor" / "start_000.jsonl", "w") as fh:
                fh.write(json.dumps({"event": "start"}) + "\n")
                fh.write(json.dumps({"n": 5, "cost": 1e-6, "best": 1e-6, "label": "smoothing", "z": zb.tolist()}) + "\n")
            with open(run / "monitor" / "start_001.jsonl", "w") as fh:
                fh.write(json.dumps({"n": 9, "cost": 9.0, "best": 9.0, "label": "fit", "z": zt.tolist()}) + "\n")
                fh.write('{"n": 10, "cost": ')                   # a partly written line of a killed run
            fit = snapshot_fit.snapshot(run, tmp / "snap.json")
            self.assertEqual(fit["snapshot"]["record"], "start_001")
            self.assertEqual(fit["family_edges_hz"], [150.0, 200.0])
            self.assertAlmostEqual(fit["couplings"]["J(C1,HC2)"]["J_at_x"][0], TRUTH["J(C1,HC2)"], places=9)
            j = p_edges.joint
            self.assertAlmostEqual(fit["scores"][0], float(np.sum(j.residual(zt) ** 2)), places=12)
            back = fj.build_problem(_args(tmp, extra + ["--from-joint", str(tmp / "snap.json")]))
            np.testing.assert_allclose(back.z0, np.clip(zt, back.lower + 1e-9, back.upper - 1e-9), atol=1e-12)


class AutoFamilyEdgesTests(unittest.TestCase):
    def test_rule_on_hand_computed_cases(self):
        lines = [100.0, 100.4, 101.0, 110.0, 110.5, 130.0]
        # gaps > 1.5 Hz: 101.0 | 110.0 -> 105.5, 110.5 | 130.0 -> 120.25; sharp 100.42 isolates 100.4: 100.2, 100.7
        self.assertEqual(fj.edges_from_lines(lines, [100.42]), [100.2, 100.7, 105.5, 120.25])
        self.assertEqual(fj.edges_from_lines(lines, [140.0]), [105.5, 120.25])       # no model line near the peak
        # two neighbouring sharp lines share the edge between them (duplicates merged)
        self.assertEqual(fj.edges_from_lines([100.0, 100.4, 100.8], [100.4, 100.8]), [100.2, 100.6])
        # a neighbour closer than min_separation_hz is not split off
        self.assertEqual(fj.edges_from_lines([100.0, 100.1, 103.0], [100.1]), [101.55])
        # edges closer than min_spacing_hz merge into their mean
        self.assertEqual(fj.edges_from_lines([100.0, 100.3, 100.7], [100.3, 100.0], min_spacing_hz=0.3), [100.15, 100.5])
        self.assertEqual(fj.edges_from_lines([120.0]), [])

    def test_auto_edges_isolate_sharp_lines(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            _problem(tmp)
            args = _args(tmp, ["--family-edges", "auto"])
            prob = fj.build_problem(args)
            edges = [float(v) for v in args.family_edges.split(",")]
            self.assertEqual(list(prob.joint.params[0].policy.family_edges_hz), edges)
            self.assertTrue(len(edges) >= 2 and edges == sorted(edges))
            first = fj.build_problem(_args(tmp))                 # the one-family problem the edges come from
            got, lines, sharp = fj.auto_family_edges(first.joint, first.z0)
            np.testing.assert_allclose(got, edges, atol=1e-3)
            fam = np.searchsorted(edges, lines, side="right")
            lines = np.asarray(lines)
            for p in sharp:
                i = int(np.argmin(np.abs(lines - p)))
                if abs(lines[i] - p) > 0.5:
                    continue
                same = lines[(fam == fam[i]) & (np.abs(lines - lines[i]) > 0.2)]
                self.assertEqual(len(same), 0, f"line {lines[i]:.2f} under the sharp peak {p:.2f} shares its family")


class RunToolsTests(unittest.TestCase):
    def test_short_run_then_figures_and_snapshot(self):
        # a real (short) fit_joint_series run on the synthetic spectrum; the tools rebuild it from its directory
        import subprocess
        import plot_components
        import plot_runs
        import run_problem
        import snapshot_fit
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            _problem(tmp)
            out = tmp / "run"
            cmd = [sys.executable, str(Path(fj.__file__)), *_argv(tmp, ["--family-edges", "150,200", "--starts", "1",
                   "--max-nfev", "8", "--component-search", "off", "--out", str(out)])]
            subprocess.run(cmd, check=True, capture_output=True, timeout=600,
                           env={**__import__("os").environ, "OMP_NUM_THREADS": "1"})
            fit = json.load(open(out / "fit.json"))
            self.assertEqual(fit["family_edges_hz"], [150.0, 200.0])
            from zulf_core.nuclei import REDUCED_COUPLING_UNIT, reduced_coupling
            pairs = {k: tuple(c["nuclei"]) for k, c in fit["couplings"].items()}
            self.assertEqual(pairs["J(HC1,HC2)"], ("1H", "1H"))           # shared H-H coupling
            self.assertIn("13C", pairs["J(C1,HC1)"])
            for c in fit["couplings"].values():                            # K = J / (h g_A g_B), per node
                self.assertAlmostEqual(c["K_at_x"][0] * REDUCED_COUPLING_UNIT,
                                       reduced_coupling(c["J_at_x"][0], *c["nuclei"]), delta=1e-9 * 1e21)
            run = run_problem.load_run(out)
            j = run.prob.joint
            for k, n in enumerate(j.coupling):
                key = run.prob.key_of.get(n, n)
                self.assertAlmostEqual(float(j.coupling_values(run.prob.z0, k)[0]), fit["couplings"][key]["J_at_x"][0],
                                       places=6)
            self.assertAlmostEqual(float(np.sum(j.residual(run.prob.z0) ** 2)), fit["scores"][0], delta=1e-6)
            rows = plot_runs.plot(tmp / "runs.png", [out, out], zooms=[(180, 210)])
            self.assertEqual(len(rows), 2)
            comps = plot_components.plot(out, tmp / "components.png", band=(120, 260))
            self.assertEqual([c[0] for c in comps], list(run.prob.model.component_labels))
            self.assertTrue((tmp / "runs.png").stat().st_size > 10000 and (tmp / "components.png").stat().st_size > 10000)
            snap = snapshot_fit.snapshot(out, tmp / "snap.json")
            self.assertLessEqual(snap["scores"][0], fit["scores"][0] * (1 + 1e-6))


class FieldOptionTests(unittest.TestCase):
    def test_fitted_field_matches_a_fixed_field_protocol(self):
        # Reference: the zero-field problem with the protocol fixed at the field (MixtureForward(protocol=...)).
        from zulf_core.physics.protocol import Protocol
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            prob, z = _problem(tmp)
            plain = prob.joint
            fitted = fj.build_problem(_args(tmp, ["--fit-field", "--field-start", "0.03,0.05"])).joint
            self.assertEqual([n for n in fitted.local if n.startswith("field_")], ["field_transverse_ut", "field_z_ut"])
            self.assertFalse(any(n.startswith("field_") for n in plain.local))
            zf = np.concatenate([z[:plain.nt], [z[plain.nt + plain.local.index(n)] if n in plain.local else 0.0
                                                for n in fitted.local]])
            i_trans, i_z = (fitted.nt + fitted.local.index(n) for n in ("field_transverse_ut", "field_z_ut"))
            zf[i_trans], zf[i_z] = 0.0, 0.0
            np.testing.assert_allclose(fitted.residual(zf), plain.residual(z), rtol=1e-9, atol=1e-12)
            zf[i_trans], zf[i_z] = 0.03, 0.05
            x = plain.spectrum_vector(z, 0)
            f0 = plain.forwards[0]
            fixed = copy.copy(f0)                        # same data and weights, protocol fixed at the field
            fixed.protocol, fixed._last = Protocol(field_ut=(0.03, 0.0, 0.05)), None
            fixed.cache = type(f0.cache)(16)
            moved = fitted.residual(zf)
            self.assertGreater(np.linalg.norm(moved - plain.residual(z)), 1e-3 * np.linalg.norm(plain.residual(z)))
            xf = fitted.spectrum_vector(zf, 0)
            vf = fitted.params[0].values(xf)
            self.assertEqual(fitted.forwards[0].protocol_for(vf).field_ut, (0.03, 0.0, 0.05))
            np.testing.assert_allclose(fitted.forwards[0].predict(xf).model, fixed.predict(x).model,
                                       rtol=1e-8, atol=1e-12)
            fitted.peak_sources(zf, 0, 130.0)            # field-aware derivatives run


class FieldRenameTests(unittest.TestCase):
    def test_old_transverse_field_name_is_read(self):
        old = {"phase_delay": -0.004, "field_perp_ut": 0.037, "field_z_ut": 0.045}
        new = fj.remap_family_rates(old, None, [])
        self.assertEqual(new, {"phase_delay": -0.004, "field_transverse_ut": 0.037, "field_z_ut": 0.045})
        self.assertEqual(fj.remap_family_rates({"field_transverse_ut": 0.01}, [1.0], [2.0, 3.0]),
                         {"field_transverse_ut": 0.01})


class TunerSpectrumParameterTests(unittest.TestCase):
    def test_spectrum_sliders_use_display_units_and_change_the_model(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            _problem(tmp)
            prob = fj.build_problem(_args(tmp, ["--fit-field", "--field-start", "0.03,0.05"]))
            s = jt.TuningSession(prob)
            params = {p["name"]: p for p in s.spectrum_parameters()}
            self.assertEqual(params["field_transverse_ut"]["unit"], "nT")
            self.assertAlmostEqual(params["field_transverse_ut"]["value"], 30.0, places=6)
            rate = next(p for p in params.values() if p["group"] == "rates")
            self.assertTrue(rate["log"])
            before = s.evaluate([110.0, 260.0], analysis=False)["spectrum"]["model_re"]
            s.set_spectrum({"field_z_ut": 80.0, rate["name"]: 2.5, "phase_delay": 1.0})
            i = lambda n: s._spectrum_index(n)
            self.assertAlmostEqual(s.z[i("field_z_ut")], 0.08)
            self.assertAlmostEqual(np.exp(s.z[i(rate["name"])]), 2.5)
            self.assertAlmostEqual(s.z[i("phase_delay")], 1e-3)
            after = s.evaluate([110.0, 260.0], analysis=False)
            self.assertGreater(np.max(np.abs(np.asarray(after["spectrum"]["model_re"]) - np.asarray(before))), 0)
            shown = {p["name"]: p["value"] for p in after["spectrum_parameters"]}
            self.assertAlmostEqual(shown["field_z_ut"], 80.0)
            s.set_spectrum({"field_z_ut": 1e6})                      # clipped to the bound (1 uT)
            self.assertLessEqual(s.z[i("field_z_ut")], 1.0)
