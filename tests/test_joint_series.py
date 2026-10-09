import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from fit_joint_series import JointSeries                                # noqa: E402

from zulf_core.render.acquisition import Acquisition
from zulf_core.solver import ObservedSpectrum, RefineSettings
from zulf_core.solver.forward import MixtureForward
from zulf_core.solver.refine import _signal_kwargs
from zulf_hypothesis import build_model, template
from zulf_hypothesis.search import _settings_for


class JointSeriesTests(unittest.TestCase):
    def test_monotone_couplings_and_jacobian_against_finite_differences(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        rng = np.random.default_rng(2)
        obs = []
        for shift in (0.0, 0.4, 0.9):
            p = settings.parameterize(model.interpretation)
            name = model.coupling_names["J(Ca,Ha)"][0]
            p.set(name, p.parameters[name].value + shift)
            sig = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                                 **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model
            obs.append(ObservedSpectrum.from_spectrum(f, np.asarray(sig) + 0.01 * rng.normal(size=len(f)),
                                                      [(100.0, 200.0)], record=acq))
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 1.0])
        table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params])
        z = joint.pack(table + np.array([0.0, 0.3, 0.2])[:, None], [p.vector() for p in joint.params])
        z[:joint.nt:joint.m] += 0.05
        z[1:joint.nt:joint.m] = np.where(np.arange(joint.nc) % 2, 0.7, -0.4)    # both directions
        joint.set_priors([0], [1.0 + joint.coupling_values(z, 0).mean()], [0.5], 2.0)
        for k in range(joint.nc):   # monotone by construction, direction = sign of A (independent check)
            steps = np.diff(joint.coupling_values(z, k))
            a = joint.theta(z, k)[1]
            self.assertTrue(np.all(steps * np.sign(a) >= -1e-12))
            self.assertAlmostEqual(joint.coupling_values(z, k)[-1] - joint.coupling_values(z, k)[0], a, places=10)
        for sigma in (0.0, 1.5):         # unsmoothed, and the coarse-to-fine smoothed objective
            joint.set_smoothing(sigma)
            jac = joint.jacobian(z)
            h = 1e-6
            numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                       for e in np.eye(len(z))])
            self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0), msg=str(sigma))
        # smoothing is the same linear kernel on data and model: a constant offset of the residual is kept
        joint.set_smoothing(1.5)
        n_res = len(joint.forwards[0].predict(joint.spectrum_vector(z, 0)).residual)   # 2F for complex data
        self.assertEqual(joint.smoothing[0].shape, (n_res, n_res))

    def test_shared_delay_is_one_value_and_jacobian_matches_finite_differences(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        rng = np.random.default_rng(8)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 200.0)], record=acq) for _ in range(3)]
        free = JointSeries(model, settings, obs, [0.1, 0.5, 1.0])
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 1.0], shared=["phase_delay"])
        self.assertEqual(joint.shared, ["phase_delay"])
        self.assertNotIn("phase_delay", joint.local)
        self.assertEqual(len(free.pack(*self._table(free))) - len(joint.pack(*self._table(joint))), 2)
        table, xs = self._table(joint)
        for s, d in enumerate((0.001, 0.002, 0.006)):
            xs[s][joint.col["phase_delay"]] = d
        z = joint.pack(table + np.array([0.0, 0.3, 0.2])[:, None], xs)
        self.assertAlmostEqual(z[joint.ntheta], 0.003)          # start: the mean of the per-spectrum delays
        for s in range(3):
            self.assertEqual(joint.spectrum_vector(z, s)[joint.col["phase_delay"]], z[joint.ntheta])
        lo, hi = joint.bounds(20.0)
        self.assertEqual(len(lo), len(z))
        jac = joint.jacobian(z)
        h = 1e-7
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_free_shape_takes_any_values_and_jacobian_matches_finite_differences(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        rng = np.random.default_rng(9)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 200.0)], record=acq) for _ in range(3)]
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 1.0], shape="free")
        table, xs = self._table(joint)
        table = table + np.array([0.0, 0.6, -0.3])[:, None]                 # up then down: not monotone
        z = joint.pack(table, xs)
        for k in range(joint.nc):
            np.testing.assert_allclose(joint.coupling_values(z, k), table[:, k])
        lo, hi = joint.bounds(20.0)
        self.assertEqual(len(lo), len(z))
        shift = joint.level_shift(np.arange(len(z), dtype=float))
        for k in range(joint.nc):                                           # one offset per coupling
            self.assertTrue(np.all(shift[k * joint.m:(k + 1) * joint.m] == k * joint.m))
        joint.set_priors([0], [1.0 + joint.coupling_values(z, 0).mean()], [0.5], 2.0)
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_one_spectrum_fits_in_free_shape_only(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        rng = np.random.default_rng(10)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 140.0), (150.0, 200.0)], record=acq)]
        with self.assertRaises(ValueError):
            JointSeries(model, settings, obs, [1.0])
        joint = JointSeries(model, settings, obs, [1.0], shape="free")
        self.assertEqual((joint.nn, joint.m), (1, 1))
        table, xs = self._table(joint)
        z = joint.pack(table + 0.3, xs)
        for k in range(joint.nc):
            np.testing.assert_allclose(joint.coupling_values(z, k), table[:, k] + 0.3)
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_rate_families_from_the_command_line_and_jacobian(self):
        import argparse
        from fit_joint_series import _rate_policy
        base = RefineSettings(band_weighting="none", background_order=-1)
        out = _rate_policy(base, argparse.Namespace(family_edges="120,160", rate_bounds="0.2,8"))
        self.assertEqual(out.policy.family_edges_hz, (120.0, 160.0))
        self.assertEqual(out.policy.rate_bounds_per_s, (0.2, 8.0))
        delayed = _rate_policy(base, argparse.Namespace(family_edges="", rate_bounds="", phase_delay_bounds="-4.6,-2.6"))
        np.testing.assert_allclose(delayed.policy.phase_delay_bounds_s, (-0.0046, -0.0026), rtol=1e-12)
        model = build_model(template("CH-CH3"))
        p = _settings_for(model, "ratios", delayed).parameterize(model.interpretation).parameters["phase_delay"]
        np.testing.assert_allclose((p.lower, p.upper), (-0.0046, -0.0026), rtol=1e-12)
        self.assertAlmostEqual(p.value, -0.0036)                       # the prior excludes 0: start in its middle
        settings = _settings_for(model, "ratios", out)
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.25)
        rng = np.random.default_rng(11)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 300.0)], record=acq)]
        joint = JointSeries(model, settings, obs, [1.0], shape="free")
        rates = [n for n in joint.local if "log_rate" in n]
        self.assertEqual(len(rates), 3 * len(model.component_labels))     # three families per isotopologue
        table, xs = self._table(joint)
        z = joint.pack(table, xs)
        for i, n in enumerate(joint.local):                                # distinct rates per family
            if "log_rate" in n:
                z[joint.nt + i] = np.log(0.5 + 0.7 * i)
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_residual_peaks_found_where_a_line_is_missing_and_rows_have_the_right_jacobian(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.05)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 300.0)], record=acq)
        p = settings.parameterize(model.interpretation)
        sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                        band_weighting="none", **_signal_kwargs(settings)).predict(
            p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model)
        tl = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                            **_signal_kwargs(settings))
        values = p.values(p.vector())
        lines = np.concatenate([np.asarray(tl.transitions(values, c, system).frequencies_hz)
                                for c, system in enumerate(p.systems(values))])
        lines = lines[(lines > 110) & (lines < 290)]
        near = float(lines[0]) + 0.2                    # an extra narrow line next to a model line: assignable
        far = float(np.clip(lines.max() + 15.0, 110, 295))  # one far from every model line: reported only
        rng = np.random.default_rng(12)
        extra = sum(0.3 * np.abs(sig).max() / (1 + ((f - c) / 0.05) ** 2) for c in (near, far))
        dip = float(lines[-1]) - 0.2                    # a negative feature (data below the model) next to a line
        extra = extra - 0.3 * np.abs(sig).max() / (1 + ((f - dip) / 0.05) ** 2)
        broad = 0.05 * np.abs(sig).max() * np.exp(-((f - 200.0) / 15.0) ** 2)     # spread-out misfit
        y = sig + extra + broad + 1e-4 * np.abs(sig).max() * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
        obs = [ObservedSpectrum.from_spectrum(f, y, [(100.0, 300.0)], record=acq)]
        joint = JointSeries(model, settings, obs, [1.0], shape="free")
        table, xs = self._table(joint)
        z = joint.pack(table, xs)
        windows, report = joint.find_residual_peaks(z, k_sigma=6.0)
        found = {round(r["frequency_hz"], 1): r["assignable"] for r in report}
        self.assertTrue(any(abs(k - near) < 0.1 and v for k, v in found.items()), msg=str(report))
        self.assertTrue(any(abs(k - far) < 0.1 and not v for k, v in found.items()), msg=str(report))
        self.assertFalse(any(abs(k - 200.0) < 5.0 for k in found), msg=str(report))       # the broad bump: no peak
        signs = {round(r["frequency_hz"], 1): r["sign"] for r in report}
        self.assertTrue(any(abs(k - dip) < 0.1 and v == "model above data" for k, v in signs.items()), msg=str(report))
        self.assertTrue(any(abs(k - near) < 0.1 and v == "model below data" for k, v in signs.items()), msg=str(report))
        self.assertGreater(len(windows[0]), 0)
        joint.res_windows, joint.res_strength = windows, 3.0
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_dip_rows_penalise_a_broad_line_over_a_splitting_and_have_the_right_jacobian(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.05)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 300.0)], record=acq)
        p = settings.parameterize(model.interpretation)
        sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                        band_weighting="none", **_signal_kwargs(settings)).predict(
            p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model)
        rng = np.random.default_rng(4)
        y = sig + 1e-4 * np.abs(sig).max() * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
        joint = JointSeries(model, settings, [ObservedSpectrum.from_spectrum(f, y, [(100.0, 300.0)], record=acq)],
                            [1.0], shape="free")
        table, xs = self._table(joint)
        z = joint.pack(table, xs)
        rates = [joint.nt + i for i, n in enumerate(joint.local) if "log_rate" in n]
        broad = z.copy()
        broad[rates] += np.log(8.0)                     # every line 8 x broader: valleys filled, tops missed
        joint.set_peak_penalty(0.0, prominence=0.02, tolerance_hz=0.15, min_sigma=3.0, dips=5.0)
        self.assertTrue(np.all(joint.peak_signs[0] < 0) and len(joint.peaks[0]) > 0)
        r_true = joint.forwards[0].predict(joint.spectrum_vector(z, 0)).residual
        r_broad = joint.forwards[0].predict(joint.spectrum_vector(broad, 0)).residual
        self.assertLess(np.abs(joint._peak_rows(0, r_true)[0]).max(), 1e-3)          # the true model reaches every valley
        self.assertGreater(np.abs(joint._peak_rows(0, r_broad)[0]).max(), 0.05)      # the broad one fills them
        joint.set_peak_penalty(5.0, prominence=0.02, tolerance_hz=0.15, min_sigma=3.0, dips=5.0)
        self.assertTrue(np.any(joint.peak_signs[0] > 0) and np.any(joint.peak_signs[0] < 0))
        n_all = len(joint.peaks[0])
        joint.set_peak_penalty(5.0, prominence=0.02, tolerance_hz=0.15, min_sigma=3.0, dips=5.0, max_width_hz=0.01)
        self.assertLess(len(joint.peaks[0]), n_all)               # nothing is narrower than 0.01 Hz on this grid
        joint.set_peak_penalty(5.0, prominence=0.02, tolerance_hz=0.15, min_sigma=3.0, dips=5.0, max_width_hz=2.0)
        zz = broad.copy()
        zz[rates] -= 0.7
        for smooth in (0.0, 0.03):
            joint.peak_smooth = smooth
            jac = joint.jacobian(zz)
            h = 1e-6
            numeric = np.column_stack([(joint.residual(zz + h * e) - joint.residual(zz - h * e)) / (2 * h)
                                       for e in np.eye(len(zz))])
            self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0), msg=str(smooth))

    def test_guard_rows_penalise_model_peaks_where_the_data_are_empty(self):
        # D57: fit range 100-200 Hz; the model also has lines at 249-264 Hz. With those lines in the data the guard
        # is quiet there (data line regions are protected); with empty data there the guard rows are active at
        # those bins only, continuous, and their Jacobian (gains solved again) matches central differences.
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.05)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 300.0)], record=acq)
        p = settings.parameterize(model.interpretation)
        sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                        band_weighting="none", **_signal_kwargs(settings)).predict(
            p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model)
        rng = np.random.default_rng(5)
        noise = 2e-3 * np.abs(sig).max() * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
        empty = np.where(f > 200.0, 0.0, sig) + noise
        for data, active in ((sig + noise, False), (empty, True)):
            obs = ObservedSpectrum.from_spectrum(f, data, [(100.0, 200.0)], record=acq)
            joint = JointSeries(model, settings, [obs], [1.0], shape="free")
            joint.set_guard_penalty(5.0, grids=[(f, data)])
            table, xs = self._table(joint)
            z = joint.pack(table, xs)
            x = joint.spectrum_vector(z, 0)
            g = joint.guard[0]
            self.assertGreater(len(g["groups"]), 50)
            for line in (129.1, 130.8, 146.8):                    # data lines inside the fit range are protected
                self.assertGreater(np.abs(g["f"] - line).min(), 2.0)
            rows = joint._guard_rows(0, x, joint.forwards[0].predict(x).gains)
            centres = np.array([g["f"][idx].mean() for idx in g["groups"]])
            hot = centres[rows > 0.05 * 5.0 / joint.forwards[0].norm]
            if not active:
                self.assertEqual(len(hot), 0, msg=str(hot))
                continue
            self.assertTrue(len(hot) > 0 and np.all((hot > 245.0) & (hot < 268.0)), msg=str(hot))
            jac = joint.jacobian(z)[-len(rows):]
            h = 1e-6
            numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e))[-len(rows):] / (2 * h)
                                       for e in np.eye(len(z))])
            self.assertLess(np.abs(jac - numeric).max(), 1e-3 * np.abs(numeric).max())
            small = joint.residual(z + 1e-9)[-len(rows):] - joint.residual(z)[-len(rows):]
            self.assertLess(np.abs(small).max(), 1e-5 * np.abs(rows).max())               # continuous

    def test_peak_family_edges_give_every_data_peak_its_own_family(self):
        # three lines, the third on the tail of the second (small prominence, tall above the noise): one edge
        # between every pair, none in the noise
        from fit_joint_series import peak_family_edges
        f = np.arange(100.0, 200.0, 0.05)
        y = sum(a / (1 + 1j * (f - c) / w) for a, c, w in ((1.0, 120.0, 0.2), (1.0, 150.0, 0.8), (0.6, 152.5, 0.3)))
        y = y + 0.003 * (np.random.default_rng(1).normal(size=len(f)) + 1j * np.random.default_rng(2).normal(size=len(f)))

        class Fw:
            pass
        fw = Fw()
        fw.f, fw.y, fw.data_signal_mask = f, y, np.abs(y) > 0.05

        class J:
            forwards = [fw]
        edges, peaks = peak_family_edges(J())
        np.testing.assert_allclose(peaks, [120.0, 150.0, 152.5], atol=0.1)
        np.testing.assert_allclose(edges, [135.0, 151.25], atol=0.1)

    def test_line_family_edges_separate_every_resolved_model_line(self):
        from fit_joint_series import edges_between_lines
        edges, groups = edges_between_lines([110.0, 100.1, 100.0, 105.0], min_distance_hz=0.25)
        self.assertEqual(groups, [[100.0, 100.1], [105.0], [110.0]])   # 0.1 Hz apart: one family
        np.testing.assert_allclose(edges, [102.55, 107.5])

    def test_line_band_and_local_fit(self):
        # a line band that holds every line changes nothing; a narrow one keeps the Jacobian exact (finite
        # differences); the local fit recovers a coupling shift from one window
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.05)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 300.0)], record=acq)
        p = settings.parameterize(model.interpretation)
        name = model.coupling_names["J(Ca,Ha)"][0]
        p.set(name, p.parameters[name].value + 0.3)
        fw = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                            **_signal_kwargs(settings))
        sig = np.asarray(fw.predict(p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model)
        obs = ObservedSpectrum.from_spectrum(f, sig, [(100.0, 300.0)], record=acq)
        q = settings.parameterize(model.interpretation)
        full = MixtureForward(q, obs, gain_model=settings.gain_model, background=-1, band_weighting="none",
                              **_signal_kwargs(settings))
        wide = MixtureForward(q, obs, gain_model=settings.gain_model, background=-1, band_weighting="none",
                              **_signal_kwargs(settings))
        wide.line_band_hz = (0.0, 1e6)
        x = q.vector()
        np.testing.assert_allclose(wide.predict(x).residual, full.predict(x).residual, atol=1e-12)
        joint = JointSeries(model, settings, [obs], [1.0], shape="free")
        values = q.values(x)
        lines = np.concatenate([np.asarray(full.transitions(values, c, sy).frequencies_hz)
                                for c, sy in enumerate(q.systems(values))])
        target = float(lines[(lines > 110) & (lines < 290)][0])
        window = (target - 1.5, target + 1.5)
        narrow = MixtureForward(q, obs.restricted([window]), gain_model=settings.gain_model, background=-1,
                                band_weighting="none", **_signal_kwargs(settings))
        narrow.line_band_hz = (window[0] - 2.0, window[1] + 2.0)
        jac = narrow.jacobian(x)
        h = 1e-6
        numeric = np.column_stack([(narrow.predict(x + h * e).residual - narrow.predict(x - h * e).residual) / (2 * h)
                                   for e in np.eye(len(x))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))
        leader = joint.params[0].ties.get(name, name)
        narrow.jacobian_only = {leader}                  # only that column, the same values
        part = narrow.jacobian(x)
        i = q.free_names.index(leader)
        np.testing.assert_allclose(part[:, i], jac[:, i], atol=1e-12)
        table, xs = self._table(joint)
        z = joint.pack(table, xs)
        lower, upper = joint.bounds(20.0)
        before = float(np.sum(narrow.predict(joint.spectrum_vector(z, 0)).residual ** 2))
        cost, zz = joint.local_fit(z, 0, window, [leader], lower=lower, upper=upper)[0]
        self.assertLess(cost, 0.1 * before)
        self.assertAlmostEqual(zz[joint.coupling.index(leader) * joint.m] - z[joint.coupling.index(leader) * joint.m],
                               0.3, delta=0.05)

    @staticmethod
    def _table(joint):
        table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params])
        return table, [p.vector().copy() for p in joint.params]

    def test_spectra_at_one_concentration_share_their_couplings(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        rng = np.random.default_rng(6)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 200.0)], record=acq) for _ in range(3)]
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 0.5])
        self.assertEqual(joint.nn, 2)
        table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params[:2]])
        z = joint.pack(table + np.array([0.0, 0.3])[:, None], [p.vector() for p in joint.params])
        a, b = joint.spectrum_vector(z, 1), joint.spectrum_vector(z, 2)
        for n in joint.coupling:
            self.assertEqual(a[joint.col[n]], b[joint.col[n]])
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_left_out_prediction_recovers_a_bracketed_spectrum(self):
        # three synthetic spectra, J(Ca,Ha) shifted 0 / 0.4 / 0.9 Hz: leave the middle out, fit the outer two with
        # their true couplings, predict the middle within the bracket; the prediction must reach the noise level
        import tempfile
        from fit_joint_series import predict_left_out
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        rng = np.random.default_rng(4)
        tmp = Path(tempfile.mkdtemp())
        name = model.coupling_names["J(Ca,Ha)"][0]
        obs, entries, noise = [], [], 0.002
        for s, shift in enumerate((0.0, 0.4, 0.9)):
            p = settings.parameterize(model.interpretation)
            p.set(name, p.parameters[name].value + shift)
            sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                            band_weighting="none", **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model).real
            y = sig + noise * rng.normal(size=len(f))
            np.save(tmp / f"f{s}.npy", f)
            np.save(tmp / f"y{s}.npy", y)
            entries.append({"id": f"s{s}", "x": [0.1, 0.5, 1.0][s], "freq": str(tmp / f"f{s}.npy"),
                            "values": str(tmp / f"y{s}.npy")})
            obs.append(ObservedSpectrum.from_spectrum(f, y, [(100.0, 200.0)], record=None, real_only=True))
        joint = JointSeries(model, settings, [obs[0], obs[2]], [0.1, 1.0])
        truth = [settings.parameterize(model.interpretation).vector() for _ in range(2)]
        for x, shift in zip(truth, (0.0, 0.9)):
            x[joint.col[name]] += shift
        z = joint.pack(np.array([[x[joint.col[n]] for n in joint.coupling] for x in truth]), truth)
        key_of = {n: k for k, names in model.coupling_names.items() for n in names}
        out = predict_left_out(joint, z, entries[1], (100.0, 200.0), key_of, settings)
        # reference: the middle spectrum fitted at its true couplings (the noise level of this measure)
        p_mid = settings.parameterize(model.interpretation)
        x_mid = p_mid.vector()
        x_mid[joint.col[name]] += 0.4
        fw = MixtureForward(p_mid, obs[1], gain_model=settings.gain_model, background=-1, band_weighting="none",
                            **_signal_kwargs(settings))
        pred = fw.predict(x_mid)
        m = np.ones(len(fw.y), bool)
        noise_level = float(np.linalg.norm(fw.mismatch(pred.model, m)) / np.linalg.norm(fw.mismatch(0 * fw.y, m)))
        self.assertLess(out["relative_residual"], 1.05 * noise_level)
        u = out["u"]["J(Ca,Ha)"]
        self.assertAlmostEqual(u, 0.4 / 0.9, delta=0.1)               # true position inside the bracket


    def test_coordinate_scan_recovers_a_displaced_coupling(self):
        from fit_joint_series import _coordinate_scan
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        obs = []
        for shift in (0.0, 0.3):
            p = settings.parameterize(model.interpretation)
            p.set(model.coupling_names["J(Ca,Hb)"][0], p.parameters[model.coupling_names["J(Ca,Hb)"][0]].value + shift)
            sig = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                                 **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model
            obs.append(ObservedSpectrum.from_spectrum(f, np.asarray(sig), [(100.0, 200.0)], record=acq))
        joint = JointSeries(model, settings, obs, [0.1, 1.0])
        truth = [p.vector() for p in joint.params]
        name = model.coupling_names["J(Ca,Hb)"][0]
        truth[1][joint.col[name]] += 0.3
        table = np.array([[x[joint.col[n]] for n in joint.coupling] for x in truth])
        z_true = joint.pack(table, truth)
        k = joint.coupling.index(name)
        z0 = z_true.copy()
        z0[k * joint.m] += 3.0                                   # far outside the local basin of a 1 Hz line
        lower, upper = joint.bounds(20.0)
        z = _coordinate_scan(joint, z0, lower, upper, 50, cycles=1, half_width=4.0, step=0.25)
        self.assertLess(abs(z[k * joint.m] - z_true[k * joint.m]), 0.3)


if __name__ == "__main__":
    unittest.main()


class SpinSystemFitTests(unittest.TestCase):
    def test_fit_of_a_typed_in_spin_system_recovers_its_variables(self):
        # PLAN 8d: a methyl (13C, 3 x 1H; variable 'a' ties the three 1J) and a 13C-1H pair with a numeric 1J, at
        # weights 1 and 0.5; the fit starts 0.4 Hz off and recovers both couplings and the amplitude ratio
        from scipy.optimize import least_squares
        from zulf_hypothesis.spin_system import spin_system_model
        truth = {"compound": "t", "spin_system": {"components": [
            {"name": "M", "isotopes": ["13C", "1H", "1H", "1H"],
             "J": [[0, "a", "a", "a"], ["a", 0, 0, 0], ["a", 0, 0, 0], ["a", 0, 0, 0]], "weight": 1.0},
            {"name": "X", "isotopes": ["13C", "1H"], "J": [[0, 141.0], [141.0, 0]], "weight": 0.5}],
            "variables": {"a": 126.0}}}
        model = spin_system_model(truth)
        self.assertEqual(sorted(model.coupling_names), ["X:J(1,2)", "a"])
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 300.0, 0.05)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 300.0)], record=acq)
        p = settings.parameterize(model.interpretation)
        sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                        band_weighting="none", **_signal_kwargs(settings)).predict(
            p.vector(), fixed_gains=np.array([1.0, 0.5], complex)).model)
        y = sig + 1e-4 * np.abs(sig).max() * np.random.default_rng(3).normal(size=len(f))
        start = spin_system_model(truth, {"a": 126.4, "X:J(1,2)": 140.6})
        joint = JointSeries(start, _settings_for(start, "ratios", RefineSettings(band_weighting="none",
                                                                                background_order=-1)),
                            [ObservedSpectrum.from_spectrum(f, y, [(100.0, 300.0)], record=acq)], [1.0], shape="free")
        table = np.array([[q.vector()[joint.col[n]] for n in joint.coupling] for q in joint.params])
        z = joint.pack(table, [q.vector().copy() for q in joint.params])
        lo, hi = joint.bounds(5.0)
        sol = least_squares(joint.residual, np.clip(z, lo + 1e-9, hi - 1e-9), jac=joint.jacobian, bounds=(lo, hi),
                            x_scale="jac")
        got = {n: joint.coupling_values(sol.x, k)[0] for k, n in enumerate(joint.coupling)}
        key = {names[0]: k for k, names in start.coupling_names.items()}
        found = {key[n]: v for n, v in got.items()}
        self.assertAlmostEqual(found["a"], 126.0, places=2)
        self.assertAlmostEqual(found["X:J(1,2)"], 141.0, places=2)
        gains = joint.forwards[0].predict(joint.spectrum_vector(sol.x, 0)).gains
        self.assertAlmostEqual(abs(gains[1]) / abs(gains[0]), 0.5, places=2)
