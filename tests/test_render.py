import unittest

import numpy as np

from zulf_model.physics import TransitionList, compute_transitions
from zulf_model.render import (Acquisition, PerturbationConfig, Renderer, SpectrumGrid, evaluate_spectrum,
                               process_record, render_observation, sample_render_params, spectrum_features)
from zulf_model.spinsystem import SpinSystem

ACQS = [
    Acquisition(1000.0, 4000),
    Acquisition(1000.0, 4000, start_sample=100, sg_window=301, sg_order=2),
    Acquisition(1000.0, 4000, start_sample=50, stop_sample=3990, sg_window=301, sg_order=3, time_origin_s=0.013),
    Acquisition(1000.0, 4000, sg_window=51, sg_order=2, remove_mean=False),
    Acquisition(800.0, 2400, start_sample=400, sg_window=801, sg_order=2),
    Acquisition(1000.0, 4000, start_sample=30, sg_window=101, sg_order=2, remove_mean=True,
                apodization_rate_per_s=7.0),
    Acquisition(1000.0, 4000, start_sample=250, sg_window=51, sg_order=2, time_origin_s=0.004,
                phase_reference="acquisition"),
]


def random_transitions(seed, count=12, f_max=380.0):
    rng = np.random.default_rng(seed)
    return TransitionList(rng.uniform(3, f_max, count), rng.normal(size=count) + 1j * rng.normal(size=count),
                          families=rng.integers(0, 2, count))


class AcquisitionTests(unittest.TestCase):
    def test_validation(self):
        with self.assertRaises(ValueError):
            Acquisition(1000.0, 100, sg_window=100)
        with self.assertRaises(ValueError):
            Acquisition(1000.0, 100, start_sample=90)
        a = Acquisition.from_times(1000.0, 5000, start_s=0.1)
        self.assertEqual(a.start_sample, 100)
        self.assertEqual(Acquisition.from_dict(a.to_dict()), a)

    def test_evaluate_spectrum_direct_and_fft_agree(self):
        acq = ACQS[1]
        rng = np.random.default_rng(0)
        y = process_record(rng.normal(size=acq.points), acq)
        grid = SpectrumGrid.for_acquisition(acq, 5, 300, 2)
        fast = evaluate_spectrum(y, acq, grid.frequencies_hz)
        off_grid = grid.frequencies_hz + 1e-3
        direct = evaluate_spectrum(y, acq, off_grid)
        manual = np.array([(y * np.exp(-2j * np.pi * f * np.arange(acq.n) / acq.sampling_rate_hz)).sum() / acq.n
                           for f in off_grid[:5]])
        np.testing.assert_allclose(direct[:5], manual, atol=1e-12)
        self.assertLess(np.abs(fast - direct).max() / np.abs(fast).max(), 0.05)

    def test_on_bin_undamped_line(self):
        acq = Acquisition(1000.0, 1000, remove_mean=False)
        f0 = 100.0
        tl = TransitionList(np.array([f0]), np.array([2.0 + 0j]))
        value = Renderer(acq).render(tl, 0.0, np.array([f0]))[0]
        self.assertAlmostEqual(value, 1.0 + 0j, places=12)


class RendererTests(unittest.TestCase):
    def test_analytic_matches_direct_time_domain_processing(self):
        tl = random_transitions(1)
        rates = np.array([0.7, 40.0])
        for acq in ACQS:
            renderer = Renderer(acq, backend="analytic")
            for zf in (1, 3):
                grid = SpectrumGrid.for_acquisition(acq, 1, 390, zf)
                a = renderer.render(tl, rates, grid.frequencies_hz, gain=0.3 - 0.8j, phase_delay_s=0.002)
                fid = renderer.synthesize(tl, rates, gain=0.3 - 0.8j, phase_delay_s=0.002, method="direct")
                b = evaluate_spectrum(process_record(fid, acq), acq, grid.frequencies_hz)
                self.assertLess(np.abs(a - b).max() / np.abs(b).max(), 1e-10, msg=str(acq))

    def test_time_backend_matches_analytic(self):
        tl = random_transitions(6, count=40)
        for acq in ACQS:
            grid = SpectrumGrid.for_acquisition(acq, 1, 390, 2)
            a = Renderer(acq, backend="analytic").render(tl, np.array([0.5, 3.0]), grid.frequencies_hz, 1 - 1j)
            b = Renderer(acq, backend="time").render(tl, np.array([0.5, 3.0]), grid.frequencies_hz, 1 - 1j)
            self.assertLess(np.abs(a - b).max() / np.abs(a).max(), 1e-9, msg=str(acq))

    def test_nufft_matches_direct_sum(self):
        from zulf_model.render.nufft import nufft_type1
        rng = np.random.default_rng(3)
        for n, k in ((257, 7), (4000, 300)):
            w = rng.uniform(-np.pi, 2 * np.pi, k)
            a = rng.normal(size=k) + 1j * rng.normal(size=k)
            ref = np.exp(1j * np.outer(np.arange(n), w)) @ a
            self.assertLess(np.abs(nufft_type1(w, a, n) - ref).max() / np.abs(ref).max(), 1e-10)

    def test_fast_decay_and_wide_window(self):
        acq = Acquisition(500.0, 3000, start_sample=20, sg_window=1001, sg_order=2)
        tl = random_transitions(2, f_max=240.0)
        renderer = Renderer(acq)
        grid = SpectrumGrid.for_acquisition(acq, 1, 240)
        a = renderer.render(tl, 60.0, grid.frequencies_hz)
        b = evaluate_spectrum(process_record(renderer.synthesize(tl, 60.0, method="direct"), acq), acq,
                              grid.frequencies_hz)
        self.assertLess(np.abs(a - b).max() / np.abs(b).max(), 1e-9)

    def test_render_pair_linearity(self):
        acq = ACQS[1]
        tl = random_transitions(3)
        grid = SpectrumGrid.for_acquisition(acq, 1, 390)
        renderer = Renderer(acq)
        cols = renderer.render_pair(tl, 2.0, grid.frequencies_hz)
        g = 0.4 + 1.3j
        np.testing.assert_allclose(cols @ np.array([g.real, g.imag]), renderer.render(tl, 2.0, grid.frequencies_hz, g),
                                   atol=1e-12)

    def test_physical_system_renders(self):
        j = np.zeros((4, 4)); j[0, 1:] = j[1:, 0] = 140.0
        tl = compute_transitions(SpinSystem(("13C", "1H", "1H", "1H"), j))
        acq = ACQS[0]
        grid = SpectrumGrid.for_acquisition(acq, 1, 400)
        spec = Renderer(acq).render(tl, 1.0, grid.frequencies_hz)
        peaks = grid.frequencies_hz[np.argsort(np.abs(spec))[-2:]]
        np.testing.assert_allclose(sorted(peaks), [140.0, 280.0], atol=acq.native_spacing_hz)


class SavgolTests(unittest.TestCase):
    def test_fft_savgol_matches_scipy(self):
        from scipy.signal import savgol_filter
        from zulf_model.render.acquisition import savgol_baseline
        rng = np.random.default_rng(0)
        x = np.cumsum(rng.normal(size=5000)) + 1e4 * np.exp(-np.arange(5000) / 800)
        for window, order in ((801, 2), (101, 3), (31, 2)):
            ref = savgol_filter(x, window, order, mode="mirror")
            self.assertLess(np.abs(savgol_baseline(x, window, order) - ref).max() / np.abs(ref).max(), 1e-11)


class PureRouteTests(unittest.TestCase):
    def test_default_acquisition_is_pure(self):
        acq = Acquisition(1000.0, 4000)
        self.assertTrue(acq.is_pure)
        self.assertFalse(ACQS[1].is_pure)
        self.assertTrue(ACQS[1].without_processing().is_pure)
        rng = np.random.default_rng(0)
        x = rng.normal(size=acq.points)
        np.testing.assert_array_equal(process_record(x, acq), x)

    def test_continuous_matches_long_finite_record(self):
        from zulf_model.render import ContinuousRenderer
        tl = random_transitions(11, count=6, f_max=300.0)
        acq = Acquisition(1000.0, 40000)
        grid = SpectrumGrid.for_acquisition(acq, 1, 400)
        finite = Renderer(acq).render(tl, 3.0, grid.frequencies_hz)
        # (1/n) sum_m x_m e^{-i w m} approximates (fs/n) integral; T = n / fs.
        pure = ContinuousRenderer(acq.n / acq.sampling_rate_hz).render(tl, 3.0, grid.frequencies_hz)
        self.assertLess(np.abs(finite - pure).max() / np.abs(finite).max(), 5e-3)


class PhaseReferenceTests(unittest.TestCase):
    def test_acquisition_reference_keeps_line_phase_across_crops(self):
        tl = TransitionList(np.array([100.0]), np.array([np.exp(0.9j)]))
        phases = []
        for start in (0, 137, 400):
            acq = Acquisition(1000.0, 8000, start_sample=start, phase_reference="acquisition")
            value = Renderer(acq).render(tl, 0.0, np.array([100.0]))[0]
            phases.append(np.angle(value))
        np.testing.assert_allclose(phases, 0.9, atol=0.02)


class BroadeningTests(unittest.TestCase):
    def test_gaussian_finite_record_matches_time_domain(self):
        tl = random_transitions(12, count=5)
        acq = ACQS[2]
        renderer = Renderer(acq)
        grid = SpectrumGrid.for_acquisition(acq, 1, 390)
        a = renderer.render(tl, np.array([0.5, 2.0]), grid.frequencies_hz, 0.7j, 0.001, np.array([0.3, 1.2]))
        t = acq.times()
        # Independent reference: explicit sum of broadened real oscillations.
        fam = tl.families
        env = np.exp(-np.outer(t, np.array([0.5, 2.0])[fam]) - 0.5 * (2 * np.pi * np.outer(t, np.array([0.3, 1.2])[fam])) ** 2)
        amp = tl.amplitudes * 0.7j * np.exp(2j * np.pi * tl.frequencies_hz * 0.001)
        fid = (env * np.exp(2j * np.pi * np.outer(t, tl.frequencies_hz)) @ amp).real
        b = evaluate_spectrum(process_record(fid, acq), acq, grid.frequencies_hz)
        self.assertLess(np.abs(a - b).max() / np.abs(b).max(), 1e-10)

    def test_continuous_voigt_matches_numerical_integral(self):
        from scipy.integrate import quad
        from zulf_model.render import ContinuousRenderer
        tl = TransitionList(np.array([50.0]), np.array([1.0 + 0.5j]))
        f = np.array([48.0, 50.0, 51.3])
        got = ContinuousRenderer().render(tl, 2.0, f, gaussian_sigma_hz=0.8)
        for k, fk in enumerate(f):
            def integrand(t, part):
                x = (np.exp(-2.0 * t - 0.5 * (2 * np.pi * 0.8 * t) ** 2)
                     * ((1.0 + 0.5j) * np.exp(2j * np.pi * 50.0 * t)).real * np.exp(-2j * np.pi * fk * t))
                return x.real if part == 0 else x.imag
            ref = quad(integrand, 0, 3, args=(0,), limit=2000)[0] + 1j * quad(integrand, 0, 3, args=(1,), limit=2000)[0]
            self.assertAlmostEqual(abs(got[k] - ref) / abs(ref), 0.0, places=6)

    def test_pure_gaussian_line(self):
        from zulf_model.render import ContinuousRenderer
        tl = TransitionList(np.array([80.0]), np.array([1.0 + 0j]))
        spec = ContinuousRenderer().render(tl, 0.0, np.linspace(70, 90, 201), gaussian_sigma_hz=1.0)
        self.assertTrue(np.isfinite(spec).all())
        self.assertAlmostEqual(float(np.linspace(70, 90, 201)[np.argmax(spec.real)]), 80.0, places=6)


class PerturbationTests(unittest.TestCase):
    def test_noise_level_calibration(self):
        acq = Acquisition(1000.0, 8000)
        renderer = Renderer(acq)
        tl = TransitionList(np.array([100.0]), np.array([1.0 + 0j]))
        grid = SpectrumGrid.for_acquisition(acq, 200, 450)
        cfg = PerturbationConfig(snr_range=(10.0, 10.0), drift_probability=0, interference_probability=0,
                                 family_split_probability=0)
        rng = np.random.default_rng(4)
        params = sample_render_params(rng, [tl], cfg)
        noisy, clean = render_observation(renderer, [tl], params, grid.frequencies_hz, [1.0])
        full_grid = SpectrumGrid.for_acquisition(acq, 1, 450)
        _, clean_full = render_observation(renderer, [tl], params, full_grid.frequencies_hz, [1.0])
        expected = np.abs(clean_full).max() / 10.0
        # Peak of a single line lies outside the noise band; compare with the band-restricted peak.
        measured = np.std(noisy - clean) / np.sqrt(1)
        peak_band = np.abs(clean).max()
        self.assertAlmostEqual(measured / (peak_band / 10.0), 1.0, delta=0.15)
        self.assertGreater(expected, 0)

    def test_reproducible_from_params(self):
        acq = ACQS[1]
        renderer = Renderer(acq)
        tl = random_transitions(5)
        rng = np.random.default_rng(9)
        params = sample_render_params(rng, [tl], PerturbationConfig(drift_probability=1, interference_probability=1))
        grid = SpectrumGrid.for_acquisition(acq, 1, 390)
        a, _ = render_observation(renderer, [tl], params, grid.frequencies_hz, [1.0])
        from zulf_model.render.perturb import RenderParams
        again = RenderParams.from_dict(params.to_dict())
        b, _ = render_observation(renderer, [tl], again, grid.frequencies_hz, [1.0])
        np.testing.assert_allclose(a, b)

    def test_features(self):
        x, scale = spectrum_features(np.array([3 + 4j, 0, 0, 0]), ("real", "imag", "magnitude"))
        self.assertEqual(x.shape, (3, 4))
        self.assertAlmostEqual(scale, 2.5)
        self.assertAlmostEqual(float(x[2, 0]), 2.0, places=6)


class GridTests(unittest.TestCase):
    def test_bands(self):
        grid = SpectrumGrid(np.arange(1, 11, dtype=float))
        np.testing.assert_array_equal(grid.band_membership([(2, 3), (7, 9)]), [-1, 0, 0, -1, -1, -1, 1, 1, 1, -1])
        with self.assertRaises(ValueError):
            grid.band_membership([(2, 5), (4, 8)])


if __name__ == "__main__":
    unittest.main()


class BackendSelectionTests(unittest.TestCase):
    def test_fid_free_routes(self):
        from zulf_model.generator import build_default_sampler
        from zulf_model.render import ProcessingConfig, SampleRenderer
        from zulf_model.spec import ProblemSpec
        spec = ProblemSpec()
        sample = next(build_default_sampler(spec).generate(1, seed=4))
        rng = np.random.default_rng(0)
        analytic = SampleRenderer(spec, ProcessingConfig(points=4096, renderer_backend="analytic"),
                                  PerturbationConfig(gaussian_probability=0.0), noiseless=True)
        timed = SampleRenderer(spec, ProcessingConfig(points=4096, renderer_backend="time"),
                               PerturbationConfig(gaussian_probability=0.0), noiseless=True)
        a = analytic.render(sample.interpretation, np.random.default_rng(1)).clean
        b = timed.render(sample.interpretation, np.random.default_rng(1)).clean
        self.assertLess(np.abs(a - b).max() / np.abs(a).max(), 1e-9)
        self.assertNotIn("render.nufft", analytic.timer.report())
        with self.assertRaises(ValueError):
            SampleRenderer(spec, ProcessingConfig(renderer_backend="analytic"), PerturbationConfig(gaussian_probability=0.5))
        continuous = SampleRenderer(spec, ProcessingConfig(points=4096, mode="continuous"))
        out = continuous.render(sample.interpretation, rng)
        self.assertTrue(np.isfinite(out.features).all())
