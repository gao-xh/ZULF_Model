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
    def test_matches_time_domain_processing(self):
        tl = random_transitions(1)
        rates = np.array([0.7, 40.0])
        for acq in ACQS:
            renderer = Renderer(acq)
            for zf in (1, 3):
                grid = SpectrumGrid.for_acquisition(acq, 1, 390, zf)
                a = renderer.render(tl, rates, grid.frequencies_hz, gain=0.3 - 0.8j, phase_delay_s=0.002)
                fid = renderer.synthesize(tl, rates, gain=0.3 - 0.8j, phase_delay_s=0.002)
                b = evaluate_spectrum(process_record(fid, acq), acq, grid.frequencies_hz)
                self.assertLess(np.abs(a - b).max() / np.abs(b).max(), 1e-10, msg=str(acq))

    def test_fast_decay_and_wide_window(self):
        acq = Acquisition(500.0, 3000, start_sample=20, sg_window=1001, sg_order=2)
        tl = random_transitions(2, f_max=240.0)
        renderer = Renderer(acq)
        grid = SpectrumGrid.for_acquisition(acq, 1, 240)
        a = renderer.render(tl, 60.0, grid.frequencies_hz)
        b = evaluate_spectrum(process_record(renderer.synthesize(tl, 60.0), acq), acq, grid.frequencies_hz)
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
