import unittest

import numpy as np

from zulf_core.diagnostics import diagnose_fid, fit_exponentials, plateau_end
from zulf_core.io import decode_dat


def synthetic_fid(fs=4000.0, points=40000, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(points) / fs
    molecular = 3.0 * np.exp(-1.5 * t) * np.cos(2 * np.pi * 126.0 * t) + 2.0 * np.exp(-2.0 * t) * np.cos(2 * np.pi * 252.0 * t)
    baseline = -1.5e4 * np.exp(-0.9 * t) + 2e3 * np.exp(-8.0 * t)
    ringing = 2e4 * np.exp(-250.0 * t) * np.cos(2 * np.pi * 480.0 * t)
    fid = molecular + baseline + ringing + rng.normal(0, 0.4, points)
    fid[1:12] = 28848.0 + rng.normal(0, 1.0, 11)      # saturation plateau
    fid[0] = -46.0                                     # anomalous first sample
    return fid, fs


class DiagnosticsTests(unittest.TestCase):
    def test_detects_first_point_plateau_and_ringing(self):
        fid, fs = synthetic_fid()
        d = diagnose_fid(fid, fs)
        self.assertTrue(d.first_point_anomaly)
        self.assertAlmostEqual(d.plateau_end_s, 12 / fs, delta=2 / fs)
        # Ringing envelope exp(-250 t) from 2e4 reaches the few-unit level near 35 ms.
        self.assertGreater(d.ringing_end_s, 0.015)
        self.assertLess(d.ringing_end_s, 0.08)
        self.assertTrue(d.candidate_recipes)
        for recipe in d.candidate_recipes:
            self.assertGreater(recipe["acquisition"]["start_sample"], 12)

    def test_no_false_plateau_on_clean_record(self):
        t = np.arange(8000) / 1000.0
        fid = np.exp(-t) * np.cos(2 * np.pi * 50 * t) + np.random.default_rng(1).normal(0, 1e-3, len(t))
        self.assertEqual(plateau_end(fid)[0], 0)
        self.assertFalse(diagnose_fid(fid, 1000.0).first_point_anomaly)

    def test_exponential_fit_recovers_rates(self):
        t = np.linspace(0, 10, 5000)
        y = -1e4 * np.exp(-0.9 * t) + 2e3 * np.exp(-8.0 * t) + 5.0
        fit = fit_exponentials(t, y, 2)
        self.assertTrue(fit["converged"])
        np.testing.assert_allclose(sorted(fit["rates_per_s"]), [0.9, 8.0], rtol=1e-3)

    def test_dat_decoder_convention(self):
        # Legacy convention: reverse the bytes, read little-endian int16, keep [20:-2], reverse samples.
        words = np.arange(-40, 60, dtype="<i2")
        payload = words.tobytes()[::-1]
        np.testing.assert_array_equal(decode_dat(payload), words[20:-2][::-1].astype(float))


if __name__ == "__main__":
    unittest.main()
