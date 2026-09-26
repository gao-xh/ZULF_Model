import unittest
import math

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_model.render import Acquisition, Renderer
from zulf_core.solver.forward import MixtureForward
from zulf_core.solver import (ObservedSpectrum, ParameterPolicy, Parameterization, RefineSettings, refine,
                               refine_candidates)
from zulf_core.spinsystem import Component, Interpretation, SpinSystem, best_permutation

ACQ = Acquisition(1000.0, 6000, start_sample=40, sg_window=101, sg_order=2)
RANGES = [(100.0, 150.0), (230.0, 275.0)]


def methyl_isotopologue(j_ch=126.3, j_far=4.1, j_hh=6.5, j_cm=-3.2):
    iso = ("1H",) * 7 + ("13C",)
    j = np.zeros((8, 8))
    j[:6, 6] = j[6, :6] = j_hh
    j[:3, 7] = j[7, :3] = j_ch
    j[3:6, 7] = j[7, 3:6] = j_far
    j[6, 7] = j[7, 6] = j_cm
    return SpinSystem(iso, j)


def methine_isotopologue(j_ch=135.0, j_hh=6.5, j_cm=-3.0):
    iso = ("1H",) * 7 + ("13C",)
    j = np.zeros((8, 8))
    j[:6, 6] = j[6, :6] = j_hh
    j[:6, 7] = j[7, :6] = j_cm
    j[6, 7] = j[7, 6] = j_ch
    return SpinSystem(iso, j)


def observe(systems, gains, rate=1.5, noise=0.0, seed=0, ranges=RANGES):
    renderer = Renderer(ACQ)
    fid = sum(renderer.synthesize(compute_transitions(s), rate, gain=g) for s, g in zip(systems, gains))
    if noise:
        fid = fid + np.random.default_rng(seed).normal(0, noise, ACQ.points)
    return ObservedSpectrum.from_fid(fid, ACQ, ranges)


def perturbed(system, delta):
    j = system.couplings_hz.copy()
    iso = system.isotopes
    rng = np.random.default_rng(1)
    groups = system.groups
    out = j.copy()
    for a in range(len(groups)):
        for b in range(a + 1, len(groups)):
            d = rng.uniform(-delta, delta)
            out[np.ix_(groups[a], groups[b])] += d
            out[np.ix_(groups[b], groups[a])] += d
    return SpinSystem(iso, out, groups)


class SolverTests(unittest.TestCase):
    def test_single_component_recovery_with_continuation(self):
        truth = methyl_isotopologue()
        obs = observe([truth], [3.0 * np.exp(0.7j)])
        res = refine(Interpretation((Component(perturbed(truth, 1.0)),)), obs, RefineSettings(starts=1))
        self.assertLess(res.relative_residual, 1e-5)
        self.assertLess(best_permutation(truth, res.interpretation.components[0].system).max_abs_error_hz, 1e-3)
        self.assertAlmostEqual(abs(res.gains[0]), 3.0, places=3)
        self.assertAlmostEqual(np.exp(res.parameters["c0.log_rate0"]), 1.5, places=3)
        self.assertEqual(res.flags, [])

    def test_two_component_mixture_shared_phase(self):
        a, b = methyl_isotopologue(), methine_isotopologue()
        phase = np.exp(-1.1j)
        obs = observe([a, b], [2.0 * phase, 1.0 * phase])
        cand = Interpretation((Component(perturbed(a, 0.4)), Component(perturbed(b, 0.4))))
        res = refine(cand, obs, RefineSettings(starts=1))
        self.assertLess(res.relative_residual, 1e-4)
        ratio = abs(res.gains[0]) / abs(res.gains[1])
        self.assertAlmostEqual(ratio, 2.0, places=2)

    def test_frozen_held_out_prediction_ranks_candidates(self):
        truth = methyl_isotopologue()
        train = observe([truth], [1.0], noise=0.02, seed=1)
        held = [observe([truth], [1.0], noise=0.02, seed=s) for s in (2, 3)]
        wrong = methyl_isotopologue(j_ch=129.0, j_far=1.0)
        results = refine_candidates([Interpretation((Component(wrong),)), Interpretation((Component(perturbed(truth, 0.3)),))],
                                    train, held, RefineSettings(starts=1, continuation_rates_per_s=(3.0, 0.0),
                                                                policy=ParameterPolicy(coupling_margin_hz=1.0,
                                                                                       coupling_margin_relative=0.0)))
        self.assertEqual(results[0].candidate_index, 1)
        self.assertLess(results[0].validation_residual, results[1].validation_residual)
        for r in results:
            self.assertEqual(len(r.validation), 2)

    def test_boundary_and_budget_flags(self):
        truth = methyl_isotopologue()
        obs = observe([truth], [1.0])
        start = methyl_isotopologue(j_ch=122.0)
        policy = ParameterPolicy(coupling_margin_hz=2.0, coupling_margin_relative=0.0)
        res = refine(Interpretation((Component(start),)), obs, RefineSettings(starts=1, policy=policy))
        self.assertIn("search_boundary", res.flags)
        res = refine(Interpretation((Component(start),)), obs, RefineSettings(starts=1, max_evaluations=5))
        self.assertIn("search_budget_exhausted", res.flags)

    def test_ties_fixes_and_no_fid(self):
        a, b = methyl_isotopologue(), methine_isotopologue()
        interp = Interpretation((Component(a), Component(b)))
        p = Parameterization.from_interpretation(interp)
        names0, names1 = p.coupling_names(0), p.coupling_names(1)
        p.tie(names0[0], names1[0]).fix(names0[1])
        self.assertNotIn(names1[0], p.free_names)
        self.assertNotIn(names0[1], p.free_names)
        values = p.values(p.vector() + 0.5)
        self.assertEqual(values[names0[0]], values[names1[0]])
        obs = observe([a], [1.0])
        spectrum_only = ObservedSpectrum(obs.frequencies_hz, obs.values, obs.acquisition, obs.band_index)
        res = refine(Interpretation((Component(a),)), spectrum_only, RefineSettings(starts=1))
        self.assertIn("continuation_unavailable_without_fid", res.flags)


class SignVariantTests(unittest.TestCase):
    @staticmethod
    def cn_system(j_nh):
        j = np.zeros((4, 4))
        for (a, b), v in {(0, 2): 140.0, (1, 3): j_nh, (0, 1): -5.0, (0, 3): 3.0, (1, 2): -2.0, (2, 3): 5.0}.items():
            j[a, b] = j[b, a] = v
        return SpinSystem(("13C", "15N", "1H", "1H"), j)

    def test_variants_skip_global_equivalents(self):
        from zulf_core.solver.variants import sign_variants
        pair = Interpretation((Component(SpinSystem(("13C", "1H"), np.array([[0.0, 140.0], [140.0, 0.0]]))),))
        self.assertEqual(sign_variants(pair, 10.0), [])
        variants = sign_variants(Interpretation((Component(self.cn_system(70.0)),)), 10.0)
        self.assertTrue(variants)
        labels = [v.metadata["sign_variant"] for v in variants]
        self.assertEqual(len(labels), len(set(labels)))

    def test_opposite_sign_is_recovered_only_with_variants(self):
        truth = self.cn_system(-70.0)
        obs = observe([truth], [2.0 * np.exp(0.3j)], rate=1.0, ranges=[(40.0, 230.0)])
        wrong = Interpretation((Component(perturbed(self.cn_system(70.0), 0.3)),))
        policy = ParameterPolicy(coupling_margin_hz=5.0, rate_bounds_per_s=(0.1, 10.0))
        plain = refine_candidates([wrong], obs, settings=RefineSettings(starts=1, policy=policy))
        with_variants = refine_candidates([wrong], obs, settings=RefineSettings(starts=1, policy=policy,
                                                                                sign_variants=True))
        def nh(result):
            system = result.interpretation.components[0].system
            match = best_permutation(truth, system)
            aligned = match.sign * system.permute(match.permutation).couplings_hz
            return aligned[1, 3], match.max_abs_error_hz
        self.assertGreater(nh(plain[0])[0], 0)
        value, error = nh(with_variants[0])
        self.assertAlmostEqual(value, -70.0, delta=0.05)
        self.assertLess(error, 0.05)
        self.assertLess(with_variants[0].relative_residual, plain[0].relative_residual)
        self.assertIn("sign_variant", with_variants[0].interpretation.metadata)


class SignalWeightingTests(unittest.TestCase):
    def test_mask_follows_lines_not_broad_background(self):
        from zulf_core.solver.forward import signal_regions
        f = np.arange(100.0, 160.0, 0.05)
        rng = np.random.default_rng(0)
        y = rng.normal(0, 1, len(f)) + 1j * rng.normal(0, 1, len(f))
        y = y + 40.0 / (1 + 1j * (f - 130.0) / 0.3)                  # a narrow line
        y = y + 6.0 * np.exp(-0.5 * ((f - 115.0) / 6.0) ** 2)         # a broad background bump
        from zulf_core.solver.forward import signal_weights
        band = np.zeros(len(f), int)
        cores, sigma = signal_regions(f, y, band)
        self.assertTrue(cores[np.argmin(np.abs(f - 130.0))])
        self.assertFalse(cores[np.argmin(np.abs(f - 115.0))])
        self.assertLess(cores.mean(), 0.1)
        self.assertGreater(sigma, 0.5)
        w = signal_weights(f, band, cores, outside=0.2, taper_hz=2.0)
        self.assertAlmostEqual(w[np.argmin(np.abs(f - 130.0))], 1.0)
        self.assertAlmostEqual(w[np.argmin(np.abs(f - 115.0))], 0.2, places=3)
        near = w[(f > 130.0) & (f < 140.0)]
        self.assertTrue(np.all(np.diff(near[np.argmax(near < 1):]) <= 1e-12))   # smooth, monotone fall-off

    def test_signal_weighting_recovers_couplings_with_background(self):
        truth = methyl_isotopologue()
        acq = Acquisition(1000.0, 6000, start_sample=30)
        t = acq.times()
        signal = Renderer(acq).synthesize(compute_transitions(truth), 1.5, gain=np.exp(0.4j))
        rng = np.random.default_rng(1)
        wobble = 0.3 * np.max(np.abs(signal)) * np.exp(-40.0 * t) * np.cos(2 * np.pi * 105.0 * t)
        obs = ObservedSpectrum.from_fid(signal + wobble + rng.normal(0, 1e-3, acq.points), acq, RANGES)
        start = Interpretation((Component(perturbed(truth, 0.3)),))
        res = refine(start, obs, RefineSettings(starts=1, band_weighting="signal", background_order=1))
        self.assertLess(best_permutation(truth, res.interpretation.components[0].system).max_abs_error_hz, 0.01)
        self.assertIsNotNone(res.signal_region_residual)
        with self.assertRaises(ValueError):
            refine(start, obs, RefineSettings(band_weighting="loud"))

    def test_model_predicted_lines_join_the_cores(self):
        # Data: a single 13C-H line near 128 Hz. Candidate: 13CH3, which also predicts a 2J line near 252 Hz where
        # the data have none; the second pass must add that region to the high-weight cores.
        acq = Acquisition(1000.0, 6000, start_sample=30)
        ch = SpinSystem.from_group_couplings(["13C", "1H"], [1, 1], np.array([[0, 128.0], [128.0, 0]]))
        fid = Renderer(acq).synthesize(compute_transitions(ch), 1.5, gain=np.exp(0.4j))
        obs = ObservedSpectrum.from_fid(fid + np.random.default_rng(1).normal(0, 2e-3, acq.points), acq, RANGES)
        ch3 = SpinSystem.from_group_couplings(["13C", "1H"], [1, 3], np.array([[0, 126.0], [126.0, 0]]))
        base = RefineSettings(starts=1, band_weighting="signal", background_order=1)
        one = refine(Interpretation((Component(ch3),)), obs, base)
        self.assertTrue(any(f.startswith("signal_mask_model_pass") for f in one.flags))
        none = refine(Interpretation((Component(ch3),)), obs, RefineSettings(**{**base.__dict__, "signal_model_passes": 0}))
        self.assertFalse(any(f.startswith("signal_mask_model_pass") for f in none.flags))
        self.assertIsNotNone(one.data_region_residual)

    def test_model_line_heights_come_from_the_transitions(self):
        # Heights read from the transition list match the rendered peak of an isolated line, and a line far below
        # the data's peak-picking threshold (weak, or of either sign) is still located.
        acq = Acquisition(1000.0, 6000, start_sample=30)
        ch = SpinSystem.from_group_couplings(["13C", "1H"], [1, 1], np.array([[0, 128.0], [128.0, 0]]))
        fid = Renderer(acq).synthesize(compute_transitions(ch), 1.5, gain=np.exp(0.4j))
        obs = ObservedSpectrum.from_fid(fid + np.random.default_rng(1).normal(0, 2e-3, acq.points), acq, RANGES)
        param = RefineSettings().parameterize(Interpretation((Component(ch),)))
        forward = MixtureForward(param, obs, band_weighting="signal", background=1)
        values = param.values()
        values = {**values, **{k: math.log(1.5) for k in values if ".log_rate" in k}}
        pred = forward.predict(values=values)
        heights = forward.model_line_heights(values, pred.gains)
        peak = int(np.argmax(heights))
        self.assertLess(abs(forward.f[peak] - 128.0), 0.5 * float(np.median(np.diff(forward.f))) + 1e-9)
        self.assertAlmostEqual(heights[peak] / np.abs(pred.component_spectra[0]).max(), 1.0, delta=0.05)
        weak = forward.model_line_points(values, pred.gains * 1e-3, 1e-3 * heights[peak] / forward.noise_sigma / 2)
        self.assertTrue(weak[peak])
        flipped = forward.model_line_heights(values, -pred.gains)
        np.testing.assert_allclose(flipped, heights)


class ProcessedSpectrumTests(unittest.TestCase):
    """Refinement directly on spectra processed elsewhere (no FID)."""

    @classmethod
    def setUpClass(cls):
        from zulf_core.render.acquisition import evaluate_spectrum, process_record
        from zulf_core.render.grid import SpectrumGrid
        cls.truth = methyl_isotopologue()
        tl = compute_transitions(cls.truth)
        cls.fid = Renderer(ACQ).synthesize(tl, 1.5, gain=3.0 * np.exp(0.7j), phase_delay_s=0.0008)
        cls.f = SpectrumGrid.for_acquisition(ACQ, 1.0, 400.0).frequencies_hz
        cls.complex_values = evaluate_spectrum(process_record(cls.fid, ACQ), ACQ, cls.f)
        cls.start = Interpretation((Component(perturbed(cls.truth, 0.5)),))
        cls.settings = RefineSettings(starts=1)

    def error(self, result):
        return best_permutation(self.truth, result.interpretation.components[0].system).max_abs_error_hz

    def test_complex_spectrum_with_record_matches_fid_route(self):
        obs = ObservedSpectrum.from_spectrum(self.f, self.complex_values, RANGES, record=ACQ)
        self.assertFalse(obs.real_only)
        res = refine(self.start, obs, self.settings)
        self.assertLess(self.error(res), 1e-3)
        self.assertIn("continuation_unavailable_without_fid", res.flags)
        fid_route = refine(self.start, ObservedSpectrum.from_fid(self.fid, ACQ, RANGES), self.settings)
        self.assertLess(abs(res.relative_residual - fid_route.relative_residual), 1e-6)

    def test_phase_corrected_real_spectrum_with_record(self):
        from zulf_core.render.phasing import phase_correct
        absorption = phase_correct(self.complex_values, self.f, 0.7, 0.0008, ACQ).real
        obs = ObservedSpectrum.from_spectrum(self.f, absorption, RANGES, record=ACQ.to_dict(),
                                             phasing={"phase0_rad": 0.7, "delay_s": 0.0008})
        self.assertTrue(obs.real_only)
        res = refine(self.start, obs, self.settings)
        self.assertLess(self.error(res), 1e-3)
        self.assertLess(res.relative_residual, 1e-4)

    def test_phased_real_spectrum_with_record_frame_background(self):
        # A smooth baseline in the record's own frame (residual baseline, ringing leftovers) becomes a ripple of
        # period 1/(crop + delay) after first-order correction; the background columns must follow the correction.
        from zulf_core.render.phasing import phase_correct
        smooth = (0.02 - 0.015j) * np.max(np.abs(self.complex_values)) * (1.0 + 0.3 * (self.f - 200.0) / 200.0)
        absorption = phase_correct(self.complex_values + smooth, self.f, 0.7, 0.0008, ACQ).real
        obs = ObservedSpectrum.from_spectrum(self.f, absorption, RANGES, record=ACQ,
                                             phasing={"phase0_rad": 0.7, "delay_s": 0.0008})
        res = refine(self.start, obs, RefineSettings(starts=1, background_order=1))
        self.assertLess(self.error(res), 1e-3)
        self.assertLess(res.relative_residual, 1e-4)

    def test_real_spectrum_without_record_uses_lorentzian_lines(self):
        from zulf_core.render.acquisition import evaluate_spectrum, process_record
        from zulf_core.render.grid import SpectrumGrid
        from zulf_core.render.phasing import phase_correct
        acq = Acquisition(1000.0, 16384)
        fid = Renderer(acq).synthesize(compute_transitions(self.truth), 1.5, gain=3.0)
        f = SpectrumGrid.for_acquisition(acq, 1.0, 400.0).frequencies_hz
        absorption = phase_correct(evaluate_spectrum(process_record(fid, acq), acq, f), f, 0.0, 0.0, acq).real
        obs = ObservedSpectrum.from_spectrum(f, absorption, RANGES)
        self.assertIsNone(obs.acquisition)
        res = refine(self.start, obs, RefineSettings(starts=1, background_order=0))
        self.assertLess(self.error(res), 0.02)

    def test_table_loader(self):
        import tempfile
        from pathlib import Path
        from zulf_core.io import load_spectrum_table
        with tempfile.TemporaryDirectory() as tmp:
            table = np.column_stack([self.f, self.complex_values.real, self.complex_values.imag])
            np.save(Path(tmp) / "s.npy", table)
            f, v = load_spectrum_table(Path(tmp) / "s.npy")
            np.testing.assert_allclose(v, self.complex_values)
            (Path(tmp) / "s.csv").write_text("freq,absorption\n" + "\n".join(f"{a},{b}" for a, b in
                                                                              zip(self.f, self.complex_values.real)))
            f, v = load_spectrum_table(Path(tmp) / "s.csv")
            self.assertFalse(np.iscomplexobj(v))
            np.testing.assert_allclose(v, self.complex_values.real)
            np.savez(Path(tmp) / "s.npz", frequency_hz=self.f, values=self.complex_values)
            f, v = load_spectrum_table(Path(tmp) / "s.npz")
            np.testing.assert_allclose(f, self.f)


if __name__ == "__main__":
    unittest.main()


class NuisanceTests(unittest.TestCase):
    def test_exponential_baseline_and_ringing_are_absorbed(self):
        truth = methyl_isotopologue()
        acq = Acquisition(1000.0, 6000, start_sample=30)
        renderer = Renderer(acq)
        t = acq.times()
        signal = renderer.synthesize(compute_transitions(truth), 1.5, gain=np.exp(0.4j))
        baseline = -2e3 * np.exp(-0.8 * t) + 400 * np.exp(-12.0 * t)
        ringing = 800 * np.exp(-60 * t) * np.cos(2 * np.pi * 180.0 * t + 0.3)
        obs = ObservedSpectrum.from_fid(signal + baseline + ringing, acq, RANGES)
        nuisance = ({"kind": "exponential", "rate_bounds_per_s": [0.2, 5.0], "initial_rate_per_s": 1.0},
                    {"kind": "exponential", "rate_bounds_per_s": [5.0, 50.0], "initial_rate_per_s": 10.0},
                    {"kind": "damped_sinusoid", "frequency_bounds_hz": [170.0, 190.0], "initial_frequency_hz": 181.0,
                     "rate_bounds_per_s": [20.0, 200.0], "initial_rate_per_s": 50.0})
        start = Interpretation((Component(perturbed(truth, 0.3)),))
        plain = refine(start, obs, RefineSettings(starts=1, continuation_rates_per_s=(3.0, 0.0)))
        modeled = refine(start, obs, RefineSettings(starts=1, continuation_rates_per_s=(3.0, 0.0),
                                                    policy=ParameterPolicy(nuisance=nuisance)))
        err_plain = best_permutation(truth, plain.interpretation.components[0].system).max_abs_error_hz
        err_model = best_permutation(truth, modeled.interpretation.components[0].system).max_abs_error_hz
        self.assertLess(modeled.relative_residual, 1e-4)
        self.assertLess(err_model, 1e-3)
        self.assertGreater(err_plain, 10 * err_model)
        # Nuisance rates are weakly identified from band tails; only a loose check.
        self.assertAlmostEqual(np.exp(modeled.parameters["n0.log_rate"]), 0.8, delta=0.04)
