import dataclasses
import unittest

import numpy as np

from zulf_model.physics import compute_transitions
from zulf_model.render import Acquisition, Renderer
from zulf_model.render.perturb import PerturbationConfig, sample_render_params
from zulf_model.render.phasing import phase1_to_delay_s, phase_correct, reference_delay_s
from zulf_model.render.pipeline import ProcessingConfig, SampleRenderer
from zulf_model.spec import GridSpec, ProblemSpec
from zulf_model.spinsystem import Component, Interpretation, SpinSystem


def ch3():
    j = np.zeros((4, 4))
    j[0, 1:] = j[1:, 0] = 126.0
    return SpinSystem(("13C", "1H", "1H", "1H"), j)


def ch():
    return SpinSystem(("13C", "1H"), np.array([[0.0, 140.0], [140.0, 0.0]]))


def perturbation(phase, delay, **extra):
    base = dict(rate_range_per_s=(1.0, 1.0), gaussian_probability=0.0, family_split_probability=0.0,
                global_phase_range_rad=(phase, phase), phase_delay_range_s=(delay, delay),
                component_gain_log10_range=(0.0, 0.0), component_phase_spread_rad=0.0)
    base.update(extra)
    return PerturbationConfig(**base)


PHASED = ProblemSpec(grid=GridSpec(f_min_hz=50.0, f_max_hz=300.0, channels=("real",), phasing="corrected"))
PROCESSING = ProcessingConfig(mode="fixed", sampling_rate_hz=1000.0, points=8192,
                              fixed={"start_sample": 60, "sg_window": 101, "sg_order": 2},
                              residual_phase0_range_rad=(0.0, 0.0), residual_delay_range_s=(0.0, 0.0))


class PhasingTests(unittest.TestCase):
    def render(self, spec, processing, pert, interp, seed=0):
        return SampleRenderer(spec, processing, pert, noiseless=True).render(interp, np.random.default_rng(seed))

    def test_correction_removes_global_phase_and_delay_exactly(self):
        interp = Interpretation((Component(ch3()),))
        a = self.render(PHASED, PROCESSING, perturbation(2.1, 0.0013), interp)
        b = self.render(PHASED, PROCESSING, perturbation(0.0, 0.0), interp)
        # Exact up to the negative-frequency mirror term of a real FID, whose phase is conjugated.
        np.testing.assert_allclose(a.spectrum, b.spectrum, atol=3e-3 * np.abs(b.spectrum).max())
        self.assertEqual(a.features.shape[0], 1)
        self.assertIn("phasing", a.params)

    def test_corrected_line_is_absorptive_with_physical_sign(self):
        interp = Interpretation((Component(ch()),))
        out = self.render(PHASED, PROCESSING, perturbation(-1.3, 0.0007), interp)
        amplitude = compute_transitions(ch()).amplitudes[0]
        peak = int(np.argmax(np.abs(out.clean)))
        value = out.clean[peak]
        # The peak grid point lies up to half a bin off the line centre, where a Lorentzian is slightly dispersive.
        self.assertLess(abs(value.imag), 0.15 * abs(value.real))
        self.assertEqual(np.sign(value.real), np.sign(amplitude.real))

    def test_residual_error_is_applied_and_recorded(self):
        interp = Interpretation((Component(ch3()),))
        noisy = dataclasses.replace(PROCESSING, residual_phase0_range_rad=(0.3, 0.3),
                                    residual_delay_range_s=(0.0, 0.0))
        exact = self.render(PHASED, PROCESSING, perturbation(0.5, 0.0), interp)
        rotated = self.render(PHASED, noisy, perturbation(0.5, 0.0), interp)
        np.testing.assert_allclose(rotated.spectrum, exact.spectrum * np.exp(0.3j),
                                   atol=1e-9 * np.abs(exact.spectrum).max())
        self.assertAlmostEqual(rotated.params["phasing"]["residual_phase0_rad"], 0.3)

    def test_unphased_default_is_unchanged(self):
        spec = ProblemSpec(grid=GridSpec(f_min_hz=50.0, f_max_hz=300.0))
        out = self.render(spec, PROCESSING, perturbation(2.1, 0.0013), Interpretation((Component(ch3()),)))
        self.assertEqual(out.features.shape[0], 2)
        self.assertNotIn("phasing", out.params)
        with self.assertRaises(ValueError):
            GridSpec(phasing="auto")

    def test_experimental_correction_matches_training_convention(self):
        # A synthetic 'experimental' FID with a known phase and delay, corrected with operator values,
        # equals the training-side corrected render.
        acq = Acquisition(1000.0, 8192, start_sample=60, sg_window=101, sg_order=2)
        tl = compute_transitions(ch3())
        from zulf_model.render.acquisition import evaluate_spectrum, process_record
        f = np.linspace(50.0, 300.0, 2001)
        fid = Renderer(acq).synthesize(tl, 1.0, gain=np.exp(0.8j), phase_delay_s=0.0009)
        fixed = phase_correct(evaluate_spectrum(process_record(fid, acq), acq, f), f, 0.8, 0.0009, acq)
        reference = evaluate_spectrum(process_record(Renderer(acq).synthesize(tl, 1.0), acq), acq, f)
        reference = phase_correct(reference, f, 0.0, 0.0, acq)
        np.testing.assert_allclose(fixed, reference, atol=3e-3 * np.abs(reference).max())
        self.assertAlmostEqual(reference_delay_s(acq), 0.06)
        self.assertAlmostEqual(phase1_to_delay_s(360.0, 1000.0), 0.001)


class PhaseEstimationTests(unittest.TestCase):
    def test_recovers_phase_and_delay_with_signed_lines_and_baseline(self):
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from test_solver import methine_isotopologue, methyl_isotopologue
        from zulf_model.render.acquisition import evaluate_spectrum, process_record
        from zulf_model.render.phasing import estimate_phase
        acq = Acquisition(1000.0, 16384, start_sample=100, sg_window=301, sg_order=2)
        renderer = Renderer(acq)
        t = np.arange(acq.points) / acq.sampling_rate_hz
        for phase0, delay in [(1.1, 0.0008), (-2.5, -0.003)]:
            fid = (renderer.synthesize(compute_transitions(methine_isotopologue()), 1.0, gain=np.exp(1j * phase0),
                                       phase_delay_s=delay)
                   + renderer.synthesize(compute_transitions(methyl_isotopologue()), 1.5, gain=2 * np.exp(1j * phase0),
                                         phase_delay_s=delay)
                   + np.random.default_rng(0).normal(0, 0.01, acq.points) + 50 * np.exp(-t / 0.5))
            f = np.arange(100.0, 300.0, acq.sampling_rate_hz / acq.points / 2)
            estimate = estimate_phase(evaluate_spectrum(process_record(fid, acq), acq, f), f, acq)
            error = (estimate["phase0_rad"] - phase0 + np.pi) % (2 * np.pi) - np.pi
            self.assertLess(abs(error), 0.15)
            self.assertLess(abs(estimate["delay_s"] - delay), 1e-4)
            self.assertAlmostEqual(estimate["reference_delay_s"], 0.1)


class ComponentRatioTests(unittest.TestCase):
    def test_free_mode_sets_rendered_weight_independent_of_abundance(self):
        tls = [compute_transitions(ch3()), compute_transitions(ch())]
        config = perturbation(0.0, 0.0, component_ratio_mode="free", free_weight_log10_range=(-1.0, -1.0))
        params = sample_render_params(np.random.default_rng(0), tls, config, noiseless=True, contributions=[0.011, 0.5])
        weights = [0.011 * abs(params.components[0].gain), 0.5 * abs(params.components[1].gain)]
        np.testing.assert_allclose(weights, [0.1, 0.1])
        with self.assertRaises(ValueError):
            sample_render_params(np.random.default_rng(0), tls, perturbation(0.0, 0.0, component_ratio_mode="x"))

    def test_targets_follow_rendered_weights(self):
        from zulf_model.codec import InterpretationCodec
        from zulf_model.training.data import make_item
        spec = ProblemSpec(grid=GridSpec(f_min_hz=50.0, f_max_hz=300.0), spin_counts=(2, 4))
        interp = Interpretation((Component(ch3(), 1.0), Component(ch(), 1.0)))
        config = perturbation(0.0, 0.0, component_gain_log10_range=(-1.0, 1.0))
        out = SampleRenderer(spec, PROCESSING, config, noiseless=True).render(interp, np.random.default_rng(3))
        gains = [abs(complex(*c["gain"])) for c in out.params["components"]]
        np.testing.assert_allclose([c.contribution for c in out.target.components], gains)
        item = make_item(out, InterpretationCodec(spec))
        top = max(gains)
        expected = sorted((np.log10(g / top) for g in gains), reverse=True)
        np.testing.assert_allclose(sorted(item["log10_contribution"][:2], reverse=True), expected, atol=1e-5)


if __name__ == "__main__":
    unittest.main()
