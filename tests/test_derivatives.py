import sys
import unittest
from pathlib import Path

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_core.physics.derivatives import transition_derivatives
from zulf_core.physics.protocol import Protocol
from zulf_core.physics.transitions import TransitionList
from zulf_core.render.acquisition import Acquisition, evaluate_spectrum, process_record
from zulf_core.render.phasing import phase_correct
from zulf_core.render.renderer import ContinuousRenderer, Renderer
from zulf_core.solver import ObservedSpectrum, ParameterPolicy, Parameterization, RefineSettings, refine
from zulf_core.solver.forward import MixtureForward
from zulf_core.spinsystem import Component, Interpretation, SpinSystem, best_permutation

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_solver import methine_isotopologue, methyl_isotopologue, observe, perturbed  # noqa: E402


def grouped(nuclei, sizes, couplings):
    return SpinSystem.from_group_couplings(nuclei, sizes, np.asarray(couplings, float))


def signal(system, protocol, t):
    tl = compute_transitions(system, protocol)
    return tl.signal(t) - np.real(tl.dc)


def numeric_jacobian(forward, x, h=1e-6):
    cols = []
    for i in range(len(x)):
        step = h * max(1.0, abs(x[i]))
        e = np.zeros(len(x))
        e[i] = step
        cols.append((forward.predict(x + e).residual - forward.predict(x - e).residual) / (2 * step))
    return np.column_stack(cols)


class TransitionDerivativeTests(unittest.TestCase):
    CASES = [
        (["13C", "1H"], [1, 3], [[0, 126], [126, 0]], Protocol()),
        (["13C", "1H", "1H"], [1, 1, 3], [[0, 133.4, -3.5], [133.4, 0, 6.4], [-3.5, 6.4, 0]], Protocol()),
        (["13C", "1H", "1H"], [1, 1, 2], [[0, 140, 0], [140, 0, 7], [0, 7, 0]], Protocol()),     # zero coupling
        (["13C", "1H", "1H"], [1, 1, 1], [[0, 140, 140], [140, 0, 0], [140, 0, 0]], Protocol()),  # degenerate
        (["15N", "1H", "13C", "1H"], [1, 2, 1, 3],
         [[0, -65, -4, 0.8], [-65, 0, 1.5, 6.5], [-4, 1.5, 0, 128], [0.8, 6.5, 128, 0]], Protocol()),
        (["13C", "1H", "1H"], [1, 1, 3], [[0, 133.4, -3.5], [133.4, 0, 6.4], [-3.5, 6.4, 0]],
         Protocol(field_ut=(0.01, 0.0, 0.02))),
    ]

    def test_matches_finite_differences_including_degenerate_levels(self):
        t = np.linspace(0.0, 0.3, 600)
        h = 1e-5
        for nuclei, sizes, couplings, protocol in self.CASES:
            system = grouped(nuclei, sizes, couplings)
            g = len(sizes)
            pairs = [(a, b) for a in range(g) for b in range(a + 1, g)]
            d = transition_derivatives(system, pairs, protocol)
            np.testing.assert_allclose(np.real(np.exp(2j * np.pi * np.outer(t, d.frequencies_hz)) @ d.amplitudes),
                                       signal(system, protocol, t), atol=1e-9)
            got = d.signal_derivative(t)
            gj = np.asarray(couplings, float)
            scale = max(np.abs(signal(system, protocol, t)).max(), 1.0)
            for p, (a, b) in enumerate(pairs):
                out = []
                for sign in (1.0, -1.0):
                    m = gj.copy()
                    m[a, b] += sign * h
                    m[b, a] += sign * h
                    out.append(signal(grouped(nuclei, sizes, m), protocol, t))
                fd = (out[0] - out[1]) / (2 * h)
                self.assertLess(np.abs(fd - got[p]).max() / scale, 1e-6, msg=f"{nuclei} {sizes} pair {(a, b)}")


class RendererDirectionTests(unittest.TestCase):
    def test_frequency_rate_and_delay_directions(self):
        tl = compute_transitions(methine_isotopologue())
        f, amp = tl.frequencies_hz, tl.amplitudes
        rates = np.linspace(0.8, 1.5, len(f))
        grid = np.linspace(100.0, 300.0, 1200)
        tau = 7e-4
        rng = np.random.default_rng(0)
        df, dr = rng.normal(size=len(f)), rng.normal(size=len(f))
        b = np.stack([2j * np.pi * tau * amp * df, np.zeros_like(amp), 2j * np.pi * f * amp])
        c = np.stack([2j * np.pi * amp * df, -amp * dr, np.zeros_like(amp)])
        acq = Acquisition(1000.0, 4096, start_sample=40, sg_window=101, sg_order=2, remove_mean=True,
                          time_origin_s=0.003, apodization_rate_per_s=2.0)
        for renderer in (Renderer(acq), ContinuousRenderer(2.0)):
            got = renderer.render_pair_directions(f, rates, np.zeros(len(f)), b, c, grid, tau)

            def pair(ff=f, rr=rates, tt=tau):
                return renderer.render_pair(TransitionList(ff, amp, families=np.arange(len(f))), rr, grid, tt)

            h = 1e-6
            numeric = [(pair(ff=f + h * df) - pair(ff=f - h * df)) / (2 * h),
                       (pair(rr=rates + h * dr) - pair(rr=rates - h * dr)) / (2 * h),
                       (pair(tt=tau + h) - pair(tt=tau - h)) / (2 * h)]
            for p, fd in enumerate(numeric):
                self.assertLess(np.abs(fd - got[:, :, p]).max() / np.abs(fd).max(), 1e-5)
        with self.assertRaises(NotImplementedError):
            ContinuousRenderer().render_pair_directions(f, rates, np.full(len(f), 0.1), b, c, grid)


class ForwardJacobianTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.truths = [methyl_isotopologue(), methine_isotopologue()]
        cls.obs = observe(cls.truths, [2.0 * np.exp(0.4j), 1.2 * np.exp(0.4j)])
        cls.start = Interpretation(tuple(Component(perturbed(t, 0.5)) for t in cls.truths))

    def check(self, forward, x, tolerance=1e-3):
        forward.predict(x)
        analytic = forward.jacobian(x)
        numeric = numeric_jacobian(forward, x)
        error = np.linalg.norm(analytic - numeric, axis=0) / np.maximum(np.linalg.norm(numeric, axis=0), 1e-9)
        self.assertLess(error.max(), tolerance, msg=str(dict(zip(forward.p.free_names, error))))

    def test_exact_away_from_the_optimum(self):
        nuisance = ({"kind": "exponential", "rate_bounds_per_s": [0.2, 5.0], "initial_rate_per_s": 1.0},)
        for gain_model, background, policy in [
                ("shared_phase", -1, ParameterPolicy()),
                ("complex", 0, ParameterPolicy()),
                # Polynomial background and exponential nuisance terms are tested separately: together they are
                # nearly collinear over a band, and the objective itself is then too noisy for a numeric reference.
                ("shared_phase", 1, ParameterPolicy(family_edges_hz=(200.0,))),
                ("shared_phase", -1, ParameterPolicy(family_edges_hz=(200.0,), nuisance=nuisance)),
                ("shared_phase", -1, ParameterPolicy(fit_sigma=True, initial_sigma_hz=0.2))]:
            param = Parameterization.from_interpretation(self.start, policy)
            param.tie("c0.J1-2", "c1.J1-2")
            forward = MixtureForward(param, self.obs, gain_model=gain_model, background=background)
            self.check(forward, param.vector() + 0.003)

    def test_exact_with_fixed_amplitude_ratios(self):
        for gain_model, background in [("shared_phase", 1), ("complex", 0)]:
            param = Parameterization.from_interpretation(self.start)
            param.tie("c0.J1-2", "c1.J1-2")
            param.tie("c0.log_rate0", "c1.log_rate0")
            forward = MixtureForward(param, self.obs, gain_model=gain_model, background=background,
                                     amplitude_ratios=(1.0, 0.6))
            self.check(forward, param.vector() + 0.003)

    def test_fixed_amplitude_ratios_in_refinement(self):
        start = Interpretation(tuple(Component(perturbed(t, 0.3)) for t in self.truths))
        res = refine(start, self.obs, RefineSettings(starts=1, amplitude_ratios=(1.0, 0.6)))
        gains = np.abs(res.gains)
        self.assertAlmostEqual(gains[1] / gains[0], 0.6, places=12)
        for truth, comp in zip(self.truths, res.interpretation.components):
            self.assertLess(best_permutation(truth, comp.system).max_abs_error_hz, 0.01)
        with self.assertRaises(ValueError):
            refine(start, self.obs, RefineSettings(starts=1, amplitude_ratios=(1.0,)))

    def test_processed_real_spectrum_and_kaufman_at_zero_residual(self):
        acq = Acquisition(1000.0, 4096, start_sample=20)
        fid = Renderer(acq).synthesize(compute_transitions(self.truths[0]), 1.5, gain=np.exp(0.7j),
                                       phase_delay_s=4e-4)
        f = np.arange(100.0, 300.0, 0.1)
        phased = phase_correct(evaluate_spectrum(process_record(fid, acq), acq, f), f, 0.7, 4e-4, acq)
        for record in (acq, None):
            obs = ObservedSpectrum.from_spectrum(f, phased.real, [(110.0, 150.0), (240.0, 270.0)], record=record,
                                                 phasing={"phase0_rad": 0.7, "delay_s": 4e-4})
            param = Parameterization.from_interpretation(Interpretation((Component(perturbed(self.truths[0], 0.3)),)))
            forward = MixtureForward(param, obs, background=0)
            self.check(forward, param.vector() + 0.002)
        truth = Parameterization.from_interpretation(Interpretation((Component(self.truths[0]),)))
        x = truth.vector()
        x[truth.free_names.index("c0.log_rate0")] = np.log(1.5)
        obs = ObservedSpectrum.from_fid(fid, acq, [(110.0, 150.0), (240.0, 270.0)])
        forward = MixtureForward(truth, obs)
        x[truth.free_names.index("phase_delay")] = 4e-4
        forward.predict(x)
        np.testing.assert_allclose(forward.jacobian(x, full=False), forward.jacobian(x), atol=1e-8)

    def test_refine_modes_agree(self):
        truth = self.truths[0]
        obs = observe([truth], [3.0 * np.exp(0.7j)], noise=1e-3)
        start = Interpretation((Component(perturbed(truth, 0.3)),))
        results = {}
        for mode in ("analytic", "kaufman", "finite_difference"):
            results[mode] = refine(start, obs, RefineSettings(starts=1, jacobian=mode))
            err = best_permutation(truth, results[mode].interpretation.components[0].system).max_abs_error_hz
            self.assertLess(err, 0.01, msg=mode)
        self.assertGreater(results["analytic"].jacobian_evaluations, 0)
        self.assertEqual(results["finite_difference"].jacobian_evaluations, 0)
        self.assertLess(results["analytic"].evaluations, results["finite_difference"].evaluations)
        with self.assertRaises(ValueError):
            refine(start, obs, RefineSettings(jacobian="numeric"))


if __name__ == "__main__":
    unittest.main()
