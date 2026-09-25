import unittest

import numpy as np

from zulf_model.physics import compute_transitions
from zulf_model.render import Acquisition, Renderer
from zulf_model.solver import (ObservedSpectrum, ParameterPolicy, Parameterization, RefineSettings, refine,
                               refine_candidates)
from zulf_model.spinsystem import Component, Interpretation, SpinSystem, best_permutation

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


if __name__ == "__main__":
    unittest.main()
