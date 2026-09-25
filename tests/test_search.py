import sys
import unittest
from pathlib import Path

import numpy as np

from zulf_model.physics import compute_transitions
from zulf_model.physics.protocol import SUDDEN_DROP, Protocol
from zulf_model.physics.transitions import reference_signal
from zulf_model.solver import (ParameterPolicy, Parameterization, PatternObjective, RefineSettings, SearchSettings,
                               global_search, refine)
from zulf_model.spinsystem import Component, Interpretation, SpinSystem, best_permutation

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_solver import methine_isotopologue, observe  # noqa: E402

POLICY = ParameterPolicy(coupling_margin_hz=8.0, rate_bounds_per_s=(0.1, 10.0))


class ResidualFieldTests(unittest.TestCase):
    def test_sectors_match_brute_force_in_any_field(self):
        j = np.zeros((5, 5))
        j[0, 1:4] = j[1:4, 0] = 125.0
        j[0, 4] = j[4, 0] = -4.0
        j[1:4, 4] = j[4, 1:4] = 6.5
        system = SpinSystem(("13C", "1H", "1H", "1H", "1H"), j)
        t = np.linspace(0, 0.05, 40)
        for field in [(0, 0, 0.05), (0.05, 0, 0), (0.02, 0.03, 0.01)]:
            protocol = Protocol(field_ut=field)
            ref = reference_signal(system, t, protocol)
            for method in ("sectors", "full"):
                sig = compute_transitions(system, protocol, method=method).signal(t)
                self.assertLess(np.max(np.abs(sig - ref)) / np.max(np.abs(ref)), 1e-10, msg=f"{field} {method}")

    def test_zero_field_default_and_serialization(self):
        self.assertFalse(SUDDEN_DROP.has_field)
        protocol = Protocol(field_ut=(0.0, 0.0, 0.1))
        again = Protocol.from_dict(protocol.to_dict())
        self.assertEqual(again, protocol)
        self.assertEqual(SUDDEN_DROP.with_field((0, 0, 0.1)), protocol)
        with self.assertRaises(ValueError):
            Protocol(field_ut=(0.0, 0.1))

    def test_field_shifts_and_splits_as_expected(self):
        # 13C-1H pair: a longitudinal field only shifts the F=1,m=0 <-> F=0 line at second order;
        # a transverse field splits it by the F=1 Zeeman frequency (gamma_H + gamma_C) / 2 * B.
        from zulf_model.nuclei import get_registry
        system = SpinSystem(("13C", "1H"), np.array([[0.0, 140.0], [140.0, 0.0]]))
        longitudinal = compute_transitions(system, Protocol(field_ut=(0.0, 0.0, 0.01)))
        self.assertEqual(len(longitudinal), 1)
        self.assertLess(abs(longitudinal.frequencies_hz[0] - 140.0), 1e-3)
        transverse = compute_transitions(system, Protocol(field_ut=(0.01, 0.0, 0.0)))
        registry = get_registry()
        zeeman = 0.5 * (registry.gamma("1H") + registry.gamma("13C")) * 0.01
        lines = np.sort(np.asarray(transverse.frequencies_hz))
        high = lines[lines > 100]
        self.assertEqual(len(high), 2)
        self.assertAlmostEqual(high[1] - high[0], 2 * zeeman, delta=0.01 * zeeman)


class GlobalSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.truth = methine_isotopologue(133.4, 6.4, -3.5)
        cls.obs = observe([cls.truth], [2.0 * np.exp(0.4j)], rate=0.8, noise=1e-4, ranges=[(120.0, 150.0)])
        cls.far = Interpretation((Component(methine_isotopologue(137.5, 5.0, -1.5)),))

    def test_pattern_objective_prefers_truth(self):
        true_param = Parameterization.from_interpretation(Interpretation((Component(self.truth),)), POLICY)
        far_param = Parameterization.from_interpretation(self.far, POLICY)
        objective = PatternObjective(true_param, self.obs)
        self.assertLess(objective(true_param.values()), 0.1)
        self.assertGreater(objective(far_param.values()), 0.5)

    def test_search_supplies_starts_that_local_refinement_misses(self):
        local_only = refine(self.far, self.obs, RefineSettings(starts=1, policy=POLICY, continuation_rates_per_s=(0.0,)))
        self.assertGreater(best_permutation(self.truth, local_only.interpretation.components[0].system).max_abs_error_hz, 0.5)
        settings = RefineSettings(policy=POLICY, continuation_rates_per_s=(0.0,),
                                  search=dict(popsize=8, maxiter=25, max_seconds=60, solutions=2))
        result = refine(self.far, self.obs, settings)
        self.assertLess(best_permutation(self.truth, result.interpretation.components[0].system).max_abs_error_hz, 1e-3)
        self.assertLess(result.relative_residual, 1e-4)
        self.assertIsNotNone(result.search)
        self.assertEqual(len(result.search["starts"]), len(result.search["costs"]))

    def test_distinct_starts_and_explicit_points(self):
        param = Parameterization.from_interpretation(self.far, POLICY)
        found = global_search(param, self.obs, SearchSettings(popsize=8, maxiter=10, solutions=3, polish=False))
        self.assertGreaterEqual(len(found.points), 1)
        self.assertEqual(found.costs, sorted(found.costs))
        couplings = [i for i, n in enumerate(param.free_names) if param.parameters[n].kind == "coupling"]
        for a in range(len(found.points)):
            for b in range(a + 1, len(found.points)):
                gap = np.max(np.abs(found.points[a][couplings] - found.points[b][couplings]))
                self.assertGreaterEqual(gap, 0.5)
        with self.assertRaises(ValueError):
            refine(self.far, self.obs, RefineSettings(policy=POLICY), initial_points=[np.zeros(2)])

    def test_dual_annealing_option(self):
        settings = RefineSettings(policy=POLICY, continuation_rates_per_s=(0.0,),
                                  search=dict(method="dual_annealing", annealing_maxiter=150, max_seconds=60,
                                              solutions=2))
        result = refine(self.far, self.obs, settings)
        self.assertLess(best_permutation(self.truth, result.interpretation.components[0].system).max_abs_error_hz, 1e-3)
        with self.assertRaises(ValueError):
            global_search(Parameterization.from_interpretation(self.far, POLICY), self.obs,
                          SearchSettings(method="annealing"))

    def test_search_requires_fid(self):
        from zulf_model.solver import ObservedSpectrum
        bare = ObservedSpectrum(self.obs.frequencies_hz, self.obs.values, self.obs.acquisition, self.obs.band_index)
        with self.assertRaises(ValueError):
            PatternObjective(Parameterization.from_interpretation(self.far, POLICY), bare)


if __name__ == "__main__":
    unittest.main()
