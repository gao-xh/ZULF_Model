import json
import tempfile
import unittest

import numpy as np

from zulf_model.evaluation import (BenchmarkConfig, CandidateProposer, PriorSearchProposer, RandomPriorProposer,
                                   basin_of_attraction, local_identifiability, quick_score, run_benchmark)
from zulf_core.evaluation.identifiability import perturb_couplings
from zulf_model.finetune import ActiveLearningLoop, FocusedSource, LoopConfig, classify_failure
from zulf_model.generator import Sample, SplitConfig, build_default_sampler, split_of
from zulf_core.physics import TransitionCache, compute_transitions
from zulf_model.render import Acquisition, PerturbationConfig, Renderer
from zulf_core.solver import ObservedSpectrum, RefineSettings
from zulf_model.spec import ProblemSpec
from zulf_core.spinsystem import Component, Interpretation, SpinSystem

ACQ = Acquisition(1000.0, 3000)
RANGES = [(100.0, 150.0), (230.0, 275.0)]


def methyl():
    iso = ("1H",) * 7 + ("13C",)
    j = np.zeros((8, 8))
    j[:6, 6] = j[6, :6] = 6.5
    j[:3, 7] = j[7, :3] = 126.0
    j[3:6, 7] = j[7, 3:6] = 4.0
    j[6, 7] = j[7, 6] = -3.0
    return SpinSystem(iso, j)


def observe(system, gain=1.0):
    fid = Renderer(ACQ).synthesize(compute_transitions(system), 2.0, gain=gain)
    return ObservedSpectrum.from_fid(fid, ACQ, RANGES)


class OracleProposer(CandidateProposer):
    name = "oracle"

    def __init__(self, truth, jitter=0.2):
        self.truth, self.jitter = truth, jitter

    def propose(self, observed, k):
        rng = np.random.default_rng(0)
        return [Interpretation(tuple(Component(perturb_couplings(c.system, self.jitter, rng), c.contribution)
                                     for c in self.truth.components)) for _ in range(min(k, 2))]


class IdentifiabilityTests(unittest.TestCase):
    def test_report_and_basin(self):
        truth = Interpretation((Component(methyl()),))
        obs = observe(methyl())
        report = local_identifiability(truth, obs)
        self.assertEqual(len(report.singular_values), len(report.parameter_names))
        self.assertTrue(all(np.diff(report.singular_values) <= 1e-12))
        self.assertAlmostEqual(float(np.linalg.norm(list(report.weakest_direction.values()))), 1.0, places=6)
        rows = basin_of_attraction(truth, obs, [0.3], trials=1,
                                   settings=RefineSettings(starts=1, continuation_rates_per_s=(3.0, 0.0)))
        self.assertEqual(rows[0]["success_rate"], 1.0)


class ProposerTests(unittest.TestCase):
    def test_quick_score_and_prior_search(self):
        truth = Interpretation((Component(methyl()),))
        obs = observe(methyl(), gain=0.7j)
        cache = TransitionCache()
        self.assertLess(quick_score(truth, obs, cache, rate=2.0), 1e-8)
        wrong = Interpretation((Component(perturb_couplings(methyl(), 3.0, np.random.default_rng(1))),))
        self.assertGreater(quick_score(wrong, obs, cache, rate=2.0), 0.1)
        spec = ProblemSpec()
        sampler = build_default_sampler(spec)
        candidates = PriorSearchProposer(sampler, pool=8).propose(obs, 3)
        self.assertEqual(len(candidates), 3)
        self.assertGreaterEqual(candidates[0].score, candidates[-1].score)
        self.assertEqual(len(RandomPriorProposer(sampler).propose(obs, 4)), 4)


class BenchmarkTests(unittest.TestCase):
    def test_oracle_beats_random(self):
        spec = ProblemSpec()
        truth = Interpretation((Component(methyl()),))
        sample = Sample(truth, "fam", "manual", [0])
        perturbation = PerturbationConfig(snr_range=(1e4, 1e4), family_split_probability=0.0, gaussian_probability=0.0,
                                          rate_range_per_s=(2.0, 2.0), drift_probability=0.0,
                                          interference_probability=0.0, phase_delay_range_s=(0.0, 0.0))
        with tempfile.TemporaryDirectory() as tmp:
            report = run_benchmark([sample], [OracleProposer(truth), RandomPriorProposer(build_default_sampler(spec))],
                                   ACQ, perturbation, BenchmarkConfig(k=2, held_out=1, ranges=tuple(RANGES),
                                                                      tolerance_hz=0.05),
                                   RefineSettings(starts=1, continuation_rates_per_s=(3.0, 0.0)), f"{tmp}/b.json")
            with open(f"{tmp}/b.json", encoding="utf-8") as handle:
                saved = json.load(handle)
        self.assertEqual(report["summary"]["oracle"]["j_hit_post"], 1.0)
        self.assertEqual(report["summary"]["random_prior"]["j_hit_post"], 0.0)
        self.assertIn("rows", saved)


class FineTuneTests(unittest.TestCase):
    def test_classification_and_focused_source(self):
        spec = ProblemSpec()
        truth = Interpretation((Component(methyl()),))
        sample = Sample(truth, "fam-dev", "manual", [0])
        near = Interpretation((Component(perturb_couplings(methyl(), 0.01, np.random.default_rng(0))),))
        far = Interpretation((Component(perturb_couplings(methyl(), 5.0, np.random.default_rng(0))),))
        other = Interpretation((Component(SpinSystem(("13C", "1H"), np.array([[0, 140.0], [140.0, 0]]))),))
        self.assertEqual(classify_failure(sample, [near], tolerance_hz=0.1).category, "success")
        self.assertEqual(classify_failure(sample, [far], tolerance_hz=0.1).category, "j_outside_basin")
        self.assertEqual(classify_failure(sample, [other]).category, "structure_miss")
        config = SplitConfig(val_fraction=0.0, test_unseen_fraction=0.0, test_seen_fraction=0.0)
        failure = classify_failure(sample, [far], tolerance_hz=0.1)
        source = FocusedSource(spec, [failure], jitter_hz=0.5, split_config=config)
        drawn = source.draw(np.random.default_rng(1))
        self.assertEqual(drawn.family_id, "fam-dev")
        self.assertEqual(drawn.interpretation.components[0].system.group_signature(), methyl().group_signature())
        frozen = SplitConfig(val_fraction=0.0, test_unseen_fraction=1.0)
        with self.assertRaises(ValueError):
            FocusedSource(spec, [failure], split_config=frozen)

    def test_loop_rounds(self):
        spec = ProblemSpec()
        base = build_default_sampler(spec, {"splits": {"val_fraction": 0.0, "test_unseen_fraction": 0.0,
                                                       "test_seen_fraction": 0.0}})
        truth = Interpretation((Component(methyl()),))
        dev = [Sample(truth, "fam-a", "manual", [0])]
        calls = []

        class Wrong(CandidateProposer):
            def propose(self, observed, k):
                return [Interpretation((Component(perturb_couplings(methyl(), 5.0, np.random.default_rng(3))),))]

        with tempfile.TemporaryDirectory() as tmp:
            loop = ActiveLearningLoop(dev, Wrong(), lambda s: (observe(methyl()), []),
                                      lambda sampler, steps, lr: calls.append((sampler, steps, lr)), base,
                                      LoopConfig(rounds=2, refine=False, output_dir=tmp))
            history = loop.run()
        self.assertEqual(len(history), 2)
        self.assertEqual(history[0]["failures"], {"j_outside_basin": 1})
        self.assertEqual(len(calls), 1)
        sampler = calls[0][0]
        self.assertIn("focused", sampler.sources)
        self.assertGreater(sampler.config.source_weights["focused"], 0)


class BlockSignMatchingTests(unittest.TestCase):
    @staticmethod
    def system(j_nh, linked):
        from zulf_core.spinsystem import SpinSystem
        j = np.zeros((4, 4))
        values = {(0, 2): 140.0, (1, 3): j_nh}
        if linked:
            values.update({(0, 1): -5.0, (0, 3): 3.0, (1, 2): -2.0, (2, 3): 5.0})
        for (a, b), v in values.items():
            j[a, b] = j[b, a] = v
        return SpinSystem(("13C", "15N", "1H", "1H"), j)

    def test_relative_sign_free_only_between_disconnected_blocks(self):
        from zulf_core.evaluation.matching import coupling_blocks, coupling_error
        from zulf_core.spinsystem import SpinSystem
        self.assertEqual(len(coupling_blocks(self.system(-70.0, False))), 2)
        self.assertEqual(len(coupling_blocks(self.system(-70.0, True))), 1)
        self.assertLess(coupling_error(self.system(-70.0, False), self.system(70.0, False))[0], 1e-9)
        self.assertAlmostEqual(coupling_error(self.system(-70.0, True), self.system(70.0, True))[0], 140.0)
        truth = self.system(-70.0, True)
        flipped = SpinSystem(truth.isotopes, -truth.couplings_hz)
        self.assertLess(coupling_error(truth, flipped)[0], 1e-9)


if __name__ == "__main__":
    unittest.main()
