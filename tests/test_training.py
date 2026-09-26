import json
import tempfile
import unittest

import numpy as np

from zulf_core.evaluation.matching import match_interpretations
from zulf_model.generator import build_default_sampler
from zulf_model.spec import ProblemSpec
from zulf_core.spinsystem import Component, Interpretation, SpinSystem
from zulf_model.training import CoverageMetrics, TrainingSetup


def tiny_run_config(output_dir, steps=12, kind="cnn_set", curriculum=()):
    return {
        "problem": {"spin_counts": [8], "max_components": 3},
        "processing": {"points": 2048, "sampling_rate_hz": 1000.0},
        "perturbation": {"snr_range": [20, 200]},
        "model": {"kind": kind,
                  "encoder": {"d_model": 32, "stages": [[16, 9, 4], [32, 7, 4]], "transformer_layers": 1,
                              "heads": 2, "feedforward": 64},
                  "set_head": {"decoder_layers": 1, "heads": 2, "feedforward": 64},
                  "sequence_head": {"decoder_layers": 1, "heads": 2, "feedforward": 64}},
        "train": {"steps": steps, "batch_size": 4, "warmup_steps": 2, "log_every": 4, "eval_every": 6,
                  "checkpoint_every": 6, "eval_propose_k": 3, "eval_propose_limit": 4, "device": "cpu",
                  "output_dir": output_dir, "renders_per_system": 2, "curriculum": list(curriculum)},
        "validation": {"count": 6, "batch_size": 3, "seed": 5},
    }


class MatchingTests(unittest.TestCase):
    def test_structure_and_j_match(self):
        j = np.zeros((4, 4)); j[0, 1:] = j[1:, 0] = 140.0; j[1, 2] = j[2, 1] = 0.0
        truth = Interpretation((Component(SpinSystem(("13C", "1H", "1H", "1H"), j)),))
        near = j.copy(); near[0, 1:] += 0.4; near[1:, 0] += 0.4
        cand = Interpretation((Component(SpinSystem(("1H", "13C", "1H", "1H"), near[[1, 0, 2, 3]][:, [1, 0, 2, 3]])),))
        m = match_interpretations(truth, cand, tolerance_hz=0.5)
        self.assertTrue(m.structure and m.j_match)
        self.assertFalse(match_interpretations(truth, cand, tolerance_hz=0.1).j_match)
        flipped = Interpretation((Component(SpinSystem(truth.components[0].system.isotopes, -j)),))
        self.assertTrue(match_interpretations(truth, flipped).j_match)

    def test_coverage_metrics(self):
        spec = ProblemSpec()
        samples = list(build_default_sampler(spec).generate(5, seed=2))
        metrics = CoverageMetrics(ks=(1, 3))
        for s in samples:
            metrics.update(s.interpretation, [samples[0].interpretation, s.interpretation])
        summary = metrics.summary()
        self.assertEqual(summary["j_coverage@3"], 1.0)
        self.assertGreaterEqual(summary["structure_coverage@1"], 0.2)


class TrainerTests(unittest.TestCase):
    def test_short_runs_log_checkpoint_and_resume(self):
        for kind in ("cnn_set", "cnn_transformer"):
            with tempfile.TemporaryDirectory() as tmp:
                curriculum = [{"until_step": 6, "name": "clean", "overrides": {"perturbation": {"snr_range": [500, 500]}}}]
                setup = TrainingSetup(tiny_run_config(tmp, kind=kind, curriculum=curriculum))
                trainer = setup.build_trainer()
                result = trainer.fit()
                self.assertEqual(result["steps"], 12)
                self.assertIn("coverage", result["last_eval"])
                with open(f"{tmp}/log.jsonl", encoding="utf-8") as handle:
                    events = [json.loads(line)["event"] for line in handle]
                self.assertIn("curriculum_stage", events)
                self.assertIn("eval", events)
                self.assertIn("train.step", result["timing"])
                again = setup.build_trainer()
                again.resume(f"{tmp}/last.pt")
                self.assertEqual(again.step, 12)
                again.config = again.config.__class__(**dict(again.config.__dict__, steps=14))
                self.assertEqual(again.fit()["steps"], 14)


if __name__ == "__main__":
    unittest.main()
