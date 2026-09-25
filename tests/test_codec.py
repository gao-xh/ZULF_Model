import math
import unittest

import numpy as np

from zulf_model.codec import BOS, EOS, GrammarState, InterpretationCodec, JBinner
from zulf_model.generator import build_default_sampler
from zulf_model.spec import JBinSpec, ProblemSpec
from zulf_model.spinsystem import Component, Interpretation, SpinSystem, best_permutation


def matched_error(a: SpinSystem, b: SpinSystem) -> float:
    return best_permutation(a, b).max_abs_error_hz


class CodecTests(unittest.TestCase):
    def setUp(self):
        self.spec = ProblemSpec()
        self.codec = InterpretationCodec(self.spec)
        self.samples = list(build_default_sampler(self.spec).generate(40, seed=11))

    def test_j_binner_roundtrip(self):
        binner = JBinner(self.spec)
        for value in (-297.3, -12.4, -0.1, 0.0, 0.13, 7.77, 42.2, 135.0, 250.3):
            index, offset = binner.encode(value)
            self.assertAlmostEqual(binner.decode(index, offset), value, places=9)
            self.assertLessEqual(abs(binner.decode(index) - value), binner.widths[index] / 2 + 1e-12)

    def test_sequence_roundtrip_within_bin_precision(self):
        for sample in self.samples:
            enc = self.codec.encode(sample.interpretation)
            self.assertEqual(enc.tokens[0], self.codec.vocab[BOS])
            self.assertEqual(enc.tokens[-1], self.codec.vocab[EOS])
            self.assertLessEqual(len(enc.tokens), self.codec.grammar.max_length())
            exact = self.codec.decode(enc.tokens, enc.j_offsets)
            coarse = self.codec.decode(enc.tokens)
            truth = sorted(sample.interpretation.components, key=lambda c: -c.contribution)
            for t, e, c in zip(truth, exact.components, coarse.components):
                self.assertLess(matched_error(t.system, e.system), 1e-6)
                self.assertLessEqual(matched_error(t.system, c.system), 1.0 + 1e-9)

    def test_grammar_masks_follow_state(self):
        g = self.codec.grammar
        v = self.codec.vocab
        state = GrammarState()
        self.assertEqual(np.flatnonzero(g.allowed(state)).tolist(), [v[BOS]])
        enc = self.codec.encode(self.samples[0].interpretation)
        for token in enc.tokens:
            self.assertTrue(g.allowed(state)[token])
            state = g.advance(state, token)
        self.assertEqual(state.phase, "done")
        self.assertFalse(g.validate(enc.tokens[:-1]))

    def test_grammar_rejects_signal_free_and_wrong_size(self):
        v = self.codec.vocab
        tokens = [v[BOS], v["<comp>"], v.w_start, v.group_tokens[("1H", 6)], v.group_tokens[("1H", 2)], v["<sep>"]]
        self.assertFalse(self.codec.grammar.validate(tokens + [v.j_start, v[EOS]]))

    def test_set_roundtrip(self):
        for sample in self.samples[:15]:
            arrays = self.codec.encode_set(sample.interpretation)
            decoded = self.codec.decode_set(arrays)
            truth = sorted(sample.interpretation.components, key=lambda c: -c.contribution)
            self.assertEqual(len(decoded), min(len(truth), self.spec.max_components))
            for t, d in zip(truth, decoded.components):
                self.assertLess(matched_error(t.system, d.system), 1e-4)

    def test_dimensions_follow_spec(self):
        spec = ProblemSpec(spin_counts=(6, 10), max_components=2, max_group_size=5,
                           j_bins=JBinSpec(((10.0, 0.5), (300.0, 2.0))))
        codec = InterpretationCodec(spec)
        self.assertEqual(len(codec.j), 2 * (20 + 145))
        interp = Interpretation((Component(SpinSystem(("13C",) + ("1H",) * 5,
                                                      np.pad(np.full((1, 5), 140.0), ((0, 5), (1, 0))) +
                                                      np.pad(np.full((1, 5), 140.0), ((0, 5), (1, 0))).T)),))
        arrays = codec.encode_set(interp)
        self.assertEqual(arrays["couplings"].shape, (2, 10, 10))
        self.assertTrue(codec.grammar.validate(codec.encode(interp).tokens))


if __name__ == "__main__":
    unittest.main()
