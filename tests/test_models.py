import tempfile
import unittest

import numpy as np
import torch

from zulf_model.codec import InterpretationCodec
from zulf_model.generator import build_default_sampler
from zulf_model.models import ModelConfig, build_model, load_model, save_model
from zulf_model.models.config import EncoderConfig, SequenceHeadConfig, SetHeadConfig
from zulf_model.render import PerturbationConfig, ProcessingConfig, SampleRenderer
from zulf_model.spec import ProblemSpec
from zulf_model.training.data import Collator, FixedDataset

TINY_ENCODER = EncoderConfig(d_model=32, stages=((16, 9, 4), (32, 7, 4)), transformer_layers=1, heads=2,
                             feedforward=64, dropout=0.0)


def tiny(kind):
    return ModelConfig(kind=kind, encoder=TINY_ENCODER, set_head=SetHeadConfig(decoder_layers=1, heads=2, feedforward=64),
                       sequence_head=SequenceHeadConfig(decoder_layers=1, heads=2, feedforward=64))


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.manual_seed(0)
        cls.spec = ProblemSpec()
        cls.codec = InterpretationCodec(cls.spec)
        samples = list(build_default_sampler(cls.spec).generate(4, seed=3))
        cls.renderer = SampleRenderer(cls.spec, ProcessingConfig(points=2048), PerturbationConfig(snr_range=(50, 50)))
        data = FixedDataset(samples, cls.renderer, cls.codec)
        cls.batch = Collator(cls.codec, cls.renderer.grid.frequencies_hz)([data[i] for i in range(len(data))])

    def check_valid(self, proposals):
        for candidates in proposals:
            self.assertGreater(len(candidates), 0)
            for interp in candidates:
                for c in interp.components:
                    self.assertTrue(self.spec.allows(c.system.isotopes))
                    self.assertGreaterEqual(len(set(c.system.isotopes)), 2)

    def test_forward_loss_backward_and_propose(self):
        for kind in ("cnn_set", "cnn_transformer"):
            model = build_model(self.spec, tiny(kind))
            out = model(self.batch)
            loss, logs = model.loss(out, self.batch)
            loss.backward()
            self.assertTrue(np.isfinite(logs["loss"]))
            self.assertTrue(all(p.grad is not None for p in model.parameters() if p.requires_grad))
            self.check_valid(model.propose(self.batch["features"][:2], self.batch["frequency_hz"], k=3))

    def test_checkpoint_roundtrip(self):
        for kind in ("cnn_set", "cnn_transformer"):
            model = build_model(self.spec, tiny(kind)).eval()
            with tempfile.TemporaryDirectory() as tmp:
                save_model(f"{tmp}/m.pt", model)
                again = load_model(f"{tmp}/m.pt").eval()
            a = model(self.batch)
            b = again(self.batch)
            for key in a:
                torch.testing.assert_close(a[key], b[key])

    def test_overfit_small_batch(self):
        for kind in ("cnn_set", "cnn_transformer"):
            torch.manual_seed(1)
            model = build_model(self.spec, tiny(kind))
            opt = torch.optim.AdamW(model.parameters(), lr=3e-3)
            first = None
            for step in range(60):
                loss, _ = model.loss(model(self.batch), self.batch)
                opt.zero_grad()
                loss.backward()
                opt.step()
                first = loss.item() if first is None else first
            self.assertLess(loss.item(), 0.6 * first, msg=kind)

    def test_beam_sequences_follow_grammar(self):
        model = build_model(self.spec, tiny("cnn_transformer")).eval()
        memory = model.encode(self.batch["features"][:1], self.batch["frequency_hz"])
        hyps = model.beam_search(memory, beam_size=3)
        self.assertTrue(hyps)
        for h in hyps:
            self.assertTrue(model.grammar.validate(h.tokens))
        samples = model.sample(memory, 3, generator=torch.Generator().manual_seed(0))
        for h in samples:
            self.assertTrue(model.grammar.validate(h.tokens))


if __name__ == "__main__":
    unittest.main()
