"""CNN + Transformer encoder-decoder over codec token sequences.

Training: teacher forcing with label smoothing; J tokens use soft targets that
spread probability over neighbouring bins (Gaussian, `j_neighbor_sigma_bins`),
so adjacent-bin errors are penalized less than distant ones. A regression head
predicts the within-bin J offset. Decoding: beam search or temperature sampling
under the codec grammar, so every returned sequence decodes to a valid
interpretation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from torch import nn

from ..codec import BOS, PAD, GrammarState
from ..spinsystem import Interpretation
from .base import CandidateModel
from .config import ModelConfig
from .encoder import SpectrumEncoder


@dataclass
class Hypothesis:
    tokens: List[int]
    offsets: List[float]
    state: GrammarState
    logp: float

    @property
    def done(self) -> bool:
        return self.state.phase == "done"


class CNNTransformerModel(CandidateModel):
    kind = "cnn_transformer"

    def __init__(self, spec, config: ModelConfig):
        super().__init__(spec, config)
        enc, head = config.encoder, config.sequence_head
        d = enc.d_model
        self.vocab = self.codec.vocab
        self.grammar = self.codec.grammar
        self.max_len = self.grammar.max_length()
        self.encoder = SpectrumEncoder(enc)
        self.embed = nn.Embedding(len(self.vocab), d, padding_idx=self.vocab[PAD])
        self.position = nn.Embedding(self.max_len + 1, d)
        layer = nn.TransformerDecoderLayer(d, head.heads, head.feedforward, enc.dropout, batch_first=True,
                                           norm_first=True)
        self.decoder = nn.TransformerDecoder(layer, head.decoder_layers)
        self.norm = nn.LayerNorm(d)
        self.logits = nn.Linear(d, len(self.vocab))
        self.offset = nn.Linear(d, 1)
        j = torch.arange(self.vocab.j_bins, dtype=torch.float32)
        self.register_buffer("j_index", j, persistent=False)

    # -- core ---------------------------------------------------------------------
    def encode(self, features: torch.Tensor, frequency_hz: torch.Tensor) -> torch.Tensor:
        return self.encoder(features, frequency_hz)

    def decode_step(self, prefix: torch.Tensor, memory: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """prefix: (B, L) -> logits (B, L, V), offsets (B, L)."""
        length = prefix.shape[1]
        positions = torch.arange(length, device=prefix.device)
        x = self.embed(prefix) + self.position(positions).unsqueeze(0)
        causal = torch.triu(torch.ones(length, length, dtype=torch.bool, device=prefix.device), 1)
        pad = prefix == self.vocab[PAD]
        h = self.norm(self.decoder(x, memory, tgt_mask=causal, tgt_key_padding_mask=pad))
        return self.logits(h), 0.5 * torch.tanh(self.offset(h).squeeze(-1))

    def forward(self, batch: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        memory = self.encode(batch["features"], batch["frequency_hz"])
        logits, offsets = self.decode_step(batch["tokens"][:, :-1], memory)
        return {"logits": logits, "offsets": offsets}

    # -- loss ---------------------------------------------------------------------
    def soft_targets(self, target: torch.Tensor) -> torch.Tensor:
        head = self.model_config.sequence_head
        v = len(self.vocab)
        eps = head.label_smoothing
        q = torch.full(target.shape + (v,), eps / v, device=target.device)
        q.scatter_(-1, target.unsqueeze(-1), 1 - eps + eps / v)
        is_j = (target >= self.vocab.j_start) & (target < self.vocab.j_start + self.vocab.j_bins)
        if head.j_neighbor_sigma_bins > 0 and is_j.any():
            centre = (target[is_j] - self.vocab.j_start).float()
            weights = torch.exp(-0.5 * ((self.j_index.unsqueeze(0) - centre.unsqueeze(1)) / head.j_neighbor_sigma_bins) ** 2)
            weights = weights / weights.sum(-1, keepdim=True)
            rows = torch.full((int(is_j.sum()), v), eps / v, device=target.device)
            rows[:, self.vocab.j_start:self.vocab.j_start + self.vocab.j_bins] += (1 - eps) * weights
            q[is_j] = rows
        return q

    def loss(self, out, batch):
        head = self.model_config.sequence_head
        target = batch["tokens"][:, 1:]
        valid = target != self.vocab[PAD]
        logp = out["logits"].float().log_softmax(-1)
        q = self.soft_targets(target)
        ce = -(q * logp).sum(-1)
        ce = (ce * valid).sum() / valid.sum().clamp(min=1)
        mask = batch["j_offset_mask"][:, 1:]
        off = ((out["offsets"].float() - batch["j_offsets"][:, 1:]) ** 2 * mask).sum() / mask.sum().clamp(min=1)
        total = ce + head.weight_offset * off
        with torch.no_grad():
            n_valid = valid.sum().clamp(min=1)
            # Entropy of the soft targets: the floor of token_ce; token_ce - floor is the KL that can be learned.
            floor = (-(q * q.clamp(min=1e-30).log()).sum(-1) * valid).sum() / n_valid
            predicted = out["logits"].argmax(-1)
            accuracy = ((predicted == target) & valid).sum() / n_valid
            is_j = (target >= self.vocab.j_start) & (target < self.vocab.j_start + self.vocab.j_bins) & valid
            structural = valid & ~is_j
            structure_accuracy = ((predicted == target) & structural).sum() / structural.sum().clamp(min=1)
            j_near = ((predicted - target).abs() <= 1) & is_j
            j_accuracy_1bin = j_near.sum() / is_j.sum().clamp(min=1)
        return total, {"loss": total.item(), "token_ce": ce.item(), "token_ce_floor": floor.item(),
                       "token_kl": (ce - floor).item(), "offset_mse": off.item(),
                       "token_accuracy": accuracy.item(), "structure_token_accuracy": structure_accuracy.item(),
                       "j_token_accuracy_1bin": j_accuracy_1bin.item()}

    # -- decoding -----------------------------------------------------------------
    def _grammar_mask(self, states: List[GrammarState], device) -> torch.Tensor:
        return torch.as_tensor(np.stack([self.grammar.allowed(s) for s in states]), device=device)

    @torch.no_grad()
    def beam_search(self, memory: torch.Tensor, beam_size: int = 8, length_penalty: float = 0.6) -> List[Hypothesis]:
        """memory: (1, T, D). Returns finished hypotheses sorted by normalized score."""
        bos = self.vocab[BOS]
        start = self.grammar.advance(GrammarState(), bos)
        beams = [Hypothesis([bos], [math.nan], start, 0.0)]
        finished: List[Hypothesis] = []
        for _ in range(self.max_len):
            if not beams:
                break
            prefix = torch.as_tensor([h.tokens for h in beams], device=memory.device)
            logits, offsets = self.decode_step(prefix, memory.expand(len(beams), -1, -1))
            logp = logits[:, -1].float().log_softmax(-1)
            logp = logp.masked_fill(~self._grammar_mask([h.state for h in beams], memory.device), float("-inf"))
            candidates = []
            width = min(beam_size, logp.shape[-1])
            top_values, top_index = logp.topk(width, dim=-1)
            for b, hyp in enumerate(beams):
                for value, index in zip(top_values[b].tolist(), top_index[b].tolist()):
                    if not math.isfinite(value):
                        continue
                    offset = float(offsets[b, -1]) if self.vocab.is_j(index) else math.nan
                    candidates.append(Hypothesis(hyp.tokens + [index], hyp.offsets + [offset],
                                                 self.grammar.advance(hyp.state, index), hyp.logp + value))
            candidates.sort(key=lambda h: -h.logp)
            beams = []
            for hyp in candidates:
                (finished if hyp.done else beams).append(hyp)
                if len(beams) >= beam_size:
                    break
            if len(finished) >= beam_size and beams and max(h.logp for h in beams) < min(
                    f.logp for f in sorted(finished, key=lambda h: -h.logp)[:beam_size]):
                break
        finished.sort(key=lambda h: -h.logp / (len(h.tokens) ** length_penalty))
        return finished

    @torch.no_grad()
    def sample(self, memory: torch.Tensor, count: int, temperature: float = 1.0,
               generator: Optional[torch.Generator] = None) -> List[Hypothesis]:
        bos = self.vocab[BOS]
        hyps = [Hypothesis([bos], [math.nan], self.grammar.advance(GrammarState(), bos), 0.0) for _ in range(count)]
        for _ in range(self.max_len):
            active = [h for h in hyps if not h.done and h.state.phase != "dead"]
            if not active:
                break
            prefix = torch.as_tensor([h.tokens for h in active], device=memory.device)
            logits, offsets = self.decode_step(prefix, memory.expand(len(active), -1, -1))
            logp = (logits[:, -1].float() / temperature).log_softmax(-1)
            logp = logp.masked_fill(~self._grammar_mask([h.state for h in active], memory.device), float("-inf"))
            probs = logp.exp()
            dead = probs.sum(-1) <= 0
            probs[dead] = 1.0  # placeholder; dead hypotheses are discarded below
            choice = torch.multinomial(probs, 1, generator=generator).squeeze(-1)
            for b, hyp in enumerate(active):
                if dead[b]:
                    hyp.state = GrammarState(phase="dead")
                    continue
                index = int(choice[b])
                hyp.tokens.append(index)
                hyp.offsets.append(float(offsets[b, -1]) if self.vocab.is_j(index) else math.nan)
                hyp.logp += float(logp[b, index])
                hyp.state = self.grammar.advance(hyp.state, index)
        return [h for h in hyps if h.done]

    @torch.no_grad()
    def propose(self, features, frequency_hz, k: int = 5, beam_size: Optional[int] = None,
                samples: int = 0, temperature: float = 1.0) -> List[List[Interpretation]]:
        was_training = self.training
        self.eval()
        memory = self.encode(features, frequency_hz)
        results = []
        for b in range(features.shape[0]):
            hyps = self.beam_search(memory[b:b + 1], beam_size or max(k, 4))
            if samples:
                hyps += self.sample(memory[b:b + 1], samples, temperature)
            candidates, seen = [], set()
            for h in hyps:
                try:
                    interp = self.codec.decode(h.tokens, h.offsets, score=h.logp)
                except ValueError:
                    continue
                interp.metadata["source"] = "cnn_transformer"
                key = interp.key(self.spec.nuclei)
                if key not in seen:
                    seen.add(key)
                    candidates.append(interp)
            results.append(candidates[:k])
        if was_training:
            self.train()
        return results
