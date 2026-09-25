"""CNN set-prediction baseline (DETR style), at the level of equivalence groups.

Learned component slots cross-attend to encoder tokens. Each slot predicts
existence, a (nucleus, group size) class for each of G group slots (with a
padding class), group-pair J means and log-variances, and a relative
contribution. Groups follow the codec's canonical group order, so magnetic
equivalence is part of the prediction rather than recovered afterwards.
Training uses Hungarian matching between slots and ground-truth components.
Decoding chooses, per slot, the most probable group sequence whose spin count
is allowed by the spec, with at least two nuclei and canonical nucleus order.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from scipy.optimize import linear_sum_assignment
from torch import nn

from ..codec import InterpretationCodec
from ..spinsystem import Component, Interpretation, SpinSystem
from .base import CandidateModel
from .config import ModelConfig
from .encoder import SpectrumEncoder


class CNNSetModel(CandidateModel):
    kind = "cnn_set"

    def __init__(self, spec, config: ModelConfig):
        super().__init__(spec, config)
        enc, head = config.encoder, config.set_head
        d = enc.d_model
        self.codec = InterpretationCodec(spec)
        self.k_slots, self.g_max = spec.max_components, spec.max_spins
        self.pad = self.codec.n_group_classes
        self.n_classes = self.pad + 1
        iu = torch.triu_indices(self.g_max, self.g_max, 1)
        self.register_buffer("pair_rows", iu[0], persistent=False)
        self.register_buffer("pair_cols", iu[1], persistent=False)
        self.n_pairs = iu.shape[1]
        self.encoder = SpectrumEncoder(enc)
        self.queries = nn.Parameter(torch.randn(self.k_slots, d) * 0.02)
        layer = nn.TransformerDecoderLayer(d, head.heads, head.feedforward, enc.dropout, batch_first=True,
                                           norm_first=True)
        self.decoder = nn.TransformerDecoder(layer, head.decoder_layers)
        self.exists = nn.Linear(d, 1)
        self.group = nn.Linear(d, self.g_max * self.n_classes)
        self.j_mean = nn.Linear(d, self.n_pairs)
        self.j_logvar = nn.Linear(d, self.n_pairs)
        self.contribution = nn.Linear(d, 1)

    def forward(self, batch: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        tokens = self.encoder(batch["features"], batch["frequency_hz"])
        b = tokens.shape[0]
        slots = self.decoder(self.queries.unsqueeze(0).expand(b, -1, -1), tokens)
        return {"exists": self.exists(slots).squeeze(-1),
                "group": self.group(slots).view(b, self.k_slots, self.g_max, self.n_classes),
                "j_mean": self.j_mean(slots),
                "j_logvar": self.j_logvar(slots).clamp(-10, 6),
                "contribution": self.contribution(slots).squeeze(-1)}

    # -- loss -------------------------------------------------------------------
    def _targets(self, batch):
        couplings = batch["group_couplings"] / self.model_config.set_head.j_scale_hz
        j = couplings[..., self.pair_rows, self.pair_cols]
        mask = batch["group_mask"]
        pair_mask = mask[..., self.pair_rows] * mask[..., self.pair_cols]
        return j, pair_mask

    def _pairwise_cost(self, out, batch, b, j_target, pair_mask):
        """Cost (slots x true components) for sample b."""
        true = torch.nonzero(batch["component_mask"][b] > 0.5).flatten()
        if len(true) == 0:
            return None, true
        logp_exist = nn.functional.logsigmoid(out["exists"][b])                        # (K,)
        logp_group = out["group"][b].log_softmax(-1)                                   # (K,G,C)
        target = batch["group_class"][b, true]                                         # (T,G)
        group_cost = -torch.gather(logp_group.unsqueeze(1).expand(-1, len(true), -1, -1), 3,
                                   target.unsqueeze(0).unsqueeze(-1).expand(self.k_slots, -1, -1, 1)).squeeze(-1).mean(-1)
        mu, logvar = out["j_mean"][b].unsqueeze(1), out["j_logvar"][b].unsqueeze(1)
        jt, pm = j_target[b, true].unsqueeze(0), pair_mask[b, true].unsqueeze(0)
        nll = 0.5 * (logvar + (jt - mu) ** 2 / logvar.exp())
        j_cost = (nll * pm).sum(-1) / pm.sum(-1).clamp(min=1)
        c_cost = (out["contribution"][b].unsqueeze(1) - batch["log10_contribution"][b, true].unsqueeze(0)).abs()
        h = self.model_config.set_head
        cost = (-h.weight_exists * logp_exist.unsqueeze(1) + h.weight_isotope * group_cost + h.weight_j * j_cost
                + h.weight_contribution * c_cost)
        return cost, true

    def loss(self, out, batch) -> Tuple[torch.Tensor, Dict[str, float]]:
        h = self.model_config.set_head
        j_target, pair_mask = self._targets(batch)
        b_size = out["exists"].shape[0]
        exist_target = torch.zeros_like(out["exists"])
        group_terms, j_terms, c_terms = [], [], []
        for b in range(b_size):
            cost, true = self._pairwise_cost(out, batch, b, j_target, pair_mask)
            if cost is None:
                continue
            rows, cols = linear_sum_assignment(cost.detach().float().cpu().numpy())
            rows_t = torch.as_tensor(rows, device=cost.device)
            cols_t = true[torch.as_tensor(cols, device=cost.device)]
            exist_target[b, rows_t] = 1.0
            group_terms.append(nn.functional.cross_entropy(out["group"][b, rows_t].reshape(-1, self.n_classes),
                                                           batch["group_class"][b, cols_t].reshape(-1)))
            mu, logvar = out["j_mean"][b, rows_t], out["j_logvar"][b, rows_t]
            pm = pair_mask[b, cols_t]
            nll = 0.5 * (logvar + (j_target[b, cols_t] - mu) ** 2 / logvar.exp())
            j_terms.append((nll * pm).sum() / pm.sum().clamp(min=1))
            c_terms.append((out["contribution"][b, rows_t] - batch["log10_contribution"][b, cols_t]).abs().mean())
        exist_loss = nn.functional.binary_cross_entropy_with_logits(out["exists"], exist_target)
        zero = out["exists"].sum() * 0
        group_loss = torch.stack(group_terms).mean() if group_terms else zero
        j_loss = torch.stack(j_terms).mean() if j_terms else zero
        c_loss = torch.stack(c_terms).mean() if c_terms else zero
        total = (h.weight_exists * exist_loss + h.weight_isotope * group_loss + h.weight_j * j_loss
                 + h.weight_contribution * c_loss)
        logs = {"loss": total.item(), "exists": exist_loss.item(), "group": group_loss.item(),
                "j_nll": j_loss.item(), "contribution": c_loss.item()}
        return total, logs

    # -- inference -------------------------------------------------------------
    def constrained_groups(self, logp: np.ndarray) -> Optional[List[int]]:
        """Most probable valid group-class sequence for one slot (dynamic programming).

        Valid: groups packed at the front (padding only after the last group),
        nucleus order non-decreasing (canonical order), total spins in
        spec.spin_counts, at least two distinct nuclei. Returns the class list of
        the groups (without padding) or None when nothing valid exists.
        """
        allowed = set(self.spec.spin_counts)
        n_nuclei = len(self.spec.nuclei)
        size_max = self.spec.max_group_size
        # state: (spins, nucleus bit mask, last nucleus) -> (score, classes)
        states = {(0, 0, -1): (0.0, [])}
        best, best_score = None, -np.inf
        for position in range(self.g_max):
            pad_tail = float(logp[position:, self.pad].sum())
            for (spins, mask, last), (score, classes) in states.items():
                if spins in allowed and bin(mask).count("1") >= 2 and score + pad_tail > best_score:
                    best_score, best = score + pad_tail, classes
            nxt: Dict[tuple, tuple] = {}
            for (spins, mask, last), (score, classes) in states.items():
                for c in range(self.pad):
                    nucleus, size = divmod(c, size_max)
                    size += 1
                    if nucleus < last or spins + size > self.g_max:
                        continue
                    key = (spins + size, mask | (1 << nucleus), nucleus)
                    value = score + float(logp[position, c])
                    if key not in nxt or value > nxt[key][0]:
                        nxt[key] = (value, classes + [c])
            states = nxt
        for (spins, mask, last), (score, classes) in states.items():
            if spins in allowed and bin(mask).count("1") >= 2 and score > best_score:
                best_score, best = score, classes
        return best

    def _slot_system(self, classes: List[int], j_upper: np.ndarray) -> Optional[SpinSystem]:
        groups = [self.codec.class_group(c) for c in classes]
        g = len(groups)
        full = np.zeros((self.g_max, self.g_max))
        full[self.pair_rows.cpu().numpy(), self.pair_cols.cpu().numpy()] = j_upper
        full = full + full.T
        try:
            return SpinSystem.from_group_couplings([n for n, _ in groups], [s for _, s in groups], full[:g, :g])
        except ValueError:
            return None

    @torch.no_grad()
    def propose(self, features, frequency_hz, k: int = 5, j_samples: int = 2,
                generator: torch.Generator = None, **_ignored) -> List[List[Interpretation]]:
        was_training = self.training
        self.eval()
        out = self.forward({"features": features, "frequency_hz": frequency_hz})
        if was_training:
            self.train()
        scale = self.model_config.set_head.j_scale_hz
        results = []
        for b in range(features.shape[0]):
            p_exist = torch.sigmoid(out["exists"][b]).cpu().numpy()
            logp = out["group"][b].log_softmax(-1).cpu().numpy()
            decoded = [self.constrained_groups(logp[i]) for i in range(self.k_slots)]
            group_score = [float(sum(logp[i][p, c] for p, c in enumerate(d))) if d else -np.inf
                           for i, d in enumerate(decoded)]
            mean = out["j_mean"][b].cpu().numpy() * scale
            std = np.exp(0.5 * out["j_logvar"][b].cpu().numpy()) * scale
            contrib = out["contribution"][b].cpu().numpy()
            order = np.argsort(-p_exist)
            subsets = [[i for i in order if p_exist[i] > 0.5] or [int(order[0])]]
            uncertain = sorted(range(self.k_slots), key=lambda i: abs(p_exist[i] - 0.5))
            for i in uncertain[:k]:
                alt = sorted(set(subsets[0]) ^ {i})
                if alt and alt not in subsets:
                    subsets.append(alt)
            rng = np.random.default_rng(0 if generator is None else int(generator.initial_seed()))
            candidates, seen = [], set()
            for subset in subsets:
                for draw in range(1 + j_samples):
                    comps, score = [], 0.0
                    for i in subset:
                        if decoded[i] is None:
                            continue
                        j = mean[i] if draw == 0 else mean[i] + rng.normal(size=mean[i].shape) * std[i]
                        system = self._slot_system(decoded[i], j)
                        if system is None:
                            continue
                        comps.append(Component(system, float(10 ** np.clip(contrib[i], -12, 0)), f"slot{i}"))
                        score += math.log(max(p_exist[i], 1e-9)) + group_score[i]
                    for i in set(range(self.k_slots)) - set(subset):
                        score += math.log(max(1 - p_exist[i], 1e-9))
                    if not comps:
                        continue
                    interp = Interpretation(tuple(comps), score - 0.5 * draw, {"source": "cnn_set"})
                    key = interp.key(self.spec.nuclei)
                    if key not in seen:
                        seen.add(key)
                        candidates.append(interp)
            candidates.sort(key=lambda c: -c.score)
            results.append(candidates[:k])
        return results
