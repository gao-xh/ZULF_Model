"""Target representations: token sequences (Transformer) and padded sets (CNN baseline).

Token grammar for one interpretation (components sorted by contribution):

    BOS ( COMP W<k> G<nucleus>x<size>+ SEP J<k>{G(G-1)/2} )+ EOS

Groups are magnetic-equivalence groups in canonical order (nucleus order of
the spec, larger groups first, then refined coupling colors). J tokens list the
upper triangle of the group coupling matrix row by row. The global J sign is
canonicalized first. All sizes come from `ProblemSpec`.

Canonical orders are training conveniences; evaluation always uses
permutation matching (`spinsystem.best_permutation`).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .spec import ProblemSpec
from .spinsystem import Component, Interpretation, SpinSystem, canonical_group_order, canonical_sign

PAD, BOS, EOS, COMP, SEP = "<pad>", "<bos>", "<eos>", "<comp>", "<sep>"
SPECIALS = (PAD, BOS, EOS, COMP, SEP)


class JBinner:
    def __init__(self, spec: ProblemSpec):
        self.edges = spec.j_bins.edges()
        self.centers = (self.edges[:-1] + self.edges[1:]) / 2
        self.widths = np.diff(self.edges)

    def __len__(self) -> int:
        return len(self.centers)

    def encode(self, value: float) -> Tuple[int, float]:
        """Bin index and offset in [-0.5, 0.5] bin widths (values beyond range are clipped)."""
        v = float(np.clip(value, self.edges[0], self.edges[-1]))
        index = int(np.clip(np.searchsorted(self.edges, v, side="right") - 1, 0, len(self) - 1))
        return index, float((v - self.centers[index]) / self.widths[index])

    def decode(self, index: int, offset: float = 0.0) -> float:
        return float(self.centers[index] + np.clip(offset, -0.5, 0.5) * self.widths[index])


class ContributionBinner:
    def __init__(self, spec: ProblemSpec):
        c = spec.contribution
        self.edges = np.linspace(c.log10_min, c.log10_max, c.bins + 1)
        self.centers = (self.edges[:-1] + self.edges[1:]) / 2

    def __len__(self) -> int:
        return len(self.centers)

    def encode(self, relative: float) -> int:
        value = math.log10(max(relative, 10 ** self.edges[0]))
        return int(np.clip(np.searchsorted(self.edges, value, side="right") - 1, 0, len(self) - 1))

    def decode(self, index: int) -> float:
        return float(10 ** self.centers[index])


class Vocabulary:
    def __init__(self, spec: ProblemSpec, j_bins: int, w_bins: int):
        tokens = list(SPECIALS)
        self.group_tokens: Dict[Tuple[str, int], int] = {}
        for nucleus in spec.nuclei:
            for size in range(1, spec.max_group_size + 1):
                self.group_tokens[(nucleus, size)] = len(tokens)
                tokens.append(f"G:{nucleus}x{size}")
        self.w_start = len(tokens)
        tokens += [f"W:{k}" for k in range(w_bins)]
        self.j_start = len(tokens)
        tokens += [f"J:{k}" for k in range(j_bins)]
        self.tokens = tokens
        self.index = {t: i for i, t in enumerate(tokens)}
        self.group_of = {v: k for k, v in self.group_tokens.items()}
        self.j_bins, self.w_bins = j_bins, w_bins

    def __len__(self) -> int:
        return len(self.tokens)

    def __getitem__(self, token: str) -> int:
        return self.index[token]

    def is_j(self, token_id: int) -> bool:
        return self.j_start <= token_id < self.j_start + self.j_bins

    def is_w(self, token_id: int) -> bool:
        return self.w_start <= token_id < self.w_start + self.w_bins

    def is_group(self, token_id: int) -> bool:
        return token_id in self.group_of


@dataclass
class GrammarState:
    phase: str = "start"      # start, comp, weight, groups, j, between, done
    components: int = 0
    groups: List[Tuple[str, int]] = field(default_factory=list)
    j_remaining: int = 0

    def copy(self) -> "GrammarState":
        return GrammarState(self.phase, self.components, list(self.groups), self.j_remaining)


class Grammar:
    """Finite-state constraints used by decoders and by sequence validation."""

    def __init__(self, spec: ProblemSpec, vocab: Vocabulary):
        self.spec, self.vocab = spec, vocab
        self.order = {n: i for i, n in enumerate(spec.nuclei)}

    def _spins(self, groups) -> int:
        return sum(size for _, size in groups)

    def _group_allowed(self, state: GrammarState, key: Tuple[str, int]) -> bool:
        if self._spins(state.groups) + key[1] > self.spec.max_spins:
            return False
        if state.groups:
            last = state.groups[-1]
            if (self.order[key[0]], -key[1]) < (self.order[last[0]], -last[1]):
                return False
        return True

    def _can_close(self, state: GrammarState) -> bool:
        nuclei = {n for n, _ in state.groups}
        return self._spins(state.groups) in self.spec.spin_counts and len(nuclei) >= 2

    def allowed(self, state: GrammarState) -> np.ndarray:
        v = self.vocab
        mask = np.zeros(len(v), dtype=bool)
        if state.phase == "start":
            mask[v[BOS]] = True
        elif state.phase == "comp":
            mask[v[COMP]] = True
        elif state.phase == "weight":
            mask[v.w_start:v.w_start + v.w_bins] = True
        elif state.phase == "groups":
            for key, tid in v.group_tokens.items():
                mask[tid] = self._group_allowed(state, key)
            mask[v[SEP]] = self._can_close(state)
        elif state.phase == "j":
            mask[v.j_start:v.j_start + v.j_bins] = True
        elif state.phase == "between":
            mask[v[EOS]] = True
            mask[v[COMP]] = state.components < self.spec.max_components
        elif state.phase == "done":
            mask[v[PAD]] = True
        return mask

    def advance(self, state: GrammarState, token_id: int) -> GrammarState:
        if not self.allowed(state)[token_id]:
            raise ValueError(f"Token {self.vocab.tokens[token_id]} not allowed in phase {state.phase}.")
        s = state.copy()
        v = self.vocab
        if s.phase == "start":
            s.phase = "comp"
        elif token_id == v[COMP]:
            s.phase, s.groups = "weight", []
            s.components += 1
        elif s.phase == "weight":
            s.phase = "groups"
        elif s.phase == "groups" and token_id == v[SEP]:
            g = len(s.groups)
            s.j_remaining = g * (g - 1) // 2
            s.phase = "j" if s.j_remaining else "between"
        elif s.phase == "groups":
            s.groups.append(v.group_of[token_id])
        elif s.phase == "j":
            s.j_remaining -= 1
            if s.j_remaining == 0:
                s.phase = "between"
        elif token_id == v[EOS]:
            s.phase = "done"
        return s

    def validate(self, tokens: Sequence[int]) -> bool:
        state = GrammarState()
        try:
            for t in tokens:
                state = self.advance(state, int(t))
        except ValueError:
            return False
        return state.phase == "done"

    def max_length(self) -> int:
        g = self.spec.max_spins
        per_component = 3 + g + g * (g - 1) // 2
        return 2 + self.spec.max_components * per_component


@dataclass
class EncodedSequence:
    tokens: List[int]
    j_offsets: List[float]   # offset in bin widths at J positions, NaN elsewhere


class InterpretationCodec:
    def __init__(self, spec: ProblemSpec):
        self.spec = spec
        self.j = JBinner(spec)
        self.w = ContributionBinner(spec)
        self.vocab = Vocabulary(spec, len(self.j), len(self.w))
        self.grammar = Grammar(spec, self.vocab)

    # -- component canonical form -------------------------------------------
    def canonical_groups(self, system: SpinSystem) -> Tuple[List[Tuple[str, int]], np.ndarray]:
        """Canonically ordered (nucleus, size) groups and group coupling matrix (sign fixed)."""
        sign = canonical_sign(system)
        groups = canonical_group_order(system, self.spec.nuclei)
        reps = [g[0] for g in groups]
        gj = sign * system.couplings_hz[np.ix_(reps, reps)]
        np.fill_diagonal(gj, 0.0)
        return [(system.isotopes[g[0]], len(g)) for g in groups], gj

    def _ordered_components(self, interpretation: Interpretation) -> List[Component]:
        return sorted(interpretation.components, key=lambda c: -c.contribution)

    # -- sequences -------------------------------------------------------------
    def encode(self, interpretation: Interpretation) -> EncodedSequence:
        v = self.vocab
        comps = self._ordered_components(interpretation)
        if len(comps) > self.spec.max_components:
            raise ValueError("Too many components for the problem spec.")
        top = max(c.contribution for c in comps) or 1.0
        tokens, offsets = [v[BOS]], [math.nan]
        for c in comps:
            groups, gj = self.canonical_groups(c.system)
            if any(size > self.spec.max_group_size for _, size in groups):
                raise ValueError("Group larger than max_group_size.")
            tokens += [v[COMP], v.w_start + self.w.encode(c.contribution / top)]
            offsets += [math.nan, math.nan]
            for key in groups:
                tokens.append(v.group_tokens[key])
                offsets.append(math.nan)
            tokens.append(v[SEP])
            offsets.append(math.nan)
            for a in range(len(groups)):
                for b in range(a + 1, len(groups)):
                    index, offset = self.j.encode(gj[a, b])
                    tokens.append(v.j_start + index)
                    offsets.append(offset)
        tokens.append(v[EOS])
        offsets.append(math.nan)
        if not self.grammar.validate(tokens):
            raise ValueError("Encoded sequence violates the grammar (check spec limits).")
        return EncodedSequence(tokens, offsets)

    def decode(self, tokens: Sequence[int], j_offsets: Optional[Sequence[float]] = None,
               score: Optional[float] = None) -> Interpretation:
        v = self.vocab
        if not self.grammar.validate(list(tokens)):
            raise ValueError("Token sequence violates the grammar.")
        offsets = list(j_offsets) if j_offsets is not None else [0.0] * len(tokens)
        components: List[Component] = []
        i = 1
        while tokens[i] == v[COMP]:
            weight = self.w.decode(tokens[i + 1] - v.w_start)
            i += 2
            groups = []
            while tokens[i] != v[SEP]:
                groups.append(v.group_of[tokens[i]])
                i += 1
            i += 1
            g = len(groups)
            gj = np.zeros((g, g))
            for a in range(g):
                for b in range(a + 1, g):
                    off = offsets[i] if offsets[i] is not None and not math.isnan(offsets[i]) else 0.0
                    gj[a, b] = gj[b, a] = self.j.decode(tokens[i] - v.j_start, off)
                    i += 1
            system = SpinSystem.from_group_couplings([n for n, _ in groups], [s for _, s in groups], gj)
            components.append(Component(system, weight, "decoded"))
        return Interpretation(tuple(components), score)

    # -- padded set targets ----------------------------------------------------------
    def encode_set(self, interpretation: Interpretation) -> Dict[str, np.ndarray]:
        """Fixed-size targets with masks: K = max_components, S = max_spins.

        Spins are expanded from canonical groups; `group_id` marks equivalence.
        Isotope class index len(nuclei) is padding.
        """
        k_max, s_max = self.spec.max_components, self.spec.max_spins
        pad = len(self.spec.nuclei)
        out = {"component_mask": np.zeros(k_max, np.float32),
               "isotope": np.full((k_max, s_max), pad, np.int64),
               "spin_mask": np.zeros((k_max, s_max), np.float32),
               "group_id": np.full((k_max, s_max), -1, np.int64),
               "couplings": np.zeros((k_max, s_max, s_max), np.float32),
               "log10_contribution": np.zeros(k_max, np.float32)}
        comps = self._ordered_components(interpretation)
        top = max(c.contribution for c in comps) or 1.0
        for k, c in enumerate(comps[:k_max]):
            groups, gj = self.canonical_groups(c.system)
            sizes = [s for _, s in groups]
            index = np.repeat(np.arange(len(groups)), sizes)
            n = len(index)
            j = gj[np.ix_(index, index)]
            j[index[:, None] == index[None, :]] = 0.0
            out["component_mask"][k] = 1
            out["isotope"][k, :n] = [self.spec.nucleus_index(groups[g][0]) for g in index]
            out["spin_mask"][k, :n] = 1
            out["group_id"][k, :n] = index
            out["couplings"][k, :n, :n] = j
            out["log10_contribution"][k] = math.log10(max(c.contribution / top, 1e-12))
        return out

    def decode_set(self, arrays: Dict[str, np.ndarray], threshold: float = 0.5,
                   score: Optional[float] = None) -> Optional[Interpretation]:
        """Decode set predictions (probabilities or hard values) into an interpretation."""
        components = []
        for k in range(self.spec.max_components):
            if arrays["component_mask"][k] < threshold:
                continue
            iso_idx = np.asarray(arrays["isotope"][k])
            keep = (iso_idx < len(self.spec.nuclei)) & (np.asarray(arrays["spin_mask"][k]) >= threshold)
            if keep.sum() < 2:
                continue
            isotopes = tuple(self.spec.nuclei[int(i)] for i in iso_idx[keep])
            j = np.asarray(arrays["couplings"][k], float)[np.ix_(keep, keep)]
            j = (j + j.T) / 2
            np.fill_diagonal(j, 0.0)
            system = SpinSystem(isotopes, j)
            contribution = float(10 ** arrays["log10_contribution"][k])
            components.append(Component(system, contribution, "decoded"))
        return Interpretation(tuple(components), score) if components else None
