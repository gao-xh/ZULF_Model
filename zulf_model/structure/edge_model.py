"""Route B of J -> structure: a bond-count classifier learned from generator pairs, used as the likelihood of
zulf_hypothesis.j_structure.rank_structures.

Each observed coupling of a generator observation (`observation_from_graph`) is one training row: the coupling
value J, its kind (13C to proton, or proton to proton), the proton counts of the carbons at the two ends, and the
label = its bond count in the true graph (1..MAX_BONDS, longer paths share the last class). A small MLP learns
p(bonds | J, kind, h). Ranking needs p(J | bonds, kind, h) up to a factor that does not depend on bonds; by Bayes it
is p(bonds | J, kind, h) / p(bonds | kind, h), where the denominator is the class frequency of the training rows
with that kind and h (`EdgeModel.prior`). `LearnedLikelihood` wraps the pair with a floor (a value the training
never saw keeps a finite log) and a cache, so rank_structures calls it once per coupling and bond count. A second
MLP learns p(hybridization | 1J, h) from the 1J rows, used the same way for the 1J terms (`LearnedLikelihood.hybrid`;
rank_structures picks it up instead of its rule `one_bond_likelihood`). zulf_model does not import zulf_hypothesis;
the observation format is shared by convention (JObservation.from_dict).

Input mode (`EdgeModel.mode`): "J" (the coupling in Hz, default) or "K" (the reduced coupling, D51, expressed as
the equivalent 13C-1H coupling J (g_13C g_1H) / (g_A g_B) so that 13C-1H rows are identical in both modes). A kind
names the elements of a pair (CH: carbon-hydrogen, HH: hydrogen-hydrogen); the isotopes default to 13C and 1H and
can be given per row (`nuclei`), e.g. 2H for deuterated samples, where only the K mode keeps the trained scale.

The classifier only knows the coupling rules it was trained on (configs/couplings_v1.json): on synthetic
observations from the same rules it has an advantage over route A; the real J networks are the fair comparison.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from ..generator.couplings import CouplingRules
from ..generator.graphs import GraphConfig, random_graph
from .observations import observation_from_graph

MAX_BONDS = 6
KIND_NUCLEI = {"CH": ("13C", "1H"), "HH": ("1H", "1H")}      # default isotopes of each element pair
K_REFERENCE = ("13C", "1H")
N_H = 4                          # proton counts 0..3 one-hot
HYBRIDS = ("sp3", "sp2", "sp")


def model_value(j: float, kind: str, mode: str = "J", nuclei: Optional[Tuple[str, str]] = None) -> float:
    """The coupling as the model sees it: J (Hz), or in mode "K" the equivalent 13C-1H coupling at equal K."""
    if mode == "J":
        return float(j)
    if mode != "K":
        raise ValueError("mode must be 'J' or 'K'")
    from zulf_core.nuclei import convert_coupling
    return convert_coupling(j, tuple(nuclei or KIND_NUCLEI[kind]), K_REFERENCE)


def features(j: float, kind: str, h: Sequence[int], mode: str = "J",
             nuclei: Optional[Tuple[str, str]] = None) -> np.ndarray:
    """Input row: J scaled, |J| on a log scale, kind, one-hot proton counts of the two ends (J replaced by the
    equivalent 13C-1H coupling in mode "K")."""
    j = model_value(j, kind, mode, nuclei)
    out = np.zeros(4 + 2 * N_H, np.float32)
    out[0] = j / 20.0
    out[1] = np.log1p(abs(j)) / 5.0
    out[2] = np.sign(j) * min(abs(j), 20.0) / 20.0
    out[3] = 1.0 if kind == "CH" else 0.0
    out[4 + min(int(h[0]), N_H - 1)] = 1.0
    out[4 + N_H + min(int(h[1]), N_H - 1)] = 1.0
    return out


def hybrid_features(j: float, h: int, mode: str = "J", nuclei: Optional[Tuple[str, str]] = None) -> np.ndarray:
    j = model_value(j, "CH", mode, nuclei)
    out = np.zeros(2 + N_H, np.float32)
    out[0] = (j - 150.0) / 30.0
    out[1] = ((j - 150.0) / 30.0) ** 2
    out[2 + min(int(h), N_H - 1)] = 1.0
    return out


def hybrid_rows(obs: dict) -> List[Tuple[float, int, int]]:
    """(1J, protons, hybridization class) for every unit whose 1J is observed."""
    rows = []
    for key, j in obs["couplings"].items():
        a, b = (t.strip() for t in key[2:-1].split(","))
        c, g = (a, b) if b in obs["protons"] else (b, a)
        if c in obs["units"] and obs["protons"].get(g, [None, None]) == [c, 0]:
            rows.append((float(j), obs["units"][c]["h"], HYBRIDS.index(obs["hybridization"][c])))
    return rows


def rows_from_observation(obs: dict) -> List[Tuple[float, str, Tuple[int, int], int]]:
    """(J, kind, h, bonds) for every coupling of an observation dict that carries bond_counts."""
    units, protons = obs["units"], obs["protons"]
    rows = []
    for key, j in obs["couplings"].items():
        a, b = (t.strip() for t in key[2:-1].split(","))
        if a in protons and b in protons:
            kind, h = "HH", tuple(sorted((units[protons[a][0]]["h"], units[protons[b][0]]["h"])))
        else:
            c, g = (a, b) if b in protons else (b, a)
            kind, h = "CH", (units[c]["h"], units[protons[g][0]]["h"])
        rows.append((float(j), kind, h, min(int(obs["bond_counts"][key]), MAX_BONDS)))
    return rows


def generate_rows(n_graphs: int, seed: int = 0, config: Optional[GraphConfig] = None,
                  rules: Optional[CouplingRules] = None, sigma: float = 0.1, max_unseen: int = 2,
                  with_hybrid: bool = False):
    """Training rows from n_graphs accepted generator graphs (and the 1J hybridization rows when with_hybrid)."""
    rng = np.random.default_rng(seed)
    config = config or GraphConfig(heavy_atoms=(2, 6))
    rules = rules or CouplingRules()
    rows, hyb, accepted, tries = [], [], 0, 0
    while accepted < n_graphs and tries < 50 * n_graphs:
        tries += 1
        out = observation_from_graph(rng, random_graph(rng, config), rules, sigma=sigma, max_unseen=max_unseen)
        if out is None:
            continue
        rows.extend(rows_from_observation(out[0]))
        hyb.extend(hybrid_rows(out[0]))
        accepted += 1
    return (rows, hyb) if with_hybrid else rows


@dataclass
class EdgeModel:
    """The trained classifier (weights as numpy arrays, so inference needs no torch) and the class priors."""
    weights: List[Tuple[np.ndarray, np.ndarray]]
    prior: Dict[Tuple[str, Tuple[int, int]], np.ndarray]
    prior_all: Dict[str, np.ndarray]
    hybrid_weights: Optional[List[Tuple[np.ndarray, np.ndarray]]] = None
    hybrid_prior: Optional[Dict[int, np.ndarray]] = None
    mode: str = "J"

    def predict(self, x: np.ndarray) -> np.ndarray:
        """Class probabilities, columns bonds 1..MAX_BONDS."""
        return _forward(self.weights, x)

    def predict_hybrid(self, x: np.ndarray) -> np.ndarray:
        """Hybridization probabilities, columns HYBRIDS."""
        return _forward(self.hybrid_weights, x)

    def class_prior(self, kind: str, h: Tuple[int, int]) -> np.ndarray:
        return self.prior.get((kind, tuple(h)), self.prior_all[kind])

    def to_dict(self) -> dict:
        out = {"weights": [[w.tolist(), b.tolist()] for w, b in self.weights],
               "prior": [[k[0], list(k[1]), v.tolist()] for k, v in self.prior.items()],
               "prior_all": {k: v.tolist() for k, v in self.prior_all.items()}, "mode": self.mode}
        if self.hybrid_weights is not None:
            out["hybrid_weights"] = [[w.tolist(), b.tolist()] for w, b in self.hybrid_weights]
            out["hybrid_prior"] = {str(k): v.tolist() for k, v in self.hybrid_prior.items()}
        return out

    @classmethod
    def from_dict(cls, data: dict) -> "EdgeModel":
        hw = data.get("hybrid_weights")
        return cls([(np.array(w), np.array(b)) for w, b in data["weights"]],
                   {(k, tuple(h)): np.array(v) for k, h, v in data["prior"]},
                   {k: np.array(v) for k, v in data["prior_all"].items()},
                   None if hw is None else [(np.array(w), np.array(b)) for w, b in hw],
                   None if hw is None else {int(k): np.array(v) for k, v in data["hybrid_prior"].items()},
                   data.get("mode", "J"))


def _forward(weights, x: np.ndarray) -> np.ndarray:
    z = np.atleast_2d(x).astype(np.float64)
    for k, (w, b) in enumerate(weights):
        z = z @ w + b
        if k < len(weights) - 1:
            z = np.maximum(z, 0.0)
    z -= z.max(axis=1, keepdims=True)
    p = np.exp(z)
    return p / p.sum(axis=1, keepdims=True)


def _priors(rows, smoothing: float = 1.0):
    counts: Dict[Tuple[str, Tuple[int, int]], np.ndarray] = {}
    totals = {"CH": np.full(MAX_BONDS, smoothing), "HH": np.full(MAX_BONDS, smoothing)}
    totals["HH"][0] = 0.0
    for _, kind, h, bonds in rows:
        c = counts.setdefault((kind, tuple(h)), np.full(MAX_BONDS, smoothing) * (totals[kind] > 0))
        c[bonds - 1] += 1
        totals[kind][bonds - 1] += 1
    norm = lambda v: v / v.sum()
    return {k: norm(v) for k, v in counts.items()}, {k: norm(v) for k, v in totals.items()}


def _train_mlp(x: np.ndarray, y: np.ndarray, n_classes: int, hidden, epochs, lr, batch, seed, verbose, threads):
    import torch

    torch.set_num_threads(threads)
    torch.manual_seed(seed)
    x = torch.tensor(np.asarray(x, np.float32))
    y = torch.tensor(np.asarray(y, np.int64))
    sizes = [x.shape[1], *hidden, n_classes]
    layers = []
    for k in range(len(sizes) - 1):
        layers.append(torch.nn.Linear(sizes[k], sizes[k + 1]))
        if k < len(sizes) - 2:
            layers.append(torch.nn.ReLU())
    net = torch.nn.Sequential(*layers)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    gen = torch.Generator().manual_seed(seed)
    for epoch in range(epochs):
        order = torch.randperm(len(y), generator=gen)
        total = 0.0
        for s in range(0, len(y), batch):
            idx = order[s:s + batch]
            loss = torch.nn.functional.cross_entropy(net(x[idx]), y[idx])
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss.detach()) * len(idx)
        if verbose:
            print(f"  epoch {epoch + 1:3d}  loss {total / len(y):.4f}")
    return [(m.weight.detach().numpy().T.astype(np.float64), m.bias.detach().numpy().astype(np.float64))
            for m in net if isinstance(m, torch.nn.Linear)]


def train_edge_model(rows, hidden: Sequence[int] = (64, 64), epochs: int = 30, lr: float = 3e-3,
                     batch: int = 512, seed: int = 0, verbose: bool = False, threads: int = 1,
                     hybrid_rows=None, mode: str = "J") -> EdgeModel:
    """Cross-entropy training of the bond-count MLP (torch, CPU; one thread is fastest for a net this small), and
    of the hybridization MLP when hybrid_rows [(1J, protons, class)] are given."""
    x = np.stack([features(j, kind, h, mode) for j, kind, h, _ in rows])
    y = np.array([bonds - 1 for *_, bonds in rows])
    weights = _train_mlp(x, y, MAX_BONDS, hidden, epochs, lr, batch, seed, verbose, threads)
    prior, prior_all = _priors(rows)
    hw = hp = None
    if hybrid_rows:
        xh = np.stack([hybrid_features(j, h, mode) for j, h, _ in hybrid_rows])
        yh = np.array([c for *_, c in hybrid_rows])
        hw = _train_mlp(xh, yh, len(HYBRIDS), (16,), epochs, lr, batch, seed, verbose, threads)
        hp = {}
        for _, h, c in hybrid_rows:
            hp.setdefault(int(h), np.ones(len(HYBRIDS)))[c] += 1
        hp = {h: v / v.sum() for h, v in hp.items()}
    return EdgeModel(weights, prior, prior_all, hw, hp, mode)


class LearnedLikelihood:
    """likelihood(j, kind, bonds, sigma, h) for rank_structures: p(bonds | J, kind, h) / p(bonds | kind, h), with a
    floor. sigma is not used (the training noise sets the width)."""

    def __init__(self, model: EdgeModel, floor: float = 1e-3):
        self.model, self.floor = model, floor
        self._cache: Dict[tuple, np.ndarray] = {}
        if model.hybrid_weights is not None:      # without the head, rank_structures uses its rule 1J likelihood
            self.hybrid = self._hybrid

    def ratios(self, j: float, kind: str, h: Tuple[int, int]) -> np.ndarray:
        key = (round(j, 6), kind, tuple(h))
        if key not in self._cache:
            p = self.model.predict(features(j, kind, h, self.model.mode))[0]
            prior = self.model.class_prior(kind, tuple(h))
            self._cache[key] = (p + self.floor) / (prior + self.floor)
        return self._cache[key]

    def __call__(self, j: float, kind: str, bonds: int, sigma: float = 0.1, h=(0, 0)) -> float:
        return float(self.ratios(j, kind, h)[min(bonds, MAX_BONDS) - 1])

    def _hybrid(self, j: float, hyb: str, sigma: float = 0.1, h=(0, 0)) -> float:
        """p(hyb | 1J, h) / p(hyb | h)."""
        p = self.model.predict_hybrid(hybrid_features(j, h[0], self.model.mode))[0]
        prior = self.model.hybrid_prior.get(int(h[0]), np.full(len(HYBRIDS), 1.0 / len(HYBRIDS)))
        k = HYBRIDS.index(hyb)
        return float((p[k] + self.floor) / (prior[k] + self.floor))


def accuracy(model: EdgeModel, rows: Iterable, nuclei: Optional[Dict[str, Tuple[str, str]]] = None) -> float:
    """Bond-count accuracy; `nuclei` {kind: (isotope, isotope)} when the rows are not 13C / 1H (e.g. 2H)."""
    rows = list(rows)
    x = np.stack([features(j, kind, h, model.mode, (nuclei or {}).get(kind)) for j, kind, h, _ in rows])
    pred = model.predict(x).argmax(axis=1) + 1
    return float(np.mean(pred == np.array([b for *_, b in rows])))
