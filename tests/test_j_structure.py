"""J -> structure (Phase 7): route A (rule likelihood ranking) and route B (learned bond-count likelihood)."""
import collections
import json
import math
import unittest
from pathlib import Path

import numpy as np

import itertools

from zulf_hypothesis.j_structure import (JObservation, _bond_orders, _distances, couplings_likelihood, describe,
                                         one_bond_likelihood, parse_bonds, rank_structures, same_structure)
from zulf_model.generator.couplings import CouplingRules
from zulf_model.generator.graphs import GraphConfig, random_graph
from zulf_model.structure.edge_model import (HYBRIDS, MAX_BONDS, EdgeModel, LearnedLikelihood, _priors, features,
                                             hybrid_rows, rows_from_observation, train_edge_model)
from zulf_model.structure.observations import observation_from_graph

ROOT = Path(__file__).resolve().parents[1]
RULES = CouplingRules.load(ROOT / "configs" / "couplings_v1.json")


def synthetic(n, seed, max_atoms=5):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        o = observation_from_graph(rng, random_graph(rng, GraphConfig(heavy_atoms=(2, 5))), RULES)
        if o is not None and len(JObservation.from_dict(o[0]).atoms()) <= max_atoms:
            out.append(o)
    return out


def truth_rank(obs, truth, **kw):
    ranked = rank_structures(obs, top=10, **kw)
    return next((k for k, c in enumerate(ranked) if same_structure(obs, c, truth)), None)


class RouteATests(unittest.TestCase):
    def test_likelihood_is_a_density(self):
        # numerical integral over J (independent of the closed-form erf convolution)
        grid = np.linspace(-60.0, 260.0, 64001)
        dj = grid[1] - grid[0]
        for kind, bonds in [("CH", 1), ("CH", 2), ("CH", 3), ("CH", 5), ("HH", 2), ("HH", 3), ("HH", 4), ("HH", 6)]:
            core = sum(couplings_likelihood(j, kind, bonds, 0.1) - 0.02 / 40.0 for j in grid) * dj
            self.assertAlmostEqual(core, 0.98, delta=0.005, msg=(kind, bonds))

    def test_bond_orders_match_brute_force(self):
        # independent: every product of orders 1..3 filtered by the valence budget and the multiple-bond count
        edges = [(0, 1), (1, 2), (2, 3), (0, 3), (1, 4)]
        free = [1, 2, 1, 2, 2]
        for max_multiple in (0, 1, 2, 3):
            got = set(_bond_orders(edges, free, 4, max_multiple))
            want = set()
            for orders in itertools.product((1, 2, 3), repeat=len(edges)):
                extra = [0] * len(free)
                for (a, b), o in zip(edges, orders):
                    extra[a] += o - 1
                    extra[b] += o - 1
                if all(e <= f for e, f in zip(extra, free)) and sum(o > 1 for o in orders) <= max_multiple:
                    want.add(orders)
            self.assertEqual(got, want, max_multiple)

    def test_one_bond_hybridization(self):
        self.assertGreater(one_bond_likelihood(126.0, "sp3"), 10 * one_bond_likelihood(126.0, "sp2"))
        self.assertGreater(one_bond_likelihood(158.0, "sp2"), 10 * one_bond_likelihood(158.0, "sp3"))
        self.assertGreater(one_bond_likelihood(250.0, "sp"), 10 * one_bond_likelihood(250.0, "sp2"))

    def test_double_bond_from_one_bond_couplings(self):
        # propene-like: CH3 (1J 127), CH (1J 155), CH2 (1J 157); the double bond goes between the sp2 carbons
        data = {"units": {"C1": {"h": 3}, "C2": {"h": 1}, "C3": {"h": 2}},
                "protons": {"HC1": ["C1", 0], "HC2": ["C2", 0], "HC3": ["C3", 0]},
                "couplings": {"J(C1,HC1)": 127.0, "J(C2,HC2)": 155.0, "J(C3,HC3)": 157.0, "J(C1,HC2)": -6.5,
                              "J(C2,HC1)": -6.8, "J(C2,HC3)": -1.0, "J(C3,HC2)": -1.2, "J(C1,HC3)": 6.0,
                              "J(C3,HC1)": 5.5, "J(HC1,HC2)": 6.5, "J(HC2,HC3)": 13.0, "J(HC1,HC3)": 1.4}}
        obs = JObservation.from_dict(data)
        best = rank_structures(obs, top=3)[0]
        self.assertTrue(same_structure(obs, best, "C1-C2, C2=C3"), describe(obs, best))

    def test_amine_networks_rank_truth_first(self):
        for path in sorted((ROOT / "configs" / "j_networks").glob("*.json")):
            data = json.loads(path.read_text())
            obs = JObservation.from_dict(data)
            self.assertEqual(truth_rank(obs, data["truth"]), 0, path.stem)

    def test_bond_counts_match_truth_graph(self):
        # bond counts recorded by the generator vs a BFS on the bond list it prints
        for obs_dict, truth in synthetic(15, seed=3):
            obs = JObservation.from_dict(obs_dict)
            nx, edges = parse_bonds(obs, truth)
            index = {a: i for i, a in enumerate(obs.atoms())}
            d = _distances(len(index) + nx, edges)
            for (a, b), _ in obs.couplings.items():
                ends = [index[obs.protons[g]] if g in obs.protons else index[(g, 0)] for g in (a, b)]
                offset = 2 if a in obs.protons and b in obs.protons else 1
                self.assertEqual(obs_dict["bond_counts"][f"J({a},{b})"], int(d[ends[0], ends[1]]) + offset)

    def test_synthetic_benchmark(self):
        ranks = [truth_rank(JObservation.from_dict(o), t) for o, t in synthetic(30, seed=11)]
        top1 = np.mean([r == 0 for r in ranks])
        top3 = np.mean([r is not None and r < 3 for r in ranks])
        self.assertGreaterEqual(top1, 0.8)
        self.assertGreaterEqual(top3, 0.95)


class RouteBTests(unittest.TestCase):
    def test_priors_are_class_frequencies(self):
        rows = [(130.0, "CH", (1, 1), 1)] * 6 + [(5.0, "CH", (1, 3), 3)] * 3 + [(7.0, "HH", (1, 3), 3)] * 4
        prior, prior_all = _priors(rows, smoothing=0.0)
        np.testing.assert_allclose(prior[("CH", (1, 1))], np.eye(MAX_BONDS)[0])
        np.testing.assert_allclose(prior_all["CH"], np.array([6, 0, 3, 0, 0, 0]) / 9.0)
        np.testing.assert_allclose(prior_all["HH"], np.eye(MAX_BONDS)[2])

    def test_rows_and_learned_likelihood(self):
        obs, _ = synthetic(1, seed=5)[0]
        rows = rows_from_observation(obs)
        self.assertEqual(len(rows), len(obs["couplings"]))
        rng = np.random.default_rng(0)
        # separable toy: 1-bond CH near 130 Hz, 3-bond CH near 5 Hz, 3-bond HH near 7 Hz
        rows = ([(rng.uniform(120, 140), "CH", (1, 3), 1) for _ in range(400)]
                + [(rng.uniform(2, 8), "CH", (1, 3), 3) for _ in range(400)]
                + [(rng.uniform(5, 9), "HH", (1, 3), 3) for _ in range(400)])
        model = train_edge_model(rows, hidden=(16,), epochs=40, seed=0)
        model = EdgeModel.from_dict(json.loads(json.dumps(model.to_dict())))
        p = model.predict(np.stack([features(131.0, "CH", (1, 3)), features(4.0, "CH", (1, 3))]))
        self.assertEqual(int(p[0].argmax()), 0)
        self.assertEqual(int(p[1].argmax()), 2)
        lik = LearnedLikelihood(model)
        prior = model.class_prior("CH", (1, 3))
        expect = (p[1, 2] + lik.floor) / (prior[2] + lik.floor)
        self.assertAlmostEqual(lik(4.0, "CH", 3, 0.1, (1, 3)), expect, places=6)
        self.assertEqual(lik(4.0, "CH", 9, 0.1, (1, 3)), lik(4.0, "CH", MAX_BONDS, 0.1, (1, 3)))
        self.assertGreater(lik(4.0, "CH", 3, 0.1, (1, 3)), lik(4.0, "CH", 1, 0.1, (1, 3)))
        self.assertFalse(hasattr(lik, "hybrid"))          # no head: rank_structures uses the rule 1J likelihood

    def test_learned_hybridization(self):
        obs, _ = synthetic(1, seed=5)[0]
        for j, h, c in hybrid_rows(obs):
            self.assertIn(c, range(len(HYBRIDS)))
        rng = np.random.default_rng(1)
        hyb = ([(rng.uniform(120, 140), 2, 0) for _ in range(300)] + [(rng.uniform(150, 170), 1, 1) for _ in range(300)]
               + [(rng.uniform(150, 170), 2, 1) for _ in range(100)])
        rows = [(130.0, "CH", (2, 2), 1)] * 50 + [(5.0, "CH", (2, 2), 3)] * 50
        model = train_edge_model(rows, hidden=(8,), epochs=60, seed=0, hybrid_rows=hyb)
        model = EdgeModel.from_dict(json.loads(json.dumps(model.to_dict())))
        lik = LearnedLikelihood(model)
        self.assertGreater(lik.hybrid(160.0, "sp2", 0.1, (2, 2)), lik.hybrid(160.0, "sp3", 0.1, (2, 2)))
        self.assertGreater(lik.hybrid(128.0, "sp3", 0.1, (2, 2)), lik.hybrid(128.0, "sp2", 0.1, (2, 2)))
        # prior of h = 2 from the rows (with the add-one smoothing of train_edge_model)
        np.testing.assert_allclose(model.hybrid_prior[2], np.array([301, 101, 1]) / 403.0)


if __name__ == "__main__":
    unittest.main()
