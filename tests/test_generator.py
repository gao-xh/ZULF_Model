import tempfile
import unittest

import numpy as np

from zulf_model.generator import (CouplingRules, GraphConfig, IsotopologueConfig, MoleculeGraph, SplitConfig,
                                  build_default_sampler, iter_samples, natural_isotopologues, random_graph,
                                  random_system, split_of, write_shards)
from zulf_model.generator.graphs import VALENCE
from zulf_model.generator.random_j import RandomSystemConfig
from zulf_model.physics import compute_transitions
from zulf_model.spec import ProblemSpec


def isopropylamine() -> MoleculeGraph:
    # Atoms: 0 methine C, 1 N, 2 methyl C, 3 methyl C, then hydrogens.
    elements = ["C", "N", "C", "C"]
    bonds = {(0, 1): 1, (0, 2): 1, (0, 3): 1}
    graph = MoleculeGraph(elements, bonds, [False] * 4)
    for heavy, count in ((0, 1), (1, 2), (2, 3), (3, 3)):
        for _ in range(count):
            graph.elements.append("H")
            graph.aromatic.append(False)
            graph.add_bond(heavy, len(graph.elements) - 1)
    graph.validate()
    return graph


class GraphTests(unittest.TestCase):
    def test_random_graphs_satisfy_valence_and_connectivity(self):
        rng = np.random.default_rng(0)
        config = GraphConfig(heavy_atoms=(2, 9), aromatic_ring_probability=0.5, ring_closure_probability=0.5)
        for _ in range(200):
            g = random_graph(rng, config)
            for i, e in enumerate(g.elements):
                self.assertEqual(g.used_valence(i), VALENCE[e])
            self.assertTrue((g.distances() < 10 ** 6).all())

    def test_skeleton_hash_ignores_labels_and_order(self):
        g = isopropylamine()
        h = g.copy()
        h.labels = {0: "13C"}
        self.assertEqual(g.skeleton_hash(), h.skeleton_hash())
        self.assertEqual(MoleculeGraph.from_dict(g.to_dict()).skeleton_hash(), g.skeleton_hash())

    def test_hybridization(self):
        g = MoleculeGraph(["C", "C"], {(0, 1): 2}, [False, False])
        self.assertEqual(g.hybridization(0), "sp2")


class IsotopologueTests(unittest.TestCase):
    def setUp(self):
        self.spec = ProblemSpec(nuclei=("1H", "13C", "15N"), spin_counts=(8,), max_components=4)

    def test_isopropylamine_natural_isotopologues(self):
        rng = np.random.default_rng(1)
        comps = natural_isotopologues(rng, isopropylamine(), CouplingRules(), self.spec,
                                      IsotopologueConfig(include_natural_n15=False))
        by_label = {c.label: c for c in comps}
        self.assertEqual(set(by_label), {"13C@0", "13C@2"})
        methine = by_label["13C@0"].system
        methyl = by_label["13C@2"].system
        self.assertEqual(methine.group_signature(), (("13C", 1), ("1H", 1), ("1H", 6)))
        self.assertEqual(methyl.group_signature(), (("13C", 1), ("1H", 1), ("1H", 3), ("1H", 3)))
        # Two equivalent methyl carbons: twice the methine contribution.
        self.assertAlmostEqual(by_label["13C@2"].contribution / by_label["13C@0"].contribution, 2.0)
        # Shared H-H couplings across isotopologues of one molecule.
        np.testing.assert_allclose(methine.couplings_hz[:7, :7], methyl.couplings_hz[:7, :7])
        self.assertGreater(len(compute_transitions(methine)), 0)

    def test_n15_isotopologue_has_negative_one_bond_coupling_when_nh_kept(self):
        rng = np.random.default_rng(2)
        spec = ProblemSpec(spin_counts=(10,))
        comps = natural_isotopologues(rng, isopropylamine(), CouplingRules(), spec,
                                      IsotopologueConfig(exchangeable_protons="keep"))
        n15 = [c for c in comps if c.label.startswith("15N")]
        self.assertEqual(len(n15), 1)
        s = n15[0].system
        n_index = s.isotopes.index("15N")
        self.assertLess(s.couplings_hz[n_index].min(), -60)

    def test_random_system_respects_spec(self):
        rng = np.random.default_rng(3)
        spec = ProblemSpec(spin_counts=(5, 7), max_group_size=3)
        for _ in range(50):
            s = random_system(rng, spec, RandomSystemConfig())
            self.assertIn(s.n_spins, (5, 7))
            self.assertLessEqual(max(len(g) for g in s.groups), 3)
            self.assertGreater(len(compute_transitions(s)), 0)


class SamplerTests(unittest.TestCase):
    def test_reproducible_and_split_consistent(self):
        spec = ProblemSpec()
        a = [s.to_dict() for s in build_default_sampler(spec).generate(15, seed=7)]
        b = [s.to_dict() for s in build_default_sampler(spec).generate(15, seed=7)]
        self.assertEqual(a, b)
        sampler = build_default_sampler(spec, {"splits": {"val_fraction": 0.3, "test_unseen_fraction": 0.3}})
        for sample in sampler.generate(10, seed=3, split="val"):
            for family in sample.metadata.get("families", [sample.family_id]):
                self.assertEqual(split_of(family, "x", sampler.split_config), "val")
            for c in sample.interpretation.components:
                self.assertTrue(spec.allows(c.system.isotopes))

    def test_split_function(self):
        cfg = SplitConfig(val_fraction=0.2, test_unseen_fraction=0.2, test_seen_fraction=0.1)
        counts = {}
        for k in range(4000):
            s = split_of(f"fam{k}", f"sample{k}", cfg)
            counts[s] = counts.get(s, 0) + 1
        self.assertAlmostEqual(counts["test_unseen"] / 4000, 0.2, delta=0.03)
        self.assertAlmostEqual(counts["val"] / 4000, 0.2, delta=0.03)
        self.assertEqual(split_of("famX", "a", cfg), split_of("famX", "a", cfg))

    def test_storage_roundtrip(self):
        spec = ProblemSpec()
        samples = list(build_default_sampler(spec).generate(12, seed=5))
        with tempfile.TemporaryDirectory() as tmp:
            manifest = write_shards(samples, tmp, shard_size=5)
            self.assertEqual(manifest["count"], 12)
            self.assertEqual(len(manifest["shards"]), 3)
            again = list(iter_samples(tmp))
        self.assertEqual([s.to_dict() for s in again], [s.to_dict() for s in samples])


if __name__ == "__main__":
    unittest.main()
