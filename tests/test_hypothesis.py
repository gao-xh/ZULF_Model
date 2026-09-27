import sys
import unittest
from pathlib import Path

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_hypothesis import (AddCoupledProton, Finding, Fragment, KnowledgeBase, Labeling, ProtonGroup,
                                  ReferenceEntry, Site, build_model, pair, propose_all, run_checks, template)

A13 = 0.0107      # 13C natural abundance in the registry
from zulf_core.solver import RefineSettings, refine
from zulf_core.spinsystem import SpinSystem

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_solver import observe  # noqa: E402

RANGES = [(100.0, 160.0), (230.0, 300.0)]


def grouped(isotopes, sizes, couplings):
    return SpinSystem.from_group_couplings(isotopes, sizes, np.asarray(couplings, float))


def same_lines(a, b, tol=1e-6):
    ta, tb = compute_transitions(a), compute_transitions(b)
    order_a, order_b = np.argsort(ta.frequencies_hz), np.argsort(tb.frequencies_hz)
    return (len(ta) == len(tb) and np.allclose(ta.frequencies_hz[order_a], tb.frequencies_hz[order_b], atol=tol)
            and np.allclose(ta.amplitudes[order_a], tb.amplitudes[order_b], atol=tol))


class BuilderTests(unittest.TestCase):
    def test_ch_ch3_matches_hand_built_isotopologues(self):
        f = template("CH-CH3", one_bond={"Ca": 147.0, "Cb": 129.0})
        m = build_model(f)
        self.assertEqual(m.component_labels, ["13C@Ca", "13C@Cb"])
        self.assertEqual(m.ratios, (1.0, 1.0))
        # Independent hand construction (the matrices used in the blind analyses).
        alpha = grouped(["13C", "1H", "1H"], [1, 1, 3], [[0, 147.0, -4.5], [147.0, 0, 7.0], [-4.5, 7.0, 0]])
        beta = grouped(["13C", "1H", "1H"], [1, 3, 1], [[0, 129.0, -4.5], [129.0, 0, 7.0], [-4.5, 7.0, 0]])
        self.assertTrue(same_lines(m.interpretation.components[0].system, alpha))
        self.assertTrue(same_lines(m.interpretation.components[1].system, beta))
        self.assertEqual(m.ties(), (("c0.J1-2", "c1.J1-2"),))
        self.assertEqual(m.ties(shared_rate=True)[-1], ("c0.log_rate0", "c1.log_rate0"))
        self.assertEqual(set(m.coupling_names), {"J(Ca,Ha)", "J(Ca,Hb)", "J(Cb,Ha)", "J(Cb,Hb)", "J(Ha,Hb)"})

    def test_isopropyl_merges_equivalent_methyls_only_where_symmetric(self):
        m = build_model(template("CH(CH3)2"))
        self.assertEqual(m.ratios, (1.0, 2.0))
        ch, methyl = (c.system for c in m.interpretation.components)
        self.assertEqual([len(g) for g in ch.groups], [1, 1, 6])
        self.assertEqual([len(g) for g in methyl.groups], [1, 1, 3, 3])
        # Magnetic equivalence: the merged 6-proton group equals two separate methyl groups with equal couplings.
        split = grouped(["13C", "1H", "1H", "1H"], [1, 1, 3, 3],
                        [[0, 146.0, -4.5, -4.5], [146.0, 0, 7.0, 7.0], [-4.5, 7.0, 0, 0.0], [-4.5, 7.0, 0.0, 0]])
        self.assertTrue(same_lines(ch, split))
        self.assertIn(("c0.J1-2", "c1.J1-2", "c1.J1-3"), m.ties())
        self.assertEqual(m.unspecified, ("J(Hb1,Hb2)",))
        self.assertEqual(m.fixed(), ("c1.J2-3",))

    def test_fifteen_n_ratio_and_omission_reasons(self):
        f = template("CH-CH3", substituent="N", heavy_neighbour="C")
        f = f.with_couplings({("X", "Ha"): -1.0, ("X", "Hb"): 0.5})
        m = build_model(f)
        self.assertEqual(m.component_labels, ["13C@Ca", "13C@Cb", "15N@X"])
        a_n = 0.00364     # exact label-set probabilities: the other sites unlabelled
        self.assertAlmostEqual(m.ratios[2], a_n * (1 - A13) / (A13 * (1 - a_n)), places=10)
        self.assertEqual([o[0] for o in m.omitted], ["13C@Cr"])
        # The 15N lines lie at a few Hz: outside the fitted ranges, so the component is omitted with a reason.
        ranged = build_model(f, ranges=RANGES)
        self.assertIn(("15N@X", "no line in the fitted ranges"), ranged.omitted)
        self.assertEqual(ranged.component_labels, ["13C@Ca", "13C@Cb"])

    def test_invalid_symmetry_is_rejected(self):
        sites = (Site("C1", "C"), Site("C2", "C"))
        protons = (ProtonGroup("H1", 3, "C1"), ProtonGroup("H2", 3, "C2"))
        swap = {"C1": "C2", "C2": "C1", "H1": "H2", "H2": "H1"}
        with self.assertRaises(ValueError):
            Fragment("bad", sites, protons, {pair("C1", "H1"): 125.0, pair("C2", "H2"): 126.0}, (swap,))
        with self.assertRaises(ValueError):
            Fragment("bad", sites, (ProtonGroup("H1", 3, "C1"), ProtonGroup("H2", 2, "C2")), {}, (swap,))

    def test_add_coupled_proton_and_fast_exchange(self):
        f0 = template("CH-CH3")
        self.assertEqual(propose_all(f0), [])            # not triggered without a finding
        names = [g.name for g in propose_all(f0, [Finding("rate_asymmetry", "info", "")])]
        self.assertEqual(names, ["CH-CH3 + HX on X", "CH-CH3 with CaH2", "CH-CH3 with CbH2"])
        f = AddCoupledProton(initial_hz=1.5).propose(f0)[0]
        m = build_model(f)
        self.assertEqual([[len(g) for g in c.system.groups] for c in m.interpretation.components],
                         [[1, 1, 3, 1], [1, 1, 3, 1]])
        self.assertEqual(len(m.ties()), 3)                # every proton-proton coupling is shared
        fast = build_model(f, include_exchangeable=False)
        self.assertEqual([[len(g) for g in c.system.groups] for c in fast.interpretation.components],
                         [[1, 1, 3], [1, 1, 3]])

    def test_change_proton_count(self):
        from zulf_hypothesis import ChangeProtonCount
        from zulf_hypothesis.motifs import MOTIFS
        move = ChangeProtonCount()
        ethyl = move.propose(template("CH-CH3"))[0]                       # CH -> CH2: the ethyl skeleton
        self.assertEqual([(p.label, p.size) for p in ethyl.protons], [("Ha", 2), ("Hb", 3)])
        self.assertEqual(ethyl.couplings, template("CH-CH3").couplings)   # couplings carried over
        iso = move.propose(template("CH(CH3)2"))
        self.assertEqual([[p.size for p in f.protons] for f in iso], [[2, 3, 3], [1, 2, 2]])   # orbit changes together
        self.assertEqual(build_model(iso[1]).component_labels[:1], ["13C@Ca"])
        self.assertEqual(move.propose(MOTIFS["benzene ring"].fragment({"A2": 158.0})), [])  # ring C-H left alone
        withx = AddCoupledProton().propose(template("CH-CH3"))[0]
        self.assertFalse(any("X" in f.name.split(" with ")[-1] for f in move.propose(withx)))  # exchangeable kept

    def test_fragment_round_trip(self):
        f = template("CH(CH3)2")
        g = Fragment.from_dict(f.to_dict())
        self.assertEqual(build_model(g).summary()["couplings"], build_model(f).summary()["couplings"])


class MinorIsotopologueTests(unittest.TestCase):
    def test_double_labels_follow_their_parent(self):
        f = template("CH-CH3", one_bond={"Ca": 146.0, "Cb": 130.0})
        m = build_model(f, max_labels=2)
        self.assertEqual(m.component_labels, ["13C@Ca", "13C@Cb", "13C@Ca+13C@Cb"])
        self.assertAlmostEqual(m.ratios[2], A13 / (1 - A13), places=10)   # both labelled vs one labelled
        self.assertEqual(m.parents, [None, None, 0])
        # Independent hand construction of the 13C-13C isotopologue.
        hand = grouped(["13C", "13C", "1H", "1H"], [1, 1, 1, 3],
                       [[0, 35.0, 146.0, -4.5], [35.0, 0, -4.5, 130.0], [146.0, -4.5, 0, 7.0], [-4.5, 130.0, 7.0, 0]])
        self.assertTrue(same_lines(m.interpretation.components[2].system, hand))
        # Free amplitudes: the minor set rides on its parent; no free amplitude or rate of its own.
        self.assertEqual(m.amplitude_map(False), ((1.0, 0.0), (0.0, 1.0), (A13 / (1 - A13), 0.0)))
        self.assertIn(("c0.log_rate0", "c2.log_rate0"), m.ties())
        # Gate: below min_ratio it is omitted with the reason.
        gated = build_model(f, max_labels=2, min_ratio=0.05)
        self.assertEqual(gated.component_labels, ["13C@Ca", "13C@Cb"])
        self.assertIn("below min_ratio", gated.omitted[0][1])

    def test_magnetically_inequivalent_pair_and_enrichment(self):
        m = build_model(template("CH(CH3)2"), max_labels=2)
        self.assertIn("13C@Cb1,Cb2 [not magnetically equivalent]", m.component_labels)
        pair_system = m.interpretation.components[m.component_labels.index(
            "13C@Cb1,Cb2 [not magnetically equivalent]")].system
        self.assertEqual([len(g) for g in pair_system.groups], [1, 1, 1, 3, 3])
        # Uniform enrichment at 99 %: the fully labelled molecule is the primary; partly labelled ones are minor.
        e = build_model(template("CH(CH3)2"), labeling=Labeling.enriched({"13C": 0.99}))
        full = [l for l in e.component_labels if l.count("13C") == 2 and "Ca" in l and "Cb1,Cb2" in l]
        self.assertEqual(len(full), 1)
        i = e.component_labels.index(full[0])
        self.assertIsNone(e.parents[i])
        self.assertEqual([p is None for p in e.parents].count(True), 1)
        self.assertTrue(all(e.ratios[j] < 0.05 * e.ratios[i] for j in range(len(e.ratios)) if j != i))
        # Site-specific enrichment ([Ca]-13C 99 %, natural elsewhere).
        s = build_model(template("CH-CH3"), labeling=Labeling.enriched(sites={"Ca": {"13C": 0.99}}), max_labels=2)
        self.assertEqual([p is None for p in s.parents], [True, False, False])
        # Site-specific labels must respect symmetry.
        with self.assertRaises(ValueError):
            build_model(template("CH(CH3)2"), labeling=Labeling.enriched(sites={"Cb1": {"13C": 0.99}}))

    def test_refinement_with_enriched_minor_isotopologue(self):
        truth = {"J(Ca,Ha)": 146.3, "J(Cb,Hb)": 129.6, "J(Ca,Hb)": -4.2, "J(Cb,Ha)": -5.1, "J(Ha,Hb)": 7.1,
                 "J(Ca,Cb)": 34.2}
        labeling = Labeling.enriched({"13C": 0.3})
        true_model = build_model(template("CH-CH3").with_couplings(truth), labeling=labeling)
        systems = [c.system for c in true_model.interpretation.components]
        obs = observe(systems, list(true_model.ratios), ranges=[(20.0, 60.0)] + RANGES)
        start = build_model(template("CH-CH3", one_bond={"Ca": 146.0, "Cb": 130.0}).with_couplings(
            {("Ca", "Cb"): 35.0}), labeling=labeling)
        res = refine(start.interpretation, obs, start.settings(RefineSettings(starts=1), shared_rate=True))
        got = start.named_couplings(res.parameters)
        for k, v in truth.items():
            self.assertLess(abs(got[k] - v), 0.02, msg=k)
        self.assertEqual(start.parents, [None, None, None])          # at 30 % every label set is a primary
        g = np.abs(res.gains)
        self.assertLess(abs(g[2] / g[0] - start.ratios[2] / start.ratios[0]), 0.05)


class RefinementTests(unittest.TestCase):
    def test_built_model_recovers_named_couplings(self):
        truth = {"J(Ca,Ha)": 146.3, "J(Cb,Hb)": 129.6, "J(Ca,Hb)": -4.2, "J(Cb,Ha)": -5.1, "J(Ha,Hb)": 7.1}
        true_model = build_model(template("CH-CH3").with_couplings(truth))
        systems = [c.system for c in true_model.interpretation.components]
        obs = observe(systems, [1.0, 1.0], ranges=RANGES)
        start = build_model(template("CH-CH3", one_bond={"Ca": 146.0, "Cb": 130.0}))
        settings = start.settings(RefineSettings(starts=1), fixed_ratios=True, shared_rate=True)
        res = refine(start.interpretation, obs, settings)
        got = start.named_couplings(res.parameters)
        for k, v in truth.items():
            self.assertLess(abs(got[k] - v), 0.01, msg=k)
        self.assertAlmostEqual(abs(res.gains[1]) / abs(res.gains[0]), 1.0, places=10)
        findings = run_checks(start, res.summary())
        self.assertEqual([f.code for f in findings], [])


class CheckTests(unittest.TestCase):
    def summary(self, gains, rates, params=None, hits=(), residual=0.2):
        p = {f"c{c}.log_rate0": float(np.log(r)) for c, r in enumerate(rates)}
        p.update(params or {"c0.J0-1": 146.0, "c0.J0-2": -4.5, "c0.J1-2": 7.0, "c1.J0-1": 130.0,
                            "c1.J0-2": -4.5, "c1.J1-2": 7.0})
        return {"gains": [[g, 0.0] for g in gains], "parameters": p, "boundary_hits": list(hits),
                "data_region_residual": residual}

    def test_rules(self):
        m = build_model(template("CH-CH3"))
        codes = lambda fs: sorted(f.code for f in fs)
        self.assertEqual(codes(run_checks(m, self.summary([1.0, 0.95], [2.0, 2.2]))), [])
        self.assertEqual(codes(run_checks(m, self.summary([1.0, 0.4], [2.0, 2.2]))), ["abundance"])
        # 4322bdfc: the CH component of a CH-CH3 fit switched off (0.061 of the methyl carbon)
        self.assertIn("collapsed_component", codes(run_checks(m, self.summary([0.061, 1.0], [2.0, 2.2]))))
        self.assertEqual(codes(run_checks(m, self.summary([1.0, 1.0], [1.9, 4.7]))), ["rate_asymmetry"])
        self.assertIn("background_component", codes(run_checks(m, self.summary([7.5, 1.0], [7.7, 1.1]))))
        self.assertEqual(codes(run_checks(m, self.summary([1.0, 1.0], [2.0, 2.0], hits=["c0.J0-1"]))), ["bounds"])
        peer = self.summary([1.0, 1.0], [2.0, 2.0], {"c0.J0-1": 146.0, "c0.J0-2": -2.4, "c0.J1-2": 7.0,
                                                    "c1.J0-1": 130.0, "c1.J0-2": -4.5, "c1.J1-2": 7.0}, residual=0.21)
        found = run_checks(m, self.summary([1.0, 1.0], [2.0, 2.0]), peers=[peer])
        self.assertEqual(codes(found), ["undetermined"])
        self.assertEqual(sorted(found[0].data["spread"]), ["J(Ca,Hb)"])


class KnowledgeTests(unittest.TestCase):
    def test_packaged_entries_and_nearest(self):
        kb = KnowledgeBase.packaged()
        self.assertEqual({e.compound for e in kb.entries_for("CH-CH3")}, {"L-alanine", "lactic acid"})
        query = {"J(Ca,Ha)": 147.0, "J(Cb,Hb)": 129.3, "J(Ha,Hb)": 6.6}
        self.assertEqual(kb.nearest("CH-CH3", query)[0][1].compound, "lactic acid")
        kb.add(ReferenceEntry("x", "CH-CH3", {"J(Ca,Ha)": 147.0}, "literature, test", source_kind="literature"))
        self.assertEqual(len(kb.nearest("CH-CH3", query)), 2)          # literature excluded by default
        self.assertEqual(len(kb.nearest("CH-CH3", query, kinds=None)), 3)
        blind = kb.nearest("CH-CH3", query, exclude=("b683220d",))           # the lactic-acid sample itself
        self.assertEqual([e.compound for _, e in blind], ["L-alanine"])
        self.assertEqual([e.compound for e in kb.entries_for("CH2-CH2")], ["ethylenediamine"])
        m = build_model(template("CH2-CH2"))
        self.assertEqual(m.component_labels, ["13C@C1 (x2)"])                # both carbons equivalent
        self.assertEqual([e.compound for e in kb.entries_for("N-ethyl (Et3N)")], ["triethylamine"])
        et3n = build_model(template("N-ethyl (Et3N)"))
        self.assertEqual(et3n.component_labels, ["13C@C1", "13C@C2"])       # N and the pseudo-site X never labelled


if __name__ == "__main__":
    unittest.main()


class Stage2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from zulf_hypothesis import propose_hypotheses
        truth = {"J(Ca,Ha)": 146.3, "J(Cb,Hb)": 129.6, "J(Ca,Hb)": -4.2, "J(Cb,Ha)": -5.1, "J(Ha,Hb)": 7.1}
        model = build_model(template("CH-CH3").with_couplings(truth))
        cls.obs = observe([c.system for c in model.interpretation.components], [1.0, 1.0], rate=1.0, noise=2.0,
                          ranges=RANGES)
        cls.propose = staticmethod(propose_hypotheses)

    def test_inventory_and_top_proposal(self):
        ps = self.propose(self.obs)
        peaks = [b.peak.frequency_hz for b in ps.inventory.bands]
        self.assertTrue(any(abs(p - 146.3) < 3 for p in peaks) and any(abs(p - 129.6) < 3 for p in peaks))
        top = ps.proposals[0]
        self.assertEqual(sorted(g.pattern for g in top.groups), ["13CH", "13CH3"])
        self.assertTrue(top.topology.startswith("bonded"))
        self.assertEqual(ps.models[0].ratios, (1.0, 1.0))
        self.assertTrue(any("2 x 13CH3" in p.topology for p in ps.proposals))      # isopropyl alternative kept

    def test_hints_steer_but_do_not_decide(self):
        from zulf_hypothesis import ProposerHints
        from zulf_core.spinsystem import Component, Interpretation

        good = build_model(template("CH-CH3", one_bond={"Ca": 146.0, "Cb": 130.0})).interpretation
        bad = Interpretation((Component(grouped(["13C", "1H"], [1, 1], [[0, 60.0], [60.0, 0]])),))
        wrong_type = Interpretation((Component(grouped(["13C", "1H"], [1, 2], [[0, 97.5], [97.5, 0]])),))

        class Fake:
            name = "fake"

            def propose(self, observed, k):
                return [good, bad, wrong_type][:k]

        ps = self.propose(self.obs, providers=[ProposerHints(Fake(), k=3)])
        self.assertIn("fake", ps.proposals[0].sources)
        self.assertEqual(sorted(g.pattern for g in ps.proposals[0].groups), ["13CH", "13CH3"])
        verdicts = [r["verdict"] for r in ps.insight]
        self.assertIn("hint predicts lines where the data show none", verdicts)
        hinted_only = [c for c in ps.candidates if "hint-only" in " ".join(c.notes)]
        self.assertTrue(all(c.score == 0.0 for c in hinted_only))
        self.assertEqual(len(ps.hinted_interpretations), 3)

    def test_trees_and_combined_models(self):
        from zulf_hypothesis.enumerate import _trees
        from zulf_hypothesis import combine_models
        self.assertEqual([len(_trees(n)) for n in (1, 2, 3, 4)], [1, 1, 3, 16])
        a = build_model(template("CH-CH3"))
        b = build_model(template("CH(CH3)2"))
        c = combine_models([a, b])
        self.assertEqual(c.blocks(), [[0, 1], [2, 3]])
        self.assertIn(("c2.J1-2", "c3.J1-2", "c3.J1-3"), c.ties())
        self.assertIn(("c2.log_rate0", "c3.log_rate0"), c.ties(shared_rate=True))
        self.assertNotIn(("c0.log_rate0", "c1.log_rate0", "c2.log_rate0", "c3.log_rate0"), c.ties(shared_rate=True))
        # Fixed ratios hold within each part only: one free amplitude per part.
        self.assertEqual(c.settings(fixed_ratios=True).amplitude_map,
                         ((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 2.0)))


class SearchTests(unittest.TestCase):
    def test_search_refines_extends_and_ranks(self):
        from zulf_hypothesis import SearchSettings, propose_hypotheses, search_hypotheses
        f = AddCoupledProton().propose(template("CH-CH3"))[0]
        truth = f.with_couplings({"J(Ca,Ha)": 146.3, "J(Cb,Hb)": 129.6, "J(Ca,Hb)": -4.2, "J(Cb,Ha)": -5.1,
                                  "J(Ha,Hb)": 7.1, "J(Ca,HX)": -2.0, "J(Cb,HX)": 1.5, "J(Ha,HX)": 5.0,
                                  "J(Hb,HX)": 0.3})
        tm = build_model(truth)
        obs = observe([c.system for c in tm.interpretation.components], [1.0, 1.0], rate=0.8, noise=2.0,
                      ranges=RANGES)
        ps = propose_hypotheses(obs)
        base = RefineSettings(starts=1, band_weighting="signal", background_order=1, max_seconds=120,
                              continuation_rates_per_s=(0.0,))
        res = search_hypotheses(obs, ps, SearchSettings(base=base, top_models=1, rounds=1, extend_top=1,
                                                        motif_screen=0,
                                                        extension_global_search={"max_seconds": 20.0, "solutions": 2,
                                                                                 "popsize": 6, "maxiter": 10}))
        steps = [l["step"] for l in res.log]
        self.assertIn("extend", steps)
        self.assertIn("accept", steps)
        self.assertIn("HX", res.best.name)                      # the extension with the coupled proton wins
        self.assertEqual(sorted(g for g in {e.variant for e in res.evaluated}), ["fixed", "free"])
        top = res.table()[0]
        self.assertEqual(top["delta"], 0.0)
        self.assertTrue(any("HX" in k for k in res.best.couplings))   # the added proton's couplings are reported
        self.assertTrue(res.best.knowledge and res.best.knowledge[0]["template"] == "CH-CH3")

    def test_demoted_fits_rank_last(self):
        from zulf_hypothesis.search import Evaluated, SearchResult
        m = build_model(template("CH-CH3"))
        summary = {"gains": [[1.0, 0.0], [1.0, 0.0]]}
        low = Evaluated("a", "free", "x", m, summary, np.zeros(1), score=10.0,
                        findings=[Finding("collapsed_component", "warn", "")])
        high = Evaluated("b", "fixed", "x", m, summary, np.zeros(1), score=50.0)
        res = SearchResult([low, high], [], None)
        self.assertEqual([e.name for e in res.ranked()], ["b", "a"])
        self.assertIs(res.best, high)
        rows = res.table()
        self.assertEqual((rows[1]["demoted"], rows[1]["delta"]), (True, -40.0))
        ext = Evaluated("b + HX", "fixed", "x", m, summary, np.zeros(1), score=48.0, status="rejected")
        res = SearchResult([low, high, ext], [], None)
        self.assertEqual([e.name for e in res.ranked()], ["b", "b + HX", "a"])   # e3d282da: rejected by 2.3
        self.assertIs(res.best, high)

    def test_motif_screen_picks_top_motifs(self):
        from zulf_hypothesis import SearchSettings, propose_hypotheses, search_hypotheses
        m = build_model(template("CH-CH3"))
        obs = observe([c.system for c in m.interpretation.components], [1.0, 1.0], rate=1.0, noise=2.0,
                      ranges=RANGES)
        ps = propose_hypotheses(obs)
        base = RefineSettings(starts=1, band_weighting="signal", background_order=1, max_seconds=60,
                              continuation_rates_per_s=(0.0,))
        res = search_hypotheses(obs, ps, SearchSettings(base=base, top_models=0, top_motifs=1, motif_screen=4,
                                                        rounds=0))
        screen = [l for l in res.log if l["step"] == "screen"]
        self.assertGreaterEqual(len(screen), 2)
        self.assertEqual(sum(l["kept"] for l in screen), 1)
        self.assertTrue(screen[0]["kept"])
        self.assertEqual({e.name for e in res.evaluated}, {screen[0]["motif"]})
        self.assertIn("CH-CH3", res.best.name)

    def test_yardstick_counts_and_noise(self):
        from zulf_hypothesis import yardstick
        m = build_model(template("CH-CH3"))
        obs = observe([c.system for c in m.interpretation.components], [1.0, 1.0], rate=1.0, noise=2.0,
                      ranges=RANGES)
        stick = yardstick(obs)
        self.assertEqual(stick.n, 2 * int(stick.mask.sum()))
        self.assertGreater(stick.sigma, 0.0)
        sel = obs.selected
        # Chi2 of the noise-free truth equals the noise contribution: of order n, not orders of magnitude off.
        clean = observe([c.system for c in m.interpretation.components], [1.0, 1.0], rate=1.0, noise=0.0,
                        ranges=RANGES)
        chi2 = stick.chi2(obs.values[sel], clean.values[clean.selected])
        self.assertLess(0.2 * stick.n, chi2)
        self.assertLess(chi2, 5.0 * stick.n)


class MotifTests(unittest.TestCase):
    def test_every_motif_builds_with_its_symmetry(self):
        from zulf_hypothesis import MOTIFS
        for name, motif in MOTIFS.items():
            model = build_model(motif.fragment({s.site: sum(s.j_range_hz) / 2 for s in motif.one_bond}))
            self.assertTrue(model.component_labels, msg=name)
        amine = build_model(MOTIFS["CH3-NH3+"].fragment({"C1": 145.0, "N1": 75.0}))
        self.assertEqual(amine.component_labels, ["13C@C1", "15N@N1"])
        n_system = amine.interpretation.components[1].system
        self.assertLess(n_system.couplings_hz[0, 1:].min(), -70.0)          # 1J(15N,H) negative
        benzene = build_model(MOTIFS["benzene ring"].fragment({"A2": 158.0}))
        self.assertEqual(len(benzene.component_labels), 1)                 # six equivalent carbons
        self.assertEqual(sum(len(g) for g in benzene.interpretation.components[0].system.groups), 7)
        diamine = build_model(MOTIFS["H2N-CH2-CH2-NH2"].fragment({"C1": 135.0, "N1": 68.0}))
        self.assertEqual(diamine.component_labels, ["13C@C1 (x2)", "15N@N1 (x2)"])  # symmetric ends merged

    def test_benchmark_cases_found_by_the_motif_scan(self):
        from zulf_hypothesis.benchmark import run_benchmark
        rows = {r["case"]: r for r in run_benchmark(["ethyl", "CH2-CH2", "isopropyl", "CH3-15NH3",
                                                          "ethylenediamine (fast exchange)"], verbose=False)}
        for case, row in rows.items():
            self.assertEqual(row["motif_rank"], 1, msg=f"{case}: {row}")
        self.assertIsNone(rows["ethyl"]["group_rank"])                      # the group path alone misses it
