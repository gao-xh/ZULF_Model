import unittest

import numpy as np

from zulf_model.spec import JBinSpec, ProblemSpec
from zulf_core.spinsystem import (Component, Interpretation, SpinSystem, best_permutation,
                                   canonical_permutation, canonicalize, detect_equivalence_groups)


def isopropyl_like(methyl_c=125.0):
    iso = ("1H",) * 7 + ("13C",)
    j = np.zeros((8, 8))
    for i in range(6):
        j[i, 6] = j[6, i] = 6.5
    for i in range(3):
        j[i, 7] = j[7, i] = methyl_c
    for i in range(3, 6):
        j[i, 7] = j[7, i] = 2.0
    j[6, 7] = j[7, 6] = -3.0
    return SpinSystem(iso, j)


class SpinSystemTests(unittest.TestCase):
    def test_equivalence_detection(self):
        s = isopropyl_like()
        self.assertEqual(s.groups, ((0, 1, 2), (3, 4, 5), (6,), (7,)))
        self.assertEqual(s.group_signature(), (("13C", 1), ("1H", 1), ("1H", 3), ("1H", 3)))

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            SpinSystem(("1H", "1H"), np.array([[0, 1], [2, 0]]))
        with self.assertRaises(ValueError):
            SpinSystem(("1H", "1H"), np.array([[1, 1], [1, 0]]))
        with self.assertRaises(KeyError):
            SpinSystem(("99X",), np.zeros((1, 1)))
        s = isopropyl_like()
        with self.assertRaises(ValueError):
            SpinSystem(s.isotopes, s.couplings_hz, ((0, 1, 2, 3), (4, 5), (6,), (7,)))

    def test_permutation_matching_recovers_relabelling(self):
        rng = np.random.default_rng(3)
        iso = ("1H",) * 5 + ("13C", "15N")
        j = np.triu(rng.uniform(-15, 150, (7, 7)), 1)
        s = SpinSystem(iso, j + j.T)
        order = rng.permutation(7).tolist()
        shuffled = s.permute(order)
        match = best_permutation(s, shuffled)
        self.assertTrue(match.exact)
        self.assertAlmostEqual(match.rms_error_hz, 0.0, places=12)
        np.testing.assert_allclose(shuffled.permute(match.permutation).couplings_hz, s.couplings_hz)

    def test_global_sign_equivalence(self):
        s = isopropyl_like()
        flipped = SpinSystem(s.isotopes, -s.couplings_hz)
        self.assertEqual(best_permutation(s, flipped).sign, -1)
        self.assertAlmostEqual(best_permutation(s, flipped).rms_error_hz, 0.0)
        self.assertGreater(best_permutation(s, flipped, allow_global_sign=False).rms_error_hz, 1.0)
        c1 = canonicalize(s, ("13C", "1H"))
        c2 = canonicalize(flipped, ("13C", "1H"))
        np.testing.assert_allclose(c1.observable_couplings(), c2.observable_couplings())

    def test_matching_rejects_different_composition(self):
        a = isopropyl_like()
        b = SpinSystem(("1H",) * 6 + ("13C", "13C"), np.zeros((8, 8)))
        self.assertFalse(best_permutation(a, b).matched)

    def test_matching_brute_force_agreement(self):
        from zulf_core.spinsystem import all_isotope_preserving_permutations
        rng = np.random.default_rng(5)
        iso = ("1H", "1H", "1H", "13C", "1H")
        a = rng.normal(size=(5, 5)); a = np.triu(a, 1); a = a + a.T
        noise = np.triu(rng.normal(scale=0.3, size=(5, 5)), 1)
        b = a + noise + noise.T
        sa, sb = SpinSystem(iso, a), SpinSystem(iso, b)
        best = min(((a - b[np.ix_(p, p)]) ** 2).sum() for p in all_isotope_preserving_permutations(iso))
        m = best_permutation(sa, sb, allow_global_sign=False)
        got = ((a - b[np.ix_(m.permutation, m.permutation)]) ** 2).sum()
        self.assertAlmostEqual(got, best, places=10)

    def test_canonical_order_is_label_invariant_and_group_contiguous(self):
        s = isopropyl_like()
        rng = np.random.default_rng(1)
        for _ in range(5):
            t = s.permute(rng.permutation(8).tolist())
            c1, c2 = canonicalize(s, ("13C", "15N", "1H")), canonicalize(t, ("13C", "15N", "1H"))
            np.testing.assert_allclose(c1.couplings_hz, c2.couplings_hz)
            self.assertEqual(c1.isotopes, c2.isotopes)
        perm = canonical_permutation(s, ("13C", "1H"))
        self.assertEqual(s.isotopes[perm[0]], "13C")

    def test_group_expansion_roundtrip(self):
        s = isopropyl_like()
        order = [g for g in s.groups]
        gj = s.group_couplings()
        rebuilt = SpinSystem.from_group_couplings([s.isotopes[g[0]] for g in order], [len(g) for g in order], gj)
        self.assertLess(best_permutation(s, rebuilt).rms_error_hz, 1e-12)

    def test_json_roundtrip(self):
        s = isopropyl_like()
        interp = Interpretation((Component(s, 0.0107, "methyl"), Component(isopropyl_like(130), 0.02)))
        again = Interpretation.from_dict(interp.to_dict())
        np.testing.assert_allclose(again.components[0].system.couplings_hz, s.couplings_hz)
        self.assertEqual(again.components[1].contribution, 0.02)
        self.assertEqual(interp.key(("13C", "1H")), again.key(("13C", "1H")))

    def test_detect_groups_for_random_matrix_are_singletons(self):
        rng = np.random.default_rng(2)
        j = np.triu(rng.normal(size=(4, 4)), 1)
        self.assertEqual(detect_equivalence_groups(("1H", "1H", "1H", "13C"), j + j.T), ((0,), (1,), (2,), (3,)))

    def test_homonuclear_system_is_one_unobservable_group(self):
        rng = np.random.default_rng(2)
        j = np.triu(rng.normal(size=(4, 4)), 1)
        self.assertEqual(detect_equivalence_groups(("1H",) * 4, j + j.T), ((0, 1, 2, 3),))

    def test_nonuniform_intragroup_couplings_are_unobservable(self):
        iso = ("1H",) * 7 + ("13C",)
        j = np.zeros((8, 8))
        for a in range(6):
            for b in range(a + 1, 6):
                j[a, b] = j[b, a] = -12.0 if (a < 3) == (b < 3) else 0.3
        j[:6, 6] = j[6, :6] = 6.5
        j[:6, 7] = j[7, :6] = -3.0
        j[6, 7] = j[7, 6] = 135.0
        s = SpinSystem(iso, j)
        self.assertEqual(s.groups, ((0, 1, 2, 3, 4, 5), (6,), (7,)))
        self.assertEqual(s.observable_couplings()[0, 1], 0.0)
        flat = SpinSystem(iso, s.observable_couplings())
        self.assertEqual(best_permutation(s, flat).rms_error_hz, 0.0)


class SpecTests(unittest.TestCase):
    def test_defaults_and_roundtrip(self):
        spec = ProblemSpec()
        self.assertEqual(spec.max_spins, 8)
        self.assertEqual(spec.max_pairs, 28)
        again = ProblemSpec.from_dict(spec.to_dict())
        self.assertEqual(again, spec)
        self.assertEqual(again.digest(), spec.digest())

    def test_spin_counts_drive_dimensions(self):
        spec = ProblemSpec(spin_counts=(6, 8, 10))
        self.assertEqual(spec.max_spins, 10)
        self.assertTrue(spec.allows(("1H",) * 9 + ("13C",)))
        self.assertFalse(spec.allows(("1H",) * 8 + ("13C",)))

    def test_j_bins(self):
        edges = JBinSpec(((2.0, 0.5), (4.0, 1.0))).edges()
        np.testing.assert_allclose(edges, [-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4])
        with self.assertRaises(ValueError):
            JBinSpec(((2.0, 0.3),)).edges()
        with self.assertRaises(ValueError):
            ProblemSpec(nuclei=("1H", "1H"))


class PackageBoundaryTests(unittest.TestCase):
    def test_core_imports_without_torch_or_model_package(self):
        import subprocess
        import sys
        from pathlib import Path
        code = ("import sys\n"
                "import zulf_core.physics, zulf_core.render, zulf_core.solver, zulf_core.evaluation\n"
                "import zulf_core.io, zulf_core.diagnostics\n"
                "bad = sorted(m for m in sys.modules if m == 'torch' or m.startswith('torch.') or m.startswith('zulf_model'))\n"
                "print(','.join(bad))\n")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300,
                             cwd=str(Path(__file__).resolve().parents[1]))
        self.assertEqual(out.returncode, 0, msg=out.stderr)
        self.assertEqual(out.stdout.strip(), "")


if __name__ == "__main__":
    unittest.main()
