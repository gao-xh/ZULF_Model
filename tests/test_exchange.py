import unittest

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_core.physics.exchange import (exchange_transition_derivatives, exchange_transitions,
                                       reference_exchange_signal)
from zulf_core.physics.transitions import TransitionList
from zulf_core.render.renderer import ContinuousRenderer
from zulf_core.solver import ObservedSpectrum, ParameterPolicy, Parameterization
from zulf_core.solver.forward import MixtureForward
from zulf_core.spinsystem import Component, Interpretation, SpinSystem


def small_system():
    # 13C with a strongly coupled proton and one exchanging proton (index 2)
    j = np.zeros((3, 3))
    j[0, 1] = j[1, 0] = 140.0
    j[0, 2] = j[2, 0] = -4.0
    j[1, 2] = j[2, 1] = 7.0
    return SpinSystem(("13C", "1H", "1H"), j, groups=[(0,), (1,), (2,)])


def full_signal(tl, t):
    m = tl.metadata
    return tl.signal(t) + sum(a * np.exp(-r * t) for r, a in zip(m["nonoscillating_rates"],
                                                                 m["nonoscillating_amplitudes"]))


class ExchangeTest(unittest.TestCase):
    def test_matches_full_liouville_propagation(self):
        system = small_system()
        t = np.linspace(0.0, 0.15, 121)
        for k in (0.0, 2.0, 25.0, 400.0):
            tl = exchange_transitions(system, {2: k})
            ref = reference_exchange_signal(system, {2: k}, t)
            np.testing.assert_allclose(full_signal(tl, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_unpolarized_solvent(self):
        system = small_system()
        t = np.linspace(0.0, 0.1, 81)
        tl = exchange_transitions(system, {2: 30.0}, solvent_weight=0.0)
        ref = reference_exchange_signal(system, {2: 30.0}, t, solvent_weight=0.0)
        np.testing.assert_allclose(full_signal(tl, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_zero_rate_is_static_system(self):
        system = small_system()
        tl = exchange_transitions(system, {2: 0.0})
        ref = compute_transitions(system)
        np.testing.assert_allclose(tl.frequencies_hz, ref.frequencies_hz, atol=1e-6)
        np.testing.assert_allclose(tl.amplitudes, ref.amplitudes, atol=1e-9 * np.abs(ref.amplitudes).max())
        self.assertLess(tl.rates_of_lines().max(), 1e-8)

    def test_fast_exchange_decouples(self):
        system = small_system()
        decoupled = compute_transitions(SpinSystem(("13C", "1H"), system.couplings_hz[:2, :2]))
        tl = exchange_transitions(system, {2: 1.0e5})
        strong = np.abs(tl.amplitudes) > 1e-3 * np.abs(tl.amplitudes).max()
        np.testing.assert_allclose(tl.frequencies_hz[strong], decoupled.frequencies_hz, atol=1e-3)
        np.testing.assert_allclose(np.abs(tl.amplitudes[strong]), np.abs(decoupled.amplitudes), rtol=1e-3)
        # residual broadening of order J^2 / k
        self.assertLess(tl.rates_of_lines()[strong].max(), 50.0 ** 2 / 1.0e5 * 10)

    def test_line_rates_add_to_fitted_rates(self):
        tl = TransitionList(np.array([10.0, 20.0]), np.array([1.0, 0.5j]), line_rates=np.array([0.5, 2.0]))
        plain = TransitionList(np.array([10.0, 20.0]), np.array([1.0, 0.5j]))
        f = np.linspace(5, 25, 200)
        r = ContinuousRenderer()
        a = r.render(tl, 1.0, f)
        b = r.render(plain.within(9, 11), 1.5, f) + r.render(plain.within(19, 21), 3.0, f)
        np.testing.assert_allclose(a, b, atol=1e-12)

    def test_derivatives_match_differences(self):
        # C with two equivalent protons (degenerate modes) and one exchanging proton; independent reference:
        # central differences of the oscillating signal of recomputed lists
        iso = ("13C", "1H", "1H", "1H")
        j = np.zeros((4, 4))
        j[0, 1] = j[1, 0] = j[0, 2] = j[2, 0] = 130.0
        j[1, 2] = j[2, 1] = -12.0
        j[0, 3] = j[3, 0] = -3.0
        j[1, 3] = j[3, 1] = j[2, 3] = j[3, 2] = 6.5
        groups = [(0,), (1, 2), (3,)]
        k0 = 15.0
        t = np.linspace(0.0, 0.3, 300)

        def osc(jm, k):
            tl = exchange_transitions(SpinSystem(iso, jm, groups=groups), {3: k})
            return tl.signal(t) - np.real(tl.dc)

        directions = [[(0, 1), (0, 2)], [(0, 3)], [(1, 3), (2, 3)]]
        d = exchange_transition_derivatives(SpinSystem(iso, j, groups=groups), {3: k0}, directions, [[3]])
        tl = d.transitions
        lam = -tl.rates_of_lines() + 2j * np.pi * tl.frequencies_hz

        def analytic(k):
            return np.real(((d.d_amplitudes[k][None, :] + d.t_terms[k][None, :] * t[:, None])
                            * np.exp(lam[None, :] * t[:, None])).sum(axis=1))
        for k, pairs in enumerate(directions):
            h = 1e-4
            jp, jm = j.copy(), j.copy()
            for p, q in pairs:
                jp[p, q] += h
                jp[q, p] += h
                jm[p, q] -= h
                jm[q, p] -= h
            num = (osc(jp, k0) - osc(jm, k0)) / (2 * h)
            np.testing.assert_allclose(analytic(k), num, atol=1e-6 * np.abs(num).max())
        h = 1e-5
        num = (osc(j, k0 * np.exp(h)) - osc(j, k0 * np.exp(-h))) / (2 * h)
        np.testing.assert_allclose(analytic(3), num, atol=1e-6 * np.abs(num).max())

    def test_torch_backend_matches_numpy(self):
        try:
            import torch  # noqa: F401
        except ImportError:
            self.skipTest("torch not installed")
        from zulf_core.physics import exchange
        system = small_system()
        before = exchange.get_backend()
        ref = exchange_transitions(system, {2: 25.0})
        try:
            exchange.set_backend("torch", "cpu")
            tl = exchange_transitions(system, {2: 25.0})
        finally:
            exchange.set_backend(before["name"], before["device"])
        np.testing.assert_allclose(tl.frequencies_hz, ref.frequencies_hz, atol=1e-8)
        np.testing.assert_allclose(tl.amplitudes, ref.amplitudes, atol=1e-10 * np.abs(ref.amplitudes).max())
        np.testing.assert_allclose(tl.rates_of_lines(), ref.rates_of_lines(), atol=1e-8)


class ExchangeSolverTest(unittest.TestCase):
    def setUp(self):
        self.system = small_system()
        self.f = np.linspace(100.0, 180.0, 1601)
        tl = exchange_transitions(self.system, {2: 20.0})
        self.values = ContinuousRenderer().render(tl, 0.8, self.f, gain=np.exp(0.4j))
        self.obs = ObservedSpectrum.from_spectrum(self.f, self.values, [(100.0, 180.0)])

    def param(self, k):
        interp = Interpretation((Component(self.system, 1.0, "x"),))
        p = Parameterization.from_interpretation(interp, ParameterPolicy(fit_phase_delay=False,
                                                                          initial_rate_per_s=0.8))
        p.add_exchange(0, 2, k)
        return p

    def test_true_rate_fits_exactly(self):
        p = self.param(20.0)
        fw = MixtureForward(p, self.obs)
        self.assertLess(np.linalg.norm(fw.predict(p.vector()).residual), 1e-8)

    def test_jacobian_matches_residual_differences(self):
        p = self.param(35.0)
        fw = MixtureForward(p, self.obs)
        x = p.vector()
        x[0] += 0.7                       # move a coupling off the truth
        jac = fw.jacobian(x)
        num = np.zeros_like(jac)
        for i in range(len(x)):
            h = 1e-6 * max(1.0, abs(x[i]))
            e = np.zeros(len(x))
            e[i] = h
            num[:, i] = (fw.predict(x + e).residual - fw.predict(x - e).residual) / (2 * h)
        np.testing.assert_allclose(jac, num, atol=1e-5 * np.abs(num).max())

    def test_refine_recovers_rate(self):
        from zulf_core.solver import RefineSettings, refine
        p = self.param(4.0)
        interp = Interpretation((Component(self.system, 1.0, "x"),))
        settings = RefineSettings(policy=ParameterPolicy(fit_phase_delay=False, initial_rate_per_s=0.8), starts=1)
        result = refine(interp, self.obs, settings, parameterization=p)
        self.assertAlmostEqual(np.exp(result.parameters["c0.log_kex2"]), 20.0, delta=0.2)


if __name__ == "__main__":
    unittest.main()


class GroupReductionTest(unittest.TestCase):
    """Non-exchanging equivalence groups replaced by their total spins (exchange_transitions(reduce=True))."""

    @staticmethod
    def methyl_nh():
        # 13C-H3 (one equivalent group) next to an exchanging N-H proton
        iso = ("13C", "1H", "1H", "1H", "1H")
        j = np.zeros((5, 5))
        for h in (1, 2, 3):
            j[0, h] = j[h, 0] = 128.0
            j[h, 4] = j[4, h] = 6.0
        j[0, 4] = j[4, 0] = -2.5
        return SpinSystem(iso, j, groups=[(0,), (1, 2, 3), (4,)])

    def test_matches_brute_force(self):
        system = self.methyl_nh()
        t = np.linspace(0.0, 0.12, 97)
        for k in (0.0, 3.0, 40.0):
            tl = exchange_transitions(system, {4: k})
            self.assertEqual(tl.metadata["method"], "exchange (group-reduced)")
            ref = reference_exchange_signal(system, {4: k}, t)
            np.testing.assert_allclose(full_signal(tl, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_matches_unreduced(self):
        # 13C, CH2, CH3 and one exchanging N-H proton (7 spins)
        iso = ("13C",) + ("1H",) * 6
        j = np.zeros((7, 7))
        for h in (1, 2):
            j[0, h] = 135.0
        for h in (3, 4, 5):
            j[0, h] = -4.5
            for g in (1, 2):
                j[g, h] = 7.1
        j[0, 6] = 1.5
        for g in (1, 2):
            j[g, 6] = 5.8
        for h in (3, 4, 5):
            j[h, 6] = 0.4
        j = j + j.T
        system = SpinSystem(iso, j, groups=[(0,), (1, 2), (3, 4, 5), (6,)])
        t = np.linspace(0.0, 0.2, 161)
        rates = {6: 8.0}
        reduced = exchange_transitions(system, rates)
        full = exchange_transitions(system, rates, reduce=False)
        self.assertLess(max(reduced.metadata["subspace_dimension"]), full.metadata["subspace_dimension"])
        ref = full_signal(full, t)
        np.testing.assert_allclose(full_signal(reduced, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_derivatives_of_reduced_blocks(self):
        system = self.methyl_nh()
        j = system.couplings_hz.copy()
        groups = system.groups
        t = np.linspace(0.0, 0.3, 300)
        k0 = 12.0

        def osc(jm, k):
            tl = exchange_transitions(SpinSystem(system.isotopes, jm, groups=groups), {4: k})
            return tl.signal(t) - np.real(tl.dc)

        directions = [[(0, 1), (0, 2), (0, 3)], [(1, 4), (2, 4), (3, 4)], [(0, 4)]]
        d = exchange_transition_derivatives(system, {4: k0}, directions, [[4]])
        tl = d.transitions
        lam = -tl.rates_of_lines() + 2j * np.pi * tl.frequencies_hz

        def analytic(k):
            return np.real(((d.d_amplitudes[k][None, :] + d.t_terms[k][None, :] * t[:, None])
                            * np.exp(lam[None, :] * t[:, None])).sum(axis=1))
        for k, pairs in enumerate(directions):
            h = 1e-4
            jp, jm = j.copy(), j.copy()
            for p, q in pairs:
                jp[p, q] += h
                jp[q, p] += h
                jm[p, q] -= h
                jm[q, p] -= h
            num = (osc(jp, k0) - osc(jm, k0)) / (2 * h)
            np.testing.assert_allclose(analytic(k), num, atol=1e-6 * np.abs(num).max())
        h = 1e-5
        num = (osc(j, k0 * np.exp(h)) - osc(j, k0 * np.exp(-h))) / (2 * h)
        np.testing.assert_allclose(analytic(3), num, atol=1e-6 * np.abs(num).max())
        with self.assertRaises(ValueError):          # a direction moving only part of a reduced group
            exchange_transition_derivatives(system, {4: k0}, [[(0, 1)]], [])


def cn_system():
    # 13C bonded to one proton and to a 14N (spin 1): 1J(C,H), 1J(C,14N), 2J(H,14N)
    j = np.zeros((3, 3))
    j[0, 1] = j[1, 0] = 140.0
    j[0, 2] = j[2, 0] = 4.0
    j[1, 2] = j[2, 1] = 2.0
    return SpinSystem(("13C", "1H", "14N"), j, groups=[(0,), (1,), (2,)])


class QuadrupolarRelaxationTest(unittest.TestCase):
    """14N (spin 1) with isotropic quadrupolar relaxation. Reference: the full Liouville space with the dense
    double commutator over the five standard rank-2 tensor operators (reference_exchange_signal)."""

    def test_matches_full_liouville_propagation(self):
        system = cn_system()
        t = np.linspace(0.0, 0.1, 81)
        for r in (0.0, 5.0, 60.0, 600.0):
            tl = exchange_transitions(system, {}, quadrupolar_per_s={2: r})
            ref = reference_exchange_signal(system, {}, t, quadrupolar_per_s={2: r})
            np.testing.assert_allclose(full_signal(tl, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_zero_rate_is_static_system(self):
        system = cn_system()
        t = np.linspace(0.0, 0.2, 101)
        tl = exchange_transitions(system, {}, quadrupolar_per_s={2: 0.0})
        static = compute_transitions(system)
        np.testing.assert_allclose(full_signal(tl, t), static.signal(t), atol=1e-9 * np.abs(static.signal(t)).max())

    def test_fast_relaxation_decouples_with_broadening_inverse_in_rate(self):
        # R >> 2 pi J: the C-H line returns to J(C,H) (14N decoupled) and its residual broadening (scalar relaxation
        # of the second kind) falls as 1 / R
        system = cn_system()
        widths = []
        for r in (3.0e3, 3.0e4):
            tl = exchange_transitions(system, {}, quadrupolar_per_s={2: r})
            i = int(np.argmin(np.abs(tl.frequencies_hz - 140.0)))
            self.assertAlmostEqual(tl.frequencies_hz[i], 140.0, delta=0.02)
            widths.append(tl.rates_of_lines()[i])
        self.assertAlmostEqual(widths[0] / widths[1], 10.0, delta=0.3)

    def test_rank1_rate_normalization(self):
        # a lone relaxing 14N: its I_z (the only rank-1 q = 0 operator of one spin) decays exactly at R
        from fractions import Fraction
        from zulf_core.physics.exchange import _structure
        st = _structure((Fraction(1), Fraction(1, 2)))
        q = st.quadrupolar(0)
        np.testing.assert_allclose(q @ st.iz[0], -st.iz[0], atol=1e-12)
        np.testing.assert_allclose(q @ st.iz[1], 0.0, atol=1e-12)
        with self.assertRaises(ValueError):
            st.quadrupolar(1)                                   # spin 1/2 has no quadrupole

    def test_with_exchange_and_group_reduction(self):
        # N-H proton exchanging with the solvent, 14N relaxing, a CH2 group reduced to its total spins
        iso = ("13C", "1H", "1H", "14N", "1H")
        j = np.zeros((5, 5))
        j[0, 1] = j[1, 0] = j[0, 2] = j[2, 0] = 132.0
        j[1, 2] = j[2, 1] = -12.0
        j[0, 3] = j[3, 0] = 4.0
        j[1, 3] = j[3, 1] = j[2, 3] = j[3, 2] = 1.5
        j[0, 4] = j[4, 0] = -3.0
        j[1, 4] = j[4, 1] = j[2, 4] = j[4, 2] = 5.0
        j[3, 4] = j[4, 3] = 50.0
        system = SpinSystem(iso, j, groups=[(0,), (1, 2), (3,), (4,)])
        t = np.linspace(0.0, 0.08, 9)                         # the dense reference costs one 2305^2 expm per time
        tl = exchange_transitions(system, {4: 20.0}, quadrupolar_per_s={3: 150.0})
        self.assertEqual(tl.metadata["method"], "exchange (group-reduced)")
        ref = reference_exchange_signal(system, {4: 20.0}, t, quadrupolar_per_s={3: 150.0})
        np.testing.assert_allclose(full_signal(tl, t), ref, atol=1e-9 * np.abs(ref).max())

    def test_derivatives_match_differences(self):
        system = cn_system()
        j0 = system.couplings_hz.copy()
        r0 = 80.0
        t = np.linspace(0.0, 0.3, 300)

        def osc(jm, r):
            tl = exchange_transitions(SpinSystem(system.isotopes, jm, groups=system.groups), {},
                                      quadrupolar_per_s={2: r})
            return tl.signal(t) - np.real(tl.dc)
        directions = [[(0, 1)], [(0, 2)], [(1, 2)]]
        d = exchange_transition_derivatives(system, {}, directions, quadrupolar_per_s={2: r0},
                                            quadrupolar_spins=[[2]])
        tl = d.transitions
        lam = -tl.rates_of_lines() + 2j * np.pi * tl.frequencies_hz

        def analytic(k):
            return np.real(((d.d_amplitudes[k][None, :] + d.t_terms[k][None, :] * t[:, None])
                            * np.exp(lam[None, :] * t[:, None])).sum(axis=1))
        for k, pairs in enumerate(directions):
            h = 1e-4
            jp, jm = j0.copy(), j0.copy()
            for p, q in pairs:
                jp[p, q] += h
                jp[q, p] += h
                jm[p, q] -= h
                jm[q, p] -= h
            num = (osc(jp, r0) - osc(jm, r0)) / (2 * h)
            np.testing.assert_allclose(analytic(k), num, atol=1e-6 * np.abs(num).max())
        h = 1e-5
        num = (osc(j0, r0 * np.exp(h)) - osc(j0, r0 * np.exp(-h))) / (2 * h)
        np.testing.assert_allclose(analytic(3), num, atol=1e-6 * np.abs(num).max())
