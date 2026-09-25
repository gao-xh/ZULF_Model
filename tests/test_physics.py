import unittest
from fractions import Fraction

import numpy as np

from zulf_model.nuclei import get_registry
from zulf_model.physics import Protocol, Pulse, compute_transitions, merge_transitions, reference_signal
from zulf_model.physics.operators import angular_momentum, collective_multiplicities
from zulf_model.spinsystem import SpinSystem


def xh_n(n_h, j=140.0):
    iso = ("13C",) + ("1H",) * n_h
    m = np.zeros((n_h + 1, n_h + 1))
    m[0, 1:] = m[1:, 0] = j
    return SpinSystem(iso, m)


def random_system(isotopes, seed, scale=60.0):
    rng = np.random.default_rng(seed)
    n = len(isotopes)
    j = np.triu(rng.uniform(-scale, scale, (n, n)), 1)
    return SpinSystem(tuple(isotopes), j + j.T)


class OperatorTests(unittest.TestCase):
    def test_commutation_and_casimir(self):
        for spin in (Fraction(1, 2), Fraction(1), Fraction(3, 2)):
            sx, sy, sz = angular_momentum(spin)
            np.testing.assert_allclose(sx @ sy - sy @ sx, 1j * sz, atol=1e-12)
            s = float(spin)
            np.testing.assert_allclose(sx @ sx + sy @ sy + sz @ sz, s * (s + 1) * np.eye(len(sz)), atol=1e-12)

    def test_collective_multiplicities(self):
        self.assertEqual(collective_multiplicities(3, Fraction(1, 2)), ((Fraction(3, 2), 1), (Fraction(1, 2), 2)))
        self.assertEqual(collective_multiplicities(2, Fraction(1)),
                         ((Fraction(2), 1), (Fraction(1), 1), (Fraction(0), 1)))


class ClosedFormTests(unittest.TestCase):
    """Known zero-field XA_n patterns (J = 1J_CH)."""

    def test_xa_lines(self):
        j = 137.0
        cases = {1: [j], 2: [1.5 * j], 3: [j, 2 * j]}
        for n_h, expected in cases.items():
            tl = compute_transitions(xh_n(n_h, j))
            np.testing.assert_allclose(tl.frequencies_hz, expected, atol=1e-9)
            self.assertTrue(np.all(np.abs(tl.amplitudes.imag) < 1e-12))

    def test_homonuclear_pair_has_no_signal(self):
        s = SpinSystem(("1H", "1H"), np.array([[0, 7.0], [7.0, 0]]))
        self.assertEqual(len(compute_transitions(s)), 0)

    def test_negative_coupling_same_frequency_magnitude(self):
        pos = compute_transitions(xh_n(1, 90.0))
        neg = compute_transitions(xh_n(1, -90.0))
        np.testing.assert_allclose(pos.frequencies_hz, neg.frequencies_hz)


class PropagationAgreementTests(unittest.TestCase):
    """Transition lists must reproduce brute-force expm propagation."""

    def check(self, system, protocol, t_max=0.04):
        t = np.linspace(0, t_max, 31)
        ref = reference_signal(system, t, protocol)
        scale = np.abs(ref).max()
        for method in ("full", "sectors"):
            tl = compute_transitions(system, protocol, method=method)
            err = np.abs(tl.signal(t) - ref).max() / scale
            self.assertLess(err, 1e-10, msg=f"{method} {protocol.name}")

    def test_random_chn_systems(self):
        for seed, iso in enumerate([("1H", "13C", "1H", "15N"), ("13C", "1H", "1H", "1H", "13C"),
                                     ("15N", "1H", "1H", "13C", "1H", "1H")]):
            self.check(random_system(iso, seed), Protocol())

    def test_equivalent_groups(self):
        iso = ("1H",) * 7 + ("13C",)
        j = np.zeros((8, 8))
        for i in range(6):
            j[i, 6] = j[6, i] = 6.5
        for i in range(3):
            j[i, 7] = j[7, i] = 125.0
        for i in range(3, 6):
            j[i, 7] = j[7, i] = 2.0
        j[6, 7] = j[7, 6] = -3.0
        system = SpinSystem(iso, j)
        self.assertEqual(len(system.groups), 4)
        self.check(system, Protocol(), t_max=0.02)

    def test_pulse_protocol(self):
        proto = Protocol(name="pulse", pulses=(Pulse(("1H",), np.pi / 2, "x"),))
        self.check(random_system(("1H", "13C", "1H"), 7), proto)
        proto_y = Protocol(name="pulse_y", pulses=(Pulse(("13C",), np.pi / 3, "y"),))
        self.check(random_system(("1H", "13C", "1H"), 8), proto_y)

    def test_nonuniform_intragroup(self):
        iso = ("1H",) * 7 + ("13C",)
        j = np.zeros((8, 8))
        for a in range(6):
            for b in range(a + 1, 6):
                j[a, b] = j[b, a] = -12.0 if (a < 3) == (b < 3) else 0.3
        j[:6, 6] = j[6, :6] = 6.5
        j[:6, 7] = j[7, :6] = -3.0
        j[6, 7] = j[7, 6] = 135.0
        system = SpinSystem(iso, j)
        self.assertEqual(len(system.groups), 3)
        self.check(system, Protocol(), t_max=0.02)

    def test_spin_one_nucleus(self):
        self.check(random_system(("14N", "1H", "1H"), 9), Protocol())

    def test_custom_weights(self):
        proto = Protocol(name="h_only", preparation={"1H": 1.0}, detection="gamma")
        self.check(random_system(("1H", "13C", "15N"), 10), proto)


class MergeTests(unittest.TestCase):
    def test_merge_sums_coincident(self):
        f, a = merge_transitions(np.array([10.0, 10.0 + 1e-9, 20.0]), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(f, [10.0 + 2e-9 / 3, 20.0], atol=1e-12)
        np.testing.assert_allclose(a, [3.0, 3.0])

    def test_merge_keeps_weak_lines(self):
        f, a = merge_transitions(np.array([10.0, 20.0]), np.array([1.0, 1e-8]))
        self.assertEqual(len(f), 2)

    def test_gamma_registry(self):
        self.assertLess(get_registry().gamma("15N"), 0)


if __name__ == "__main__":
    unittest.main()
