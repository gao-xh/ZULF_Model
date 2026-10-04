import unittest

import numpy as np

from zulf_core.physics.lines import line_table, predict_lines
from zulf_core.physics.transitions import compute_transitions
from zulf_core.spinsystem import SpinSystem


def system(j):
    """13C (group 0), one proton (group 1), three methyl protons (group 2)."""
    g = np.array([[0.0, j[0], j[1]], [j[0], 0.0, j[2]], [j[1], j[2], 0.0]])
    return SpinSystem.from_group_couplings(["13C", "1H", "1H"], [1, 1, 3], g)


class LineTableTests(unittest.TestCase):
    def test_closed_forms_of_xa_and_xa3(self):
        # zero field: 13C-H gives one line at J; 13C-H3 gives lines at J and 2J
        xa = SpinSystem.from_group_couplings(["13C", "1H"], [1, 1], np.array([[0.0, 140.0], [140.0, 0.0]]))
        lines = line_table(xa)
        self.assertEqual(len(lines), 1)
        self.assertAlmostEqual(lines[0].frequency_hz, 140.0, places=6)
        self.assertAlmostEqual(lines[0].df_dj["J(0,1)"], 1.0, places=6)
        xa3 = SpinSystem.from_group_couplings(["13C", "1H"], [1, 3], np.array([[0.0, 125.0], [125.0, 0.0]]))
        f = sorted(round(line.frequency_hz, 6) for line in line_table(xa3))
        self.assertEqual(f, [125.0, 250.0])
        self.assertEqual(sorted(round(line.df_dj["J(0,1)"], 6) for line in line_table(xa3)), [1.0, 2.0])

    def test_contributions_sum_to_the_frequency(self):
        # Euler: the Hamiltonian is homogeneous of degree 1 in J, so f = sum_k J_k df/dJ_k
        for line in line_table(system([132.0, -4.5, 6.8])):
            self.assertAlmostEqual(sum(line.contribution_hz.values()), line.frequency_hz, places=6)

    def test_derivatives_against_finite_differences_of_the_transition_list(self):
        j0 = np.array([132.0, -4.5, 6.8])
        lines = [line for line in line_table(system(j0), min_relative_amplitude=0.05) if line.frequency_hz > 50]
        h = 1e-4
        names = ["J(0,1)", "J(0,2)", "J(1,2)"]
        for k, name in enumerate(names):
            jp, jm = j0.copy(), j0.copy()
            jp[k] += h
            jm[k] -= h
            fp = np.asarray(compute_transitions(system(jp)).frequencies_hz)
            fm = np.asarray(compute_transitions(system(jm)).frequencies_hz)
            for line in lines:
                ip, im = np.argmin(np.abs(fp - line.frequency_hz)), np.argmin(np.abs(fm - line.frequency_hz))
                numeric = (fp[ip] - fm[im]) / (2 * h)
                self.assertAlmostEqual(line.df_dj[name], numeric, places=4, msg=(name, line.frequency_hz))

    def test_zero_and_tied_couplings(self):
        # a coupling at 0 still has a df/dJ (and contributes 0 Hz); two pairs under one label are summed
        j0 = np.array([132.0, -4.5, 0.0])
        names = {(0, 1): "A", (0, 2): "B", (1, 2): "C"}
        lines = [line for line in line_table(system(j0), names, min_relative_amplitude=0.05) if line.frequency_hz > 50]
        self.assertTrue(all("C" in line.df_dj and line.contribution_hz["C"] == 0.0 for line in lines))
        h = 1e-4
        fp = np.asarray(compute_transitions(system(j0 + [0, 0, h])).frequencies_hz)
        fm = np.asarray(compute_transitions(system(j0 - [0, 0, h])).frequencies_hz)
        for line in lines:
            numeric = (fp[np.argmin(np.abs(fp - line.frequency_hz))] - fm[np.argmin(np.abs(fm - line.frequency_hz))]) / (2 * h)
            self.assertAlmostEqual(line.df_dj["C"], numeric, places=4)
        tied = line_table(system(np.array([132.0, -4.5, 6.8])), {(0, 1): "A", (0, 2): "T", (1, 2): "T"})
        single = line_table(system(np.array([132.0, -4.5, 6.8])), names)
        for a, b in zip(tied, single):
            self.assertAlmostEqual(a.df_dj["T"], b.df_dj["B"] + b.df_dj["C"], places=9)
            self.assertAlmostEqual(sum(a.contribution_hz.values()), a.frequency_hz, places=6)

    def test_second_order_and_cross_terms(self):
        j0 = np.array([132.0, -4.5, 6.8])
        names = {(0, 1): "A", (0, 2): "B", (1, 2): "C"}
        lines = [line for line in line_table(system(j0), names, min_relative_amplitude=0.05, second_order=True)
                 if line.frequency_hz > 50]
        for line in lines:
            hj = [sum(line.hessian[(a, b)] * v for b, v in zip("ABC", j0)) for a in "ABC"]
            self.assertLess(max(abs(v) for v in hj), 1e-4)            # homogeneity: H J = 0
            self.assertGreater(line.trust_step_hz, 0.0)
        # independent check: second differences of the transition list for a cross term (B, C)
        h = 1e-2

        def freqs(db, dc):
            return np.asarray(compute_transitions(system(j0 + np.array([0.0, db, dc]))).frequencies_hz)
        for line in lines:
            vals = []
            for sb, sc in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                fr = freqs(sb * h, sc * h)
                vals.append(fr[np.argmin(np.abs(fr - line.frequency_hz))])
            numeric = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * h * h)
            self.assertAlmostEqual(line.hessian[("B", "C")], numeric, places=2, msg=line.frequency_hz)
        # second order improves on first order for a 0.3 Hz change of B and C
        dj = {"B": 0.3, "C": -0.3}
        fr = freqs(0.3, -0.3)
        for line in lines:
            exact = fr[np.argmin(np.abs(fr - line.shifted(dj, second_order=True)))]
            self.assertLessEqual(abs(line.shifted(dj, second_order=True) - exact),
                                 abs(line.shifted(dj) - exact) + 1e-6)

    def test_first_order_prediction(self):
        j0 = np.array([132.0, -4.5, 6.8])
        lines = [line for line in line_table(system(j0), min_relative_amplitude=0.05) if line.frequency_hz > 50]
        moved = np.asarray(compute_transitions(system(j0 + np.array([0.0, 0.02, 0.0]))).frequencies_hz)
        pred = predict_lines(lines, {"J(0,2)": 0.02})
        for p in pred:
            self.assertLess(np.min(np.abs(moved - p)), 1e-4)


if __name__ == "__main__":
    unittest.main()
