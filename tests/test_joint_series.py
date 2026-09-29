import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from fit_joint_series import JointSeries                                # noqa: E402

from zulf_core.render.acquisition import Acquisition
from zulf_core.solver import ObservedSpectrum, RefineSettings
from zulf_core.solver.forward import MixtureForward
from zulf_core.solver.refine import _signal_kwargs
from zulf_hypothesis import build_model, template
from zulf_hypothesis.search import _settings_for


class JointSeriesTests(unittest.TestCase):
    def test_jacobian_matches_finite_differences_and_linear_couplings(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        rng = np.random.default_rng(2)
        obs = []
        for shift in (0.0, 0.4, 0.9):
            p = settings.parameterize(model.interpretation)
            name = model.coupling_names["J(Ca,Ha)"][0]
            p.set(name, p.parameters[name].value + shift)
            sig = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                                 **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model
            obs.append(ObservedSpectrum.from_spectrum(f, np.asarray(sig) + 0.01 * rng.normal(size=len(f)),
                                                      [(100.0, 200.0)], record=acq))
        for shape in ("linear", "monotone"):
            joint = JointSeries(model, settings, obs, [0.1, 0.5, 1.0], shape=shape,
                                signs=[1.0 if k % 2 else -1.0 for k in range(len(model.coupling_names))])
            table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params])
            z = joint.pack(table + np.arange(3)[:, None] * 0.2, [p.vector() for p in joint.params]) + 0.05
            joint.set_priors([0], [1.0 + joint.coupling_values(z, 0).mean()], [0.5], 2.0)
            i = joint.coupling.index(model.coupling_names["J(Ca,Ha)"][0])
            vals = np.array([joint.spectrum_vector(z, s)[joint.col[joint.coupling[i]]] for s in range(3)])
            steps = np.diff(vals)
            if shape == "linear":      # exactly linear in x
                self.assertAlmostEqual(steps[1], steps[0] * 0.5 / 0.4, places=10)
            else:                      # one direction only (independent check of the design matrix)
                self.assertTrue(np.all(steps >= -1e-12) or np.all(steps <= 1e-12))
            jac = joint.jacobian(z)
            h = 1e-6
            numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                       for e in np.eye(len(z))])
            self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0), msg=shape)

if __name__ == "__main__":
    unittest.main()
