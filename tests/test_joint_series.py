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
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 1.0])
        z = joint.pack([p.vector() for p in joint.params]) + 0.05
        joint.set_priors([0], [z[0] + 1.0], [0.5], 2.0)
        # couplings are exactly linear in x (independent check of spectrum_vector)
        i = joint.coupling.index(model.coupling_names["J(Ca,Ha)"][0])
        vals = [joint.spectrum_vector(z, s)[joint.col[joint.coupling[i]]] for s in range(3)]
        self.assertAlmostEqual(vals[2] - vals[1], (vals[1] - vals[0]) * 0.5 / 0.4, places=10)
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))


if __name__ == "__main__":
    unittest.main()
