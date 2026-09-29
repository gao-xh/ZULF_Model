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
    def test_monotone_couplings_and_jacobian_against_finite_differences(self):
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
        table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params])
        z = joint.pack(table + np.array([0.0, 0.3, 0.2])[:, None], [p.vector() for p in joint.params])
        z[:joint.nt:joint.m] += 0.05
        z[1:joint.nt:joint.m] = np.where(np.arange(joint.nc) % 2, 0.7, -0.4)    # both directions
        joint.set_priors([0], [1.0 + joint.coupling_values(z, 0).mean()], [0.5], 2.0)
        for k in range(joint.nc):   # monotone by construction, direction = sign of A (independent check)
            steps = np.diff(joint.coupling_values(z, k))
            a = joint.theta(z, k)[1]
            self.assertTrue(np.all(steps * np.sign(a) >= -1e-12))
            self.assertAlmostEqual(joint.coupling_values(z, k)[-1] - joint.coupling_values(z, k)[0], a, places=10)
        for sigma in (0.0, 1.5):         # unsmoothed, and the coarse-to-fine smoothed objective
            joint.set_smoothing(sigma)
            jac = joint.jacobian(z)
            h = 1e-6
            numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                       for e in np.eye(len(z))])
            self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0), msg=str(sigma))
        # smoothing is the same linear kernel on data and model: a constant offset of the residual is kept
        joint.set_smoothing(1.5)
        n_res = len(joint.forwards[0].predict(joint.spectrum_vector(z, 0)).residual)   # 2F for complex data
        self.assertEqual(joint.smoothing[0].shape, (n_res, n_res))

    def test_spectra_at_one_concentration_share_their_couplings(self):
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        rng = np.random.default_rng(6)
        obs = [ObservedSpectrum.from_spectrum(f, rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)),
                                              [(100.0, 200.0)], record=acq) for _ in range(3)]
        joint = JointSeries(model, settings, obs, [0.1, 0.5, 0.5])
        self.assertEqual(joint.nn, 2)
        table = np.array([[p.vector()[joint.col[n]] for n in joint.coupling] for p in joint.params[:2]])
        z = joint.pack(table + np.array([0.0, 0.3])[:, None], [p.vector() for p in joint.params])
        a, b = joint.spectrum_vector(z, 1), joint.spectrum_vector(z, 2)
        for n in joint.coupling:
            self.assertEqual(a[joint.col[n]], b[joint.col[n]])
        jac = joint.jacobian(z)
        h = 1e-6
        numeric = np.column_stack([(joint.residual(z + h * e) - joint.residual(z - h * e)) / (2 * h)
                                   for e in np.eye(len(z))])
        self.assertLess(np.abs(jac - numeric).max(), 1e-5 * max(np.abs(numeric).max(), 1.0))

    def test_left_out_prediction_recovers_a_bracketed_spectrum(self):
        # three synthetic spectra, J(Ca,Ha) shifted 0 / 0.4 / 0.9 Hz: leave the middle out, fit the outer two with
        # their true couplings, predict the middle within the bracket; the prediction must reach the noise level
        import tempfile
        from fit_joint_series import predict_left_out
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        rng = np.random.default_rng(4)
        tmp = Path(tempfile.mkdtemp())
        name = model.coupling_names["J(Ca,Ha)"][0]
        obs, entries, noise = [], [], 0.002
        for s, shift in enumerate((0.0, 0.4, 0.9)):
            p = settings.parameterize(model.interpretation)
            p.set(name, p.parameters[name].value + shift)
            sig = np.asarray(MixtureForward(p, clean, gain_model=settings.gain_model, background=-1,
                                            band_weighting="none", **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model).real
            y = sig + noise * rng.normal(size=len(f))
            np.save(tmp / f"f{s}.npy", f)
            np.save(tmp / f"y{s}.npy", y)
            entries.append({"id": f"s{s}", "x": [0.1, 0.5, 1.0][s], "freq": str(tmp / f"f{s}.npy"),
                            "values": str(tmp / f"y{s}.npy")})
            obs.append(ObservedSpectrum.from_spectrum(f, y, [(100.0, 200.0)], record=None, real_only=True))
        joint = JointSeries(model, settings, [obs[0], obs[2]], [0.1, 1.0])
        truth = [settings.parameterize(model.interpretation).vector() for _ in range(2)]
        for x, shift in zip(truth, (0.0, 0.9)):
            x[joint.col[name]] += shift
        z = joint.pack(np.array([[x[joint.col[n]] for n in joint.coupling] for x in truth]), truth)
        key_of = {n: k for k, names in model.coupling_names.items() for n in names}
        out = predict_left_out(joint, z, entries[1], (100.0, 200.0), key_of, settings)
        # reference: the middle spectrum fitted at its true couplings (the noise level of this measure)
        p_mid = settings.parameterize(model.interpretation)
        x_mid = p_mid.vector()
        x_mid[joint.col[name]] += 0.4
        fw = MixtureForward(p_mid, obs[1], gain_model=settings.gain_model, background=-1, band_weighting="none",
                            **_signal_kwargs(settings))
        pred = fw.predict(x_mid)
        m = np.ones(len(fw.y), bool)
        noise_level = float(np.linalg.norm(fw.mismatch(pred.model, m)) / np.linalg.norm(fw.mismatch(0 * fw.y, m)))
        self.assertLess(out["relative_residual"], 1.05 * noise_level)
        u = out["u"]["J(Ca,Ha)"]
        self.assertAlmostEqual(u, 0.4 / 0.9, delta=0.1)               # true position inside the bracket


    def test_coordinate_scan_recovers_a_displaced_coupling(self):
        from fit_joint_series import _coordinate_scan
        model = build_model(template("CH-CH3"))
        settings = _settings_for(model, "ratios", RefineSettings(band_weighting="none", background_order=-1))
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        obs = []
        for shift in (0.0, 0.3):
            p = settings.parameterize(model.interpretation)
            p.set(model.coupling_names["J(Ca,Hb)"][0], p.parameters[model.coupling_names["J(Ca,Hb)"][0]].value + shift)
            sig = MixtureForward(p, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                                 **_signal_kwargs(settings)).predict(
                p.vector(), fixed_gains=np.ones(len(model.component_labels), complex)).model
            obs.append(ObservedSpectrum.from_spectrum(f, np.asarray(sig), [(100.0, 200.0)], record=acq))
        joint = JointSeries(model, settings, obs, [0.1, 1.0])
        truth = [p.vector() for p in joint.params]
        name = model.coupling_names["J(Ca,Hb)"][0]
        truth[1][joint.col[name]] += 0.3
        table = np.array([[x[joint.col[n]] for n in joint.coupling] for x in truth])
        z_true = joint.pack(table, truth)
        k = joint.coupling.index(name)
        z0 = z_true.copy()
        z0[k * joint.m] += 3.0                                   # far outside the local basin of a 1 Hz line
        lower, upper = joint.bounds(20.0)
        z = _coordinate_scan(joint, z0, lower, upper, 50, cycles=1, half_width=4.0, step=0.25)
        self.assertLess(abs(z[k * joint.m] - z_true[k * joint.m]), 0.3)


if __name__ == "__main__":
    unittest.main()
