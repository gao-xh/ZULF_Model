import unittest

import numpy as np

from zulf_core.render.acquisition import Acquisition
from zulf_core.solver import ObservedSpectrum, RefineSettings, refine
from zulf_core.solver.forward import MixtureForward
from zulf_core.solver.refine import _signal_kwargs
from zulf_hypothesis import build_model, template
from zulf_hypothesis.scoring import free_parameter_count
from zulf_hypothesis.search import Evaluated, _settings_for
from zulf_hypothesis.uncertainty import coupling_uncertainties, uncertainty_markdown


class CouplingUncertaintyTests(unittest.TestCase):
    def test_standard_errors_match_monte_carlo_scatter(self):
        # Independent reference: the scatter of refits over noise realisations of one synthetic spectrum
        # (white complex noise at native spacing) against the linearised standard errors of one fit.
        model = build_model(template("CH-CH3"))
        base = RefineSettings(starts=1, band_weighting="none", background_order=-1,
                              continuation_rates_per_s=(0.0,))
        settings = _settings_for(model, "fixed", base)
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        param = settings.parameterize(model.interpretation)
        truth = MixtureForward(param, clean, gain_model=settings.gain_model, background=-1, band_weighting="none",
                               **_signal_kwargs(settings)).predict(
            param.vector(), fixed_gains=np.asarray(settings.amplitude_map, float)[:, 0].astype(complex))
        signal = np.asarray(truth.model)
        sigma = 0.02 * np.abs(signal).max()
        rng = np.random.default_rng(3)
        fits, predicted = [], None
        for trial in range(16):
            noise = sigma * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f))) / np.sqrt(2)
            obs = ObservedSpectrum.from_spectrum(f, signal + noise, [(100.0, 200.0)], record=acq)
            result = refine(model.interpretation, obs, settings)
            e = Evaluated("CH-CH3", "fixed", "test", model, result.summary(), result.prediction,
                          k=free_parameter_count(model, settings))
            e.couplings = model.named_couplings(result.parameters)
            if predicted is None:
                predicted = coupling_uncertainties(e, obs, base)
            fits.append([e.couplings[k] for k in predicted["keys"]])
        scatter = np.std(np.array(fits), axis=0, ddof=1)
        for key, s, p in zip(predicted["keys"], scatter, predicted["std_hz"]):
            self.assertGreater(p / s, 0.5, msg=f"{key}: predicted {p} vs scatter {s}")
            self.assertLess(p / s, 2.0, msg=f"{key}: predicted {p} vs scatter {s}")
        self.assertIn("J(Ca,Ha)", uncertainty_markdown(predicted))


class CouplingPriorTests(unittest.TestCase):
    def test_prior_shift_matches_linearised_closed_form(self):
        # Gaussian prior on one small coupling, centred 1.5 Hz off the data optimum. Independent reference: the
        # linearised posterior mode x = x_d - (H + P)^-1 P (x_d - mu), H = J^T J at the prior-free optimum x_d,
        # P = diag(weight / (sigma norm)^2) on the prior parameter.
        import dataclasses
        model = build_model(template("CH-CH3"))
        base = RefineSettings(starts=1, band_weighting="none", background_order=-1,
                              continuation_rates_per_s=(0.0,))
        settings = _settings_for(model, "fixed", base)
        acq = Acquisition.pure(1000.0, 4000)
        f = np.arange(100.0, 200.0, 0.25)
        clean = ObservedSpectrum.from_spectrum(f, np.zeros(len(f), complex), [(100.0, 200.0)], record=acq)
        param = settings.parameterize(model.interpretation)
        kw = _signal_kwargs(settings)
        signal = np.asarray(MixtureForward(param, clean, gain_model=settings.gain_model, background=-1,
                                           band_weighting="none", **kw).predict(
            param.vector(), fixed_gains=np.asarray(settings.amplitude_map, float)[:, 0].astype(complex)).model)
        noise = 0.02 * np.abs(signal).max() * (np.array([1.0, 1j]) @ np.random.default_rng(5).normal(size=(2, len(f))))
        obs = ObservedSpectrum.from_spectrum(f, signal + noise, [(100.0, 200.0)], record=acq)
        free = refine(model.interpretation, obs, settings)
        name = model.coupling_names["J(Ca,Hb)"][0]
        names = settings.parameterize(model.interpretation).free_names
        i = names.index(name)
        x_d = np.array([free.parameters[n] for n in names])
        fw = MixtureForward(settings.parameterize(model.interpretation), obs, gain_model=settings.gain_model,
                            background=-1, band_weighting="none", **kw)
        jac = fw.jacobian(x_d)
        mu, sigma, weight = x_d[i] + 1.5, 0.5, 1.0          # moderate prior: shift of a few 0.1 Hz (linear regime)
        p = np.zeros((len(x_d), len(x_d)))
        p[i, i] = weight / (sigma * fw.norm) ** 2
        m = x_d.copy()
        m[i] = mu
        expected = x_d - np.linalg.solve(jac.T @ jac + p, p @ (x_d - m))
        with_prior = refine(model.interpretation, obs, dataclasses.replace(
            settings, priors=((name, mu, sigma),), prior_weight=weight))
        got = np.array([with_prior.parameters[n] for n in names])
        shift = got[i] - x_d[i]
        self.assertGreater(abs(shift), 0.05)                              # the prior acts
        self.assertAlmostEqual(got[i], expected[i], delta=0.1 * abs(expected[i] - x_d[i]))

        def objective(x):                                                 # data chi2 plus the prior row
            extra = (x[i] - mu) * np.sqrt(weight) / sigma / fw.norm
            return fw.predict(x).score + extra ** 2
        self.assertLessEqual(objective(got), objective(expected) + 1e-9)


if __name__ == "__main__":
    unittest.main()
