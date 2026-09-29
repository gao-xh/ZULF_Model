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


if __name__ == "__main__":
    unittest.main()
