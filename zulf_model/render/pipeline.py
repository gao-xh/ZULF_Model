"""From a ground-truth sample to a model-ready spectrum.

`ProcessingConfig.mode` selects the route:

* `pure`: sampling only, no processing (default);
* `fixed`: one explicit experimental recipe;
* `randomized`: recipe drawn per sample from configured ranges;
* `continuous`: infinite-record Lorentzian/Voigt lines, no sampling at all.

Noise, drift and interference come from `PerturbationConfig` and can all be
disabled (`noiseless=True`).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..physics.protocol import SUDDEN_DROP, Protocol
from ..physics.transitions import TransitionCache, TransitionList
from ..spec import ProblemSpec
from ..spinsystem import Interpretation
from ..timing import Timer
from .acquisition import Acquisition
from .features import spectrum_features
from .grid import SpectrumGrid
from .perturb import PerturbationConfig, RenderParams, render_observation, render_signal, sample_render_params
from .renderer import ContinuousRenderer, Renderer

MODES = ("pure", "fixed", "randomized", "continuous")


@dataclass(frozen=True)
class ProcessingConfig:
    mode: str = "pure"
    sampling_rate_hz: float = 1000.0
    points: int = 16384
    fixed: Dict[str, object] = field(default_factory=dict)
    start_s_range: Tuple[float, float] = (0.0, 0.2)
    sg_window_choices: Tuple[int, ...] = (0, 101, 301)
    sg_order_choices: Tuple[int, ...] = (2,)
    remove_mean_probability: float = 0.5
    time_origin_range_s: Tuple[float, float] = (0.0, 0.0)
    renderer_backend: str = "auto"  # auto | analytic (never builds an FID) | time (NUFFT FID)

    def __post_init__(self):
        if self.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}.")
        if self.renderer_backend not in Renderer.BACKENDS:
            raise ValueError(f"renderer_backend must be one of {Renderer.BACKENDS}.")

    @classmethod
    def from_dict(cls, data: dict) -> "ProcessingConfig":
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items()
                      if k in cls.__dataclass_fields__})

    @classmethod
    def load(cls, path) -> "ProcessingConfig":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def base(self) -> Acquisition:
        return Acquisition.pure(self.sampling_rate_hz, self.points)

    def acquisition(self, rng: Optional[np.random.Generator] = None) -> Acquisition:
        base = self.base()
        if self.mode in ("pure", "continuous"):
            return base
        if self.mode == "fixed":
            return base.with_processing(**self.fixed)
        rng = rng or np.random.default_rng()
        start = int(round(rng.uniform(*self.start_s_range) * self.sampling_rate_hz))
        window = int(rng.choice(self.sg_window_choices))
        order = int(rng.choice(self.sg_order_choices))
        if window and window <= order:
            window = 0
        return base.with_processing(start_sample=start, sg_window=window, sg_order=order,
                                    remove_mean=bool(rng.random() < self.remove_mean_probability),
                                    time_origin_s=float(rng.uniform(*self.time_origin_range_s)))


@dataclass
class RenderedSample:
    features: np.ndarray
    scale: float
    spectrum: np.ndarray
    clean: np.ndarray
    acquisition: dict
    params: dict
    interpretation: Interpretation
    timings_ms: dict = field(default_factory=dict)


class SampleRenderer:
    """Renders interpretations to model inputs on a fixed grid."""

    def __init__(self, spec: ProblemSpec, processing: ProcessingConfig = ProcessingConfig(),
                 perturbation: PerturbationConfig = PerturbationConfig(), protocol: Protocol = SUDDEN_DROP,
                 noiseless: bool = False, cache_entries: int = 8192, timer: Optional[Timer] = None):
        self.spec = spec
        self.processing = processing
        self.perturbation = perturbation
        self.protocol = protocol
        self.noiseless = noiseless
        if (processing.renderer_backend == "analytic" and processing.mode != "continuous"
                and perturbation.gaussian_probability > 0):
            raise ValueError("renderer_backend 'analytic' cannot render finite-record Gaussian broadening; set "
                             "perturbation.gaussian_probability to 0, use backend 'auto'/'time', or mode 'continuous'.")
        self.grid = SpectrumGrid.from_spec(spec.grid, processing.base())
        self.transitions = TransitionCache(cache_entries)
        self._renderers: Dict[Acquisition, Renderer] = {}
        self.timer = timer or Timer()

    def renderer(self, acquisition: Acquisition) -> Renderer:
        if acquisition not in self._renderers:
            if len(self._renderers) > 64:
                self._renderers.clear()
            self._renderers[acquisition] = Renderer(acquisition, timer=self.timer,
                                                    backend=self.processing.renderer_backend)
        return self._renderers[acquisition]

    def component_transitions(self, interpretation: Interpretation) -> List[TransitionList]:
        return [self.transitions.get(c.system, self.protocol) for c in interpretation.components]

    def render(self, interpretation: Interpretation, rng: np.random.Generator) -> RenderedSample:
        import time
        stamps = {}
        t0 = time.perf_counter()
        with self.timer.section("sample.transitions"):
            transitions = self.component_transitions(interpretation)
        stamps["transitions"] = time.perf_counter() - t0
        contributions = [c.contribution for c in interpretation.components]
        params = sample_render_params(rng, transitions, self.perturbation, noiseless=self.noiseless)
        f = self.grid.frequencies_hz
        if self.processing.mode == "continuous":
            acquisition = self.processing.base()
            record = acquisition.points / acquisition.sampling_rate_hz
            clean = self._continuous(transitions, params, contributions, record)
            spectrum = clean + self._spectral_noise(params, clean, len(f))
        elif self.processing.renderer_backend == "analytic":
            acquisition = self.processing.acquisition(rng)
            with self.timer.section("sample.observation"):
                renderer = self.renderer(acquisition)
                clean = render_signal(renderer, transitions, params, f, contributions)
                spectrum = clean + self._spectral_noise(params, clean, len(f))
        else:
            acquisition = self.processing.acquisition(rng)
            with self.timer.section("sample.observation"):
                spectrum, clean = render_observation(self.renderer(acquisition), transitions, params, f, contributions)
        stamps["render"] = time.perf_counter() - t0 - stamps["transitions"]
        features, scale = spectrum_features(spectrum, self.spec.grid.channels)
        stamps["total"] = time.perf_counter() - t0
        self.timer.add("sample.total", stamps["total"])
        return RenderedSample(features, scale, spectrum, clean, acquisition.to_dict(), params.to_dict(),
                              interpretation, {k: 1000 * v for k, v in stamps.items()})

    @staticmethod
    def _spectral_noise(params: RenderParams, clean: np.ndarray, size: int) -> np.ndarray:
        """White complex noise added directly in the frequency domain (no FID).

        Used by the FID-free routes; drift and interference lines are not modelled
        there. The time-domain route models them through the acquisition operator.
        """
        if params.snr is None:
            return np.zeros(size, complex)
        peak = float(np.abs(clean).max()) or 1.0
        rng = np.random.default_rng(params.noise_seed)
        std = peak / params.snr / np.sqrt(2)
        return rng.normal(0, std, size) + 1j * rng.normal(0, std, size)

    def _continuous(self, transitions, params: RenderParams, contributions, record_s: float) -> np.ndarray:
        renderer = ContinuousRenderer(record_s)
        total = np.zeros(len(self.grid), complex)
        phase = np.exp(1j * params.global_phase_rad)
        for tl, comp, weight in zip(transitions, params.components, contributions):
            labelled = tl.split_families(comp.family_edges_hz)
            total += renderer.render(labelled, comp.rates_per_s, self.grid.frequencies_hz,
                                     comp.gain * phase * weight, params.phase_delay_s, comp.gaussian_sigma_hz)
        return total
