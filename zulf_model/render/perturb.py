"""Nuisance processes for training spectra.

Signal nuisances (linewidth, gain, phase, delay) act on the analytic render.
Additive nuisances (white noise, slow drift, interference lines) are generated
in the time domain and pass through `process_record` like experimental data.
All randomness comes from an explicit `numpy.random.Generator`.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..physics.transitions import TransitionList
from .acquisition import Acquisition, evaluate_spectrum, process_record
from .renderer import Renderer


@dataclass(frozen=True)
class PerturbationConfig:
    rate_range_per_s: Tuple[float, float] = (0.3, 6.0)
    gaussian_probability: float = 0.5
    gaussian_sigma_range_hz: Tuple[float, float] = (0.02, 1.0)
    family_rate_spread: float = 0.3
    family_split_probability: float = 0.3
    global_phase_range_rad: Tuple[float, float] = (-np.pi, np.pi)
    phase_delay_range_s: Tuple[float, float] = (-0.002, 0.002)
    component_gain_log10_range: Tuple[float, float] = (-0.3, 0.3)
    component_phase_spread_rad: float = 0.2
    snr_range: Tuple[float, float] = (5.0, 500.0)
    drift_probability: float = 0.3
    drift_relative_amplitude: Tuple[float, float] = (0.1, 3.0)
    interference_probability: float = 0.2
    mains_hz: Tuple[float, ...] = (50.0, 60.0)
    interference_max_lines: int = 3
    interference_relative_amplitude: Tuple[float, float] = (0.05, 1.0)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "PerturbationConfig":
        data = {k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items() if k in cls.__dataclass_fields__}
        return cls(**data)

    @classmethod
    def load(cls, path) -> "PerturbationConfig":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


@dataclass
class ComponentRender:
    rates_per_s: np.ndarray
    family_edges_hz: np.ndarray
    gain: complex
    gaussian_sigma_hz: float = 0.0


@dataclass
class RenderParams:
    components: List[ComponentRender]
    global_phase_rad: float = 0.0
    phase_delay_s: float = 0.0
    snr: Optional[float] = None
    drift: Optional[dict] = None
    interference: List[dict] = field(default_factory=list)
    noise_seed: int = 0

    def to_dict(self) -> dict:
        return {"components": [{"rates_per_s": c.rates_per_s.tolist(), "family_edges_hz": c.family_edges_hz.tolist(),
                                "gain": [c.gain.real, c.gain.imag], "gaussian_sigma_hz": c.gaussian_sigma_hz}
                               for c in self.components],
                "global_phase_rad": self.global_phase_rad, "phase_delay_s": self.phase_delay_s,
                "snr": self.snr, "drift": self.drift, "interference": self.interference,
                "noise_seed": self.noise_seed}

    @classmethod
    def from_dict(cls, data: dict) -> "RenderParams":
        comps = [ComponentRender(np.asarray(c["rates_per_s"], float), np.asarray(c["family_edges_hz"], float),
                                 complex(*c["gain"]), float(c.get("gaussian_sigma_hz", 0.0)))
                 for c in data["components"]]
        return cls(comps, data["global_phase_rad"], data["phase_delay_s"], data["snr"], data["drift"],
                   list(data["interference"]), int(data["noise_seed"]))


def _log_uniform(rng, lo, hi, size=None):
    return np.exp(rng.uniform(np.log(lo), np.log(hi), size))


def sample_render_params(rng: np.random.Generator, transitions: Sequence[TransitionList],
                         config: PerturbationConfig, noiseless: bool = False) -> RenderParams:
    """Draw nuisance parameters for a list of component transition lists."""
    comps = []
    for tl in transitions:
        base = _log_uniform(rng, *config.rate_range_per_s)
        edges = np.zeros(0)
        if len(tl) > 1 and rng.random() < config.family_split_probability:
            lo, hi = tl.frequencies_hz.min(), tl.frequencies_hz.max()
            edges = np.array([rng.uniform(lo, hi)])
        rates = base * np.exp(rng.normal(0, config.family_rate_spread, len(edges) + 1))
        rates = np.clip(rates, *config.rate_range_per_s)
        gain = 10 ** rng.uniform(*config.component_gain_log10_range) * np.exp(
            1j * rng.normal(0, config.component_phase_spread_rad))
        sigma = 0.0
        if rng.random() < config.gaussian_probability:
            sigma = float(_log_uniform(rng, *config.gaussian_sigma_range_hz))
        comps.append(ComponentRender(rates, edges, complex(gain), sigma))
    params = RenderParams(comps, float(rng.uniform(*config.global_phase_range_rad)),
                          float(rng.uniform(*config.phase_delay_range_s)), noise_seed=int(rng.integers(2**31)))
    if noiseless:
        return params
    params.snr = float(_log_uniform(rng, *config.snr_range))
    if rng.random() < config.drift_probability:
        params.drift = {"relative_amplitude": float(rng.uniform(*config.drift_relative_amplitude)),
                        "order": int(rng.integers(1, 4)), "seed": int(rng.integers(2**31))}
    if rng.random() < config.interference_probability:
        base = float(rng.choice(config.mains_hz))
        for _ in range(int(rng.integers(1, config.interference_max_lines + 1))):
            params.interference.append({"frequency_hz": base * int(rng.integers(1, 9)),
                                        "relative_amplitude": float(rng.uniform(*config.interference_relative_amplitude)),
                                        "phase_rad": float(rng.uniform(-np.pi, np.pi))})
    return params


def render_signal(renderer: Renderer, transitions: Sequence[TransitionList], params: RenderParams,
                  frequencies_hz: np.ndarray, contributions: Sequence[float]) -> np.ndarray:
    """Noiseless processed spectrum of a mixture."""
    total = np.zeros(len(frequencies_hz), complex)
    phase = np.exp(1j * params.global_phase_rad)
    for tl, comp, weight in zip(transitions, params.components, contributions):
        labelled = tl.split_families(comp.family_edges_hz)
        total += renderer.render(labelled, comp.rates_per_s, frequencies_hz, comp.gain * phase * weight,
                                 params.phase_delay_s, comp.gaussian_sigma_hz)
    return total


def additive_time_domain(acquisition: Acquisition, params: RenderParams, signal_peak: float,
                         spectral_noise_std: float) -> np.ndarray:
    """White noise, drift and interference as a real FID over the full record."""
    rng = np.random.default_rng(params.noise_seed)
    n_full = acquisition.points
    # White noise: spectral std of (1/n) sum of n iid samples is sigma / sqrt(n).
    sigma = spectral_noise_std * np.sqrt(acquisition.n)
    fid = rng.normal(0, sigma, n_full) if spectral_noise_std > 0 else np.zeros(n_full)
    t = np.linspace(-1, 1, n_full)
    if params.drift:
        drng = np.random.default_rng(params.drift["seed"])
        coeffs = drng.normal(size=params.drift["order"] + 1)
        poly = np.polynomial.legendre.legval(t, coeffs)
        poly /= max(np.abs(poly).max(), 1e-30)
        fid += poly * params.drift["relative_amplitude"] * sigma * 10
    times = acquisition.times()
    for line in params.interference:
        amp = line["relative_amplitude"] * signal_peak * 2
        fid += amp * np.cos(2 * np.pi * line["frequency_hz"] * times + line["phase_rad"])
    return fid


def mixture_fid(renderer: Renderer, transitions: Sequence[TransitionList], params: RenderParams,
                contributions: Sequence[float]) -> np.ndarray:
    """Noiseless real FID of the whole mixture over the full record."""
    total = np.zeros(renderer.acq.points)
    phase = np.exp(1j * params.global_phase_rad)
    for tl, comp, weight in zip(transitions, params.components, contributions):
        labelled = tl.split_families(comp.family_edges_hz)
        total += renderer.synthesize(labelled, comp.rates_per_s, comp.gain * phase * weight,
                                     params.phase_delay_s, gaussian_sigma_hz=comp.gaussian_sigma_hz)
    return total


def render_observation(renderer: Renderer, transitions: Sequence[TransitionList], params: RenderParams,
                       frequencies_hz: np.ndarray, contributions: Sequence[float]) -> Tuple[np.ndarray, np.ndarray]:
    """Return (noisy spectrum, noiseless spectrum) on `frequencies_hz`.

    Signal and nuisances are summed as one real FID, which then passes once
    through the acquisition operator, exactly as experimental data would.
    """
    acq = renderer.acq
    timer = renderer.timer
    with timer.section("observation.synthesize"):
        clean_fid = mixture_fid(renderer, transitions, params, contributions)
    with timer.section("observation.process"):
        clean = evaluate_spectrum(process_record(clean_fid, acq), acq, frequencies_hz)
    if params.snr is None:
        return clean.copy(), clean
    peak = float(np.abs(clean).max()) if len(clean) else 0.0
    noise_std = peak / params.snr if peak > 0 else 1.0
    with timer.section("observation.nuisance"):
        fid = clean_fid + additive_time_domain(acq, params, peak, noise_std)
    with timer.section("observation.process"):
        noisy = evaluate_spectrum(process_record(fid, acq), acq, frequencies_hz)
    return noisy, clean
