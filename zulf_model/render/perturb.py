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
    # "jitter": weight = abundance-based contribution times 10**U(component_gain_log10_range);
    # "free": weight = 10**U(free_weight_log10_range) independent of abundance (relaxation, polarization
    # transfer and detection can change isotopologue ratios arbitrarily).
    component_ratio_mode: str = "jitter"
    free_weight_log10_range: Tuple[float, float] = (-1.5, 0.0)
    component_phase_spread_rad: float = 0.2
    snr_range: Tuple[float, float] = (5.0, 500.0)
    drift_probability: float = 0.3
    drift_relative_amplitude: Tuple[float, float] = (0.1, 3.0)
    interference_probability: float = 0.2
    mains_hz: Tuple[float, ...] = (50.0, 60.0)
    interference_max_lines: int = 3
    interference_relative_amplitude: Tuple[float, float] = (0.05, 1.0)
    # Instrument response seen in experiments: slow multi-exponential baseline,
    # fast ringing after the pulse and ADC saturation at the start of the record.
    # Amplitudes are relative to the peak of the noiseless molecular FID.
    instrument_probability: float = 0.0
    baseline_exponentials: Tuple[int, int] = (1, 3)
    baseline_rate_per_s: Tuple[float, float] = (0.3, 30.0)
    baseline_relative_amplitude: Tuple[float, float] = (10.0, 1e4)
    ringing_probability: float = 0.7
    ringing_frequency_hz: Tuple[float, float] = (200.0, 1500.0)
    ringing_rate_per_s: Tuple[float, float] = (80.0, 500.0)
    ringing_relative_amplitude: Tuple[float, float] = (10.0, 1e4)
    saturation_probability: float = 0.5
    saturation_duration_s: Tuple[float, float] = (0.001, 0.005)

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
    instrument: Optional[dict] = None

    def to_dict(self) -> dict:
        return {"components": [{"rates_per_s": c.rates_per_s.tolist(), "family_edges_hz": c.family_edges_hz.tolist(),
                                "gain": [c.gain.real, c.gain.imag], "gaussian_sigma_hz": c.gaussian_sigma_hz}
                               for c in self.components],
                "global_phase_rad": self.global_phase_rad, "phase_delay_s": self.phase_delay_s,
                "snr": self.snr, "drift": self.drift, "interference": self.interference,
                "noise_seed": self.noise_seed, "instrument": self.instrument}

    @classmethod
    def from_dict(cls, data: dict) -> "RenderParams":
        comps = [ComponentRender(np.asarray(c["rates_per_s"], float), np.asarray(c["family_edges_hz"], float),
                                 complex(*c["gain"]), float(c.get("gaussian_sigma_hz", 0.0)))
                 for c in data["components"]]
        return cls(comps, data["global_phase_rad"], data["phase_delay_s"], data["snr"], data["drift"],
                   list(data["interference"]), int(data["noise_seed"]), data.get("instrument"))


def _log_uniform(rng, lo, hi, size=None):
    return np.exp(rng.uniform(np.log(lo), np.log(hi), size))


def sample_render_params(rng: np.random.Generator, transitions: Sequence[TransitionList],
                         config: PerturbationConfig, noiseless: bool = False,
                         contributions: Optional[Sequence[float]] = None) -> RenderParams:
    """Draw nuisance parameters for a list of component transition lists.

    The rendered weight of component c is contributions[c] * |gain_c|. In "free" ratio mode the gain
    magnitude is chosen so that the weight is 10**U(free_weight_log10_range) regardless of abundance.
    """
    if config.component_ratio_mode not in ("jitter", "free"):
        raise ValueError("component_ratio_mode must be 'jitter' or 'free'.")
    comps = []
    for index, tl in enumerate(transitions):
        base = _log_uniform(rng, *config.rate_range_per_s)
        edges = np.zeros(0)
        if len(tl) > 1 and rng.random() < config.family_split_probability:
            lo, hi = tl.frequencies_hz.min(), tl.frequencies_hz.max()
            edges = np.array([rng.uniform(lo, hi)])
        rates = base * np.exp(rng.normal(0, config.family_rate_spread, len(edges) + 1))
        rates = np.clip(rates, *config.rate_range_per_s)
        if config.component_ratio_mode == "free":
            nominal = float(contributions[index]) if contributions is not None else 1.0
            magnitude = 10 ** rng.uniform(*config.free_weight_log10_range) / max(nominal, 1e-300)
        else:
            magnitude = 10 ** rng.uniform(*config.component_gain_log10_range)
        gain = magnitude * np.exp(1j * rng.normal(0, config.component_phase_spread_rad))
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
    if rng.random() < config.instrument_probability:
        count = int(rng.integers(config.baseline_exponentials[0], config.baseline_exponentials[1] + 1))
        instrument = {"baseline": [{"rate_per_s": float(_log_uniform(rng, *config.baseline_rate_per_s)),
                                    "relative_amplitude": float(_log_uniform(rng, *config.baseline_relative_amplitude))
                                    * float(rng.choice([-1.0, 1.0]))} for _ in range(count)]}
        if rng.random() < config.ringing_probability:
            instrument["ringing"] = {"frequency_hz": float(rng.uniform(*config.ringing_frequency_hz)),
                                     "rate_per_s": float(_log_uniform(rng, *config.ringing_rate_per_s)),
                                     "relative_amplitude": float(_log_uniform(rng, *config.ringing_relative_amplitude)),
                                     "phase_rad": float(rng.uniform(-np.pi, np.pi))}
        if rng.random() < config.saturation_probability:
            instrument["saturation_s"] = float(rng.uniform(*config.saturation_duration_s))
        params.instrument = instrument
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


def instrument_fid(acquisition: Acquisition, instrument: dict, signal_peak: float) -> np.ndarray:
    """Additive instrument response (baseline exponentials and ringing) over the full record."""
    t = acquisition.times()
    out = np.zeros(acquisition.points)
    for term in instrument.get("baseline", []):
        out += term["relative_amplitude"] * signal_peak * np.exp(-term["rate_per_s"] * t)
    ring = instrument.get("ringing")
    if ring:
        out += (ring["relative_amplitude"] * signal_peak * np.exp(-ring["rate_per_s"] * t)
                * np.cos(2 * np.pi * ring["frequency_hz"] * t + ring["phase_rad"]))
    return out


def apply_saturation(fid: np.ndarray, acquisition: Acquisition, duration_s: float) -> np.ndarray:
    """Hold the first `duration_s` of the record at the extreme value reached there (ADC/amplifier limit)."""
    count = int(round(duration_s * acquisition.sampling_rate_hz))
    if count <= 0:
        return fid
    out = fid.copy()
    head = out[:count]
    out[:count] = head[np.argmax(np.abs(head))]
    return out


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
        if params.instrument:
            fid = fid + instrument_fid(acq, params.instrument, float(np.abs(clean_fid).max()))
            if params.instrument.get("saturation_s"):
                fid = apply_saturation(fid, acq, params.instrument["saturation_s"])
    with timer.section("observation.process"):
        noisy = evaluate_spectrum(process_record(fid, acq), acq, frequencies_hz)
    return noisy, clean
