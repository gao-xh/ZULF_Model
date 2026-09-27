"""One dataset from raw FID to processed spectra, with a record of every choice."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from zulf_core.render.acquisition import evaluate_spectrum, process_record
from zulf_core.render.phasing import phase_correct
from zulf_core.solver import ObservedSpectrum

from .diagnostics import RawDiagnostics, diagnose_raw
from .phase import PhaseResult, calibrated_phase, phase_dataset
from .plan import ProcessingPlan, plan_for_dataset


@dataclass
class ProcessedDataset:
    sample_id: str
    diagnostics: RawDiagnostics
    plan: ProcessingPlan
    observed: ObservedSpectrum                 # complex, fitted ranges (what the solver uses)
    frequencies_hz: np.ndarray = field(repr=False)   # continuous grid over the fitted span
    spectrum: np.ndarray = field(repr=False)          # complex processed spectrum on that grid
    phase: Optional[PhaseResult] = None
    phased: Optional[ObservedSpectrum] = None  # real (absorption) part with the phase recorded for the solver

    def record(self) -> dict:
        return {"sample_id": self.sample_id, "diagnostics": self.diagnostics.to_dict(), "plan": self.plan.to_dict(),
                "phase": self.phase.to_dict() if self.phase else None}

    def save_record(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(self.record(), fh, indent=1, default=float)

    def figure(self, path: str) -> Optional[str]:
        """Raw FID start (edge, ringing, crop), magnitude overview, phased real part over the fitted span."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            return None
        fs = self.plan.sampling_rate_hz
        fig, ax = plt.subplots(3, 1, figsize=(12, 10))
        n0 = min(int(0.1 * fs), self.plan.points)
        ax[0].plot(np.arange(n0) / fs * 1e3, self._fid[:n0], lw=0.8, color="#333333")
        if self.diagnostics.edge_time_s:
            ax[0].axvline(self.diagnostics.edge_time_s * 1e3, color="#d1495b", ls="--", lw=0.9, label="switching edge")
        ax[0].axvline(self.plan.start_sample / fs * 1e3, color="#00798c", ls=":", lw=0.9, label="crop start")
        ax[0].set_title("raw FID, first 100 ms")
        ax[0].set_xlabel("time (ms)")
        ax[0].legend(frameon=False)
        f, v = self.frequencies_hz, self.spectrum
        noise = 1.4826 * float(np.median(np.abs(v.real - np.median(v.real)))) or 1.0
        ax[1].plot(f, np.abs(v) / noise, lw=0.7, color="#00798c")
        ax[1].set_title("magnitude / noise")
        if self.phase is not None:
            real = phase_correct(v, f, self.phase.phase0_rad, self.phase.delay_s, self.plan.acquisition()).real
            ax[2].plot(f, real / noise, lw=0.7, color="#333333")
            ax[2].axhline(0, color="#bbbbbb", lw=0.6)
            ax[2].set_title(f"phased real part ({self.phase.criterion}, {self.phase.delay_mode}): phase0 "
                            f"{np.degrees(self.phase.phase0_rad):.1f} deg, delay {self.phase.delay_s * 1e3:.3f} ms, "
                            f"edge {self.phase.edge_delay_s * 1e3 if self.phase.edge_delay_s else float('nan'):.3f} ms")
        for a in ax[1:]:
            for line in self.plan.instrument_lines_hz:
                if f[0] <= line <= f[-1]:
                    a.axvline(line, color="#dddddd", ls=":", lw=0.8, zorder=0)
            a.set_xlabel("frequency (Hz)")
        fig.suptitle(f"{self.sample_id}: crop {self.plan.start_sample}-{self.plan.stop_sample}, SG "
                     f"{self.plan.sg_window}/{self.plan.sg_order}, apodization {self.plan.apodization_rate_per_s} 1/s",
                     x=0.01, ha="left", fontweight="bold")
        fig.tight_layout()
        fig.savefig(path, dpi=105)
        plt.close(fig)
        return path


def process_dataset(fid: np.ndarray, sampling_rate_hz: float = 4000.0, sample_id: str = "", plan=None,
                    defaults: Optional[dict] = None, phase_criterion: Optional[str] = "entropy",
                    delay_mode: str = "edge_prior", phase_calibration: Optional[dict] = None,
                    **phase_options) -> ProcessedDataset:
    """Raw FID -> per-dataset diagnostics -> plan (or the given plan) -> complex observation on the fitted ranges ->
    phase by global search and fine-tune on the continuous span of the ranges -> phased real observation.
    `phase_criterion=None` skips phasing (complex fits need no phase correction); "calibration" uses the instrument
    phase calibration with this dataset's switching edge (calibrated_phase); other names are model-free criteria
    searched globally and fine-tuned (phase_dataset)."""
    x = np.asarray(fid, float)
    diagnostics = diagnose_raw(x, sampling_rate_hz)
    plan = plan or plan_for_dataset(len(x), sampling_rate_hz, diagnostics, defaults)
    acq = plan.acquisition()
    observed = ObservedSpectrum.from_fid(x, acq, [tuple(r) for r in plan.ranges], zero_fill=plan.zero_fill,
                                         label=sample_id)
    step = sampling_rate_hz / acq.n / plan.zero_fill
    lo = min(r[0] for r in plan.ranges)
    hi = max(r[1] for r in plan.ranges)
    f = np.arange(int(np.ceil(lo / step)), int(np.floor(hi / step)) + 1) * step
    spectrum = evaluate_spectrum(process_record(x, acq), acq, f)
    out = ProcessedDataset(sample_id, diagnostics, plan, observed, f, spectrum)
    out._fid = x
    if phase_criterion:
        edge = -(diagnostics.edge_time_s + acq.time_origin_s) if diagnostics.edge_time_s else None
        if phase_criterion == "calibration":
            calibration = phase_calibration or (defaults or {}).get("phase_calibration")
            if not calibration:
                raise ValueError("phase_criterion 'calibration' needs phase_calibration (or defaults['phase_calibration']).")
            out.phase = calibrated_phase(edge, calibration)
        else:
            out.phase = phase_dataset(spectrum, f, acq, edge_delay_s=edge, criterion=phase_criterion,
                                      delay_mode=delay_mode, exclude_hz=plan.instrument_lines_hz, **phase_options)
        sel = observed.selected
        real = phase_correct(observed.values[sel], observed.frequencies_hz[sel], out.phase.phase0_rad,
                             out.phase.delay_s, acq).real
        out.phased = ObservedSpectrum.from_spectrum(observed.frequencies_hz[sel], real,
                                                    [tuple(r) for r in plan.ranges], record=acq,
                                                    phasing={"phase0_rad": out.phase.phase0_rad,
                                                             "delay_s": out.phase.delay_s}, label=sample_id)
        out.phased.metadata["zero_fill"] = plan.zero_fill
    return out
