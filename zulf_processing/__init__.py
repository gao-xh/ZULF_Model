"""zulf_processing: experimental ZULF data processing, one dataset at a time.

A separate package (like zulf_hypothesis): it turns a raw averaged FID into processed spectra and never
interprets them. It uses the single processing operator of zulf_core (`Acquisition`, `process_record`,
`evaluate_spectrum`; never re-implemented) and the phase convention of `zulf_core.render.phasing`.

Layers:
* `raw`         - reading FIDs and instrument settings (NMRduino .ini);
* `diagnostics` - raw-FID diagnostics per dataset: switching edge, plateau, ringing end, noise
                  (wraps zulf_core.diagnostics.diagnose_fid);
* `plan`        - processing parameters chosen per dataset (crop start and stop, SG window, apodization,
                  zero filling, ranges): `ProcessingPlan`, `plan_for_dataset`;
* `phase`       - per-dataset phase: delay from the switching edge, then a global search over
                  (phase0, delay) with a registered criterion and a local fine-tune (`phase_dataset`);
* `dataset`     - `process_dataset`: raw FID -> plan -> complex ObservedSpectrum -> phase -> phased real
                  spectrum, with an overview figure and a JSON record of every choice.
"""
from .diagnostics import RawDiagnostics, diagnose_raw, switching_edge
from .plan import ProcessingPlan, plan_for_dataset
from .phase import PHASE_CRITERIA, PhaseResult, phase_dataset, register_phase_criterion
from .dataset import ProcessedDataset, process_dataset
from .raw import load_fid, read_settings

__all__ = ["RawDiagnostics", "diagnose_raw", "switching_edge", "ProcessingPlan", "plan_for_dataset",
           "PHASE_CRITERIA", "PhaseResult", "phase_dataset", "register_phase_criterion", "ProcessedDataset",
           "process_dataset", "load_fid", "read_settings"]
