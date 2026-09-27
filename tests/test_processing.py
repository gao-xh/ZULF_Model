import sys
import unittest
from pathlib import Path

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_core.render import Acquisition, Renderer
from zulf_core.render.acquisition import evaluate_spectrum, process_record
from zulf_core.render.phasing import phase_correct

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_solver import methine_isotopologue, methyl_isotopologue  # noqa: E402


def synthetic_fid(phase0, delay, points=16384, fs=1000.0, start=100, noise=0.01, edge=None):
    acq = Acquisition(fs, points, start_sample=start, sg_window=301, sg_order=2)
    renderer = Renderer(acq)
    fid = (renderer.synthesize(compute_transitions(methine_isotopologue()), 1.0, gain=np.exp(1j * phase0),
                               phase_delay_s=delay)
           + renderer.synthesize(compute_transitions(methyl_isotopologue()), 1.5, gain=2 * np.exp(1j * phase0),
                                 phase_delay_s=delay)
           + np.random.default_rng(0).normal(0, noise, points))
    return acq, fid


class SwitchingEdgeTests(unittest.TestCase):
    def test_edge_time_between_samples(self):
        from zulf_processing import switching_edge
        t = np.arange(4096) / 4000.0
        edge = 0.00346
        x = np.where(t < edge - 0.0005, 28840.0, 0.0)
        ramp = (t >= edge - 0.0005) & (t < edge + 0.0005)          # linear 1 ms edge centred on `edge`
        x[ramp] = 28840.0 * (1 - (t[ramp] - (edge - 0.0005)) / 0.001)
        after = t >= edge + 0.0005
        x[after] = -8000.0 * np.exp(-(t[after] - edge) / 0.01)
        x += np.random.default_rng(0).normal(0, 5.0, len(x))
        out = switching_edge(x, 4000.0)
        # half height between the plateau and the first extremum lies on the ramp (independent closed form)
        extreme = -8000.0 * np.exp(-0.0005 / 0.01)
        expected = edge - 0.0005 + 0.001 * (28840.0 - 0.5 * (28840.0 + extreme)) / 28840.0
        self.assertAlmostEqual(out["edge_time_s"], expected, delta=2e-5)


class PlanTests(unittest.TestCase):
    def test_crop_after_ringing_and_reasons(self):
        from zulf_processing import ProcessingPlan, plan_for_dataset
        from zulf_processing.diagnostics import RawDiagnostics
        d = RawDiagnostics(4000.0, 65516, 0.0035, 0.00275, 0.05075, 0.5, True)
        plan = plan_for_dataset(65516, 4000.0, d)
        self.assertEqual(plan.start_sample, 203)                       # ceil(50.75 ms * 4 kHz)
        self.assertEqual(plan.stop_sample - plan.start_sample, 4000)    # default 1 s window kept
        self.assertIn("ringing", plan.reasons["start_sample"])
        early = plan_for_dataset(65516, 4000.0, RawDiagnostics(4000.0, 65516, None, 0.0, 0.02, 0.5, False))
        self.assertEqual(early.start_sample, 200)                      # never before the default
        self.assertEqual(ProcessingPlan.from_dict(plan.to_dict()), plan)
        self.assertEqual(plan.acquisition().start_sample, 203)


class PhaseSearchTests(unittest.TestCase):
    def test_global_search_then_fine_tune_on_resolved_lines(self):
        # edge_prior: the delay is known to +-0.3 ms from the raw FID (here: truth + 0.2 ms); the zero-order
        # phase is searched over its whole range. Without the edge the entropy optimum of this two-band spectrum
        # (lines at J and 2J) moves by about 4 ms (documented ambiguity), so the global mode is not tested here.
        from zulf_processing import phase_dataset
        for criterion in ("entropy", "lines"):
            for phase0, delay in [(1.1, 0.0008), (-2.5, -0.003)]:
                acq, fid = synthetic_fid(phase0, delay)
                f = np.arange(100.0, 300.0, acq.sampling_rate_hz / acq.points / 2)
                v = evaluate_spectrum(process_record(fid, acq), acq, f)
                r = phase_dataset(v, f, acq, edge_delay_s=delay + 0.0002, criterion=criterion,
                                  delay_span_s=0.0005)
                a = phase_correct(v, f, r.phase0_rad, r.delay_s, acq).real
                b = phase_correct(v, f, phase0, delay, acq).real        # truth (reference)
                corr = abs(np.dot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))
                self.assertGreater(corr, 0.97, msg=f"{criterion} {phase0} {delay}: {corr}")
                self.assertTrue(r.candidates and r.candidates[0]["cost"] == r.cost)

    def test_edge_fixed_mode_keeps_the_delay(self):
        from zulf_processing import phase_dataset
        acq, fid = synthetic_fid(1.1, 0.0008)
        f = np.arange(100.0, 300.0, acq.sampling_rate_hz / acq.points / 2)
        v = evaluate_spectrum(process_record(fid, acq), acq, f)
        r = phase_dataset(v, f, acq, edge_delay_s=0.0008, delay_mode="edge_fixed")
        self.assertEqual(r.delay_s, 0.0008)
        err = (r.phase0_rad - 1.1 + np.pi / 2) % np.pi - np.pi / 2       # sign is a convention
        self.assertLess(abs(err), 0.1)


class CalibratedPhaseTests(unittest.TestCase):
    def test_edge_plus_offset_and_calibrated_phase_recovers_truth(self):
        from zulf_processing import calibrated_phase
        acq, fid = synthetic_fid(1.1, 0.0008)
        r = calibrated_phase(0.00085, {"phase0_deg": np.degrees(1.1), "delay_offset_s": -5e-5})
        self.assertAlmostEqual(r.delay_s, 0.0008, places=12)            # edge + offset (arithmetic reference)
        self.assertAlmostEqual(r.phase0_rad, 1.1, places=12)
        f = np.arange(100.0, 300.0, acq.sampling_rate_hz / acq.points / 2)
        v = evaluate_spectrum(process_record(fid, acq), acq, f)
        a = phase_correct(v, f, r.phase0_rad, r.delay_s, acq)
        # truth phase: the corrected spectrum is the absorption spectrum (imaginary part small next to real)
        b = phase_correct(v, f, 1.1, 0.0008, acq)
        self.assertLess(np.linalg.norm(a - b), 1e-9 * np.linalg.norm(b))
        with self.assertRaises(ValueError):
            calibrated_phase(None, {"phase0_deg": 0.0})

    def test_diagnostics_are_fast(self):
        # regression: the multi-exponential fits of zulf_core.diagnose_fid took minutes per 64k-point FID
        import time
        from zulf_processing import diagnose_raw
        _, fid = synthetic_fid(1.1, 0.0008, points=65536, fs=4000.0)
        t = time.time()
        diagnose_raw(fid, 4000.0)
        self.assertLess(time.time() - t, 20.0)


class DatasetTests(unittest.TestCase):
    def test_process_dataset_records_every_choice(self):
        from zulf_processing import process_dataset
        acq, fid = synthetic_fid(1.1, 0.0008, points=8192, fs=1000.0)
        d = process_dataset(fid, 1000.0, "synthetic", defaults={"start_sample": 100, "stop_sample": 4100,
                                                                   "sg_window": 301, "apodization_rate_per_s": 0.0,
                                                                   "ranges": [[100.0, 300.0]],
                                                                   "instrument_lines_hz": []},
                            delay_mode="global", delay_span_s=0.005)
        rec = d.record()
        self.assertEqual(set(rec), {"sample_id", "diagnostics", "plan", "phase"})
        self.assertIn("start_sample", rec["plan"]["reasons"])
        self.assertTrue(d.phased.real_only)
        self.assertEqual(d.phased.phasing["delay_s"], d.phase.delay_s)
        with self.assertRaises(ValueError):                           # calibration needs a calibration
            process_dataset(fid, 1000.0, "synthetic", defaults={"ranges": [[100.0, 300.0]]},
                            phase_criterion="calibration")


if __name__ == "__main__":
    unittest.main()
