"""End-to-end smoke test of the whole pipeline through the agent tool registry.

generate -> train (tiny) -> simulate an observation -> diagnose -> propose ->
refine. Runs on CPU in a few minutes. Proves that every stage connects; the
tiny model is not expected to propose correct candidates.

Usage: python scripts/smoke_pipeline.py [WORKSPACE]
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="zulf_smoke_"))
    os.environ["ZULF_MODEL_WORKSPACE"] = str(root)
    from zulf_model.agent import registry
    from zulf_model.generator import iter_samples
    from zulf_core.physics import compute_transitions
    from zulf_model.render import Acquisition, Renderer

    reg = registry()
    repo = Path(__file__).resolve().parents[1]
    report = {"workspace": str(root)}
    report["describe"] = reg.call("describe_project", {})["tools"]
    gen = reg.call("generate_samples", {"count": 6, "seed": 3, "split": None})
    report["generate"] = gen
    config = json.loads((repo / "configs" / "run_smoke.json").read_text())
    config["train"]["output_dir"] = str(root / "smoke_train")
    trained = reg.call("train_model", {"run_config": config})
    report["train"] = {k: trained[k] for k in ("steps", "device", "model_path")}

    sample = next(iter_samples(gen["directory"]))
    acq = Acquisition(1000.0, 2048)
    renderer = Renderer(acq)
    fid = sum(renderer.synthesize(compute_transitions(c.system), 2.0, gain=c.contribution)
              for c in sample.interpretation.components)
    fid = fid + np.random.default_rng(0).normal(0, 1e-4 * np.abs(fid).max(), acq.points)
    np.save(root / "observed.npy", fid)
    source = {"fid_npy": str(root / "observed.npy"), "sampling_rate_hz": 1000.0}
    report["diagnose"] = {k: v for k, v in reg.call("diagnose_fid", {"source": source}).items()
                          if k in ("plateau_end_s", "ringing_end_s", "first_point_anomaly")}
    proposed = reg.call("propose_candidates", {"model_path": trained["model_path"], "source": source, "k": 2})
    candidates = proposed["candidates"] + [sample.interpretation.to_dict()]
    refined = reg.call("refine_candidates", {"candidates": candidates, "source": source, "ranges": [[1.0, 450.0]],
                                             "settings": {"starts": 1, "max_seconds": 60,
                                                          "continuation_rates_per_s": [3.0, 0.0]}})
    report["refine"] = [{"candidate_index": r["candidate_index"], "relative_residual": r["relative_residual"],
                         "flags": r["flags"]} for r in refined["results"]]
    print(json.dumps(report, indent=2, default=float))


if __name__ == "__main__":
    main()
