"""Optional cross-check against ZULF_Analysis_Tools (skipped when that repository is absent).

Set ZULF_ANALYSIS_TOOLS to its path, or place it next to this repository.
"""
import importlib.util
import os
import sys
import unittest
from pathlib import Path

import numpy as np

from zulf_model.physics import compute_transitions

REPO = Path(__file__).resolve().parents[1]
CANDIDATES = [os.environ.get("ZULF_ANALYSIS_TOOLS", ""), str(REPO.parent / "ZULF_Analysis_Tools")]
TOOLS = next((p for p in CANDIDATES if p and (Path(p) / "zulf_tools" / "jfit.py").exists()), None)


@unittest.skipIf(TOOLS is None, "ZULF_Analysis_Tools not available")
class AnalysisToolsCrossCheck(unittest.TestCase):
    def test_isopropylamine_transitions_match(self):
        sys.path.insert(0, TOOLS)
        from zulf_tools import jfit
        spec = importlib.util.spec_from_file_location("iso", REPO / "scripts" / "refine_isopropylamine.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        p = dict(jfit.DEFAULT)
        systems = {"methine": module.methine(p["J_CH_methine"], p["J_HH_vicinal"], p["J_Cmethine_Hmethyl"]),
                   "methyl": module.methyl(p["J_CH_methyl"], p["J_HH_vicinal"], p["J_Cmethyl_Hmethine"],
                                           p["J_Cmethyl_Hother_methyl"])}
        for name, system in systems.items():
            f_old, w_old = jfit.transitions(p, name)
            tl = compute_transitions(system)
            keep = np.abs(tl.amplitudes) > 1e-9 * np.abs(tl.amplitudes).max()
            f_new, w_new = tl.frequencies_hz[keep], tl.amplitudes[keep].real
            order_old, order_new = np.argsort(f_old), np.argsort(f_new)
            np.testing.assert_allclose(f_new[order_new], f_old[order_old], atol=1e-6)
            np.testing.assert_allclose(w_new[order_new] / np.abs(w_new).sum(), w_old[order_old] / np.abs(w_old).sum(),
                                       atol=1e-12)


if __name__ == "__main__":
    unittest.main()
