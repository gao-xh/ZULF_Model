import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import figure_tools  # noqa: E402


class FigureToolsTests(unittest.TestCase):
    def test_palette_is_distinct_hex(self):
        self.assertEqual(len(set(figure_tools.PALETTE)), len(figure_tools.PALETTE))
        for c in figure_tools.PALETTE:
            self.assertRegex(c, r"^#[0-9a-f]{6}$")

    def test_isotopologue_png(self):
        try:
            import rdkit  # noqa: F401
        except ImportError:
            self.skipTest("RDKit not installed")
        png = figure_tools.isotopologue_png("CC(N)C", 1, figure_tools.PALETTE[1], size=(200, 160))
        self.assertTrue(png.startswith(b"\x89PNG"))
        import io
        import matplotlib.image as mpimg
        img = mpimg.imread(io.BytesIO(png), format="png")
        self.assertEqual(img.shape[:2], (160, 200))
        self.assertGreater(float(img[..., 3].max()), 0.5)        # something drawn on a transparent background


if __name__ == "__main__":
    unittest.main()
