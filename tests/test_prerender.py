import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_training import tiny_run_config  # noqa: E402

from zulf_model.training import TrainingSetup  # noqa: E402
from zulf_model.training.data import Collator  # noqa: E402
from zulf_model.training.prerender import (PrerenderRequest, PrerenderedDataset, as_torch_dataset,  # noqa: E402
                                           prerender, read_manifest)


class PrerenderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        cls.config_path = root / "run.json"
        cls.config_path.write_text(json.dumps(tiny_run_config(str(root / "out"), steps=4)), encoding="utf-8")
        cls.shards = root / "shards"
        cls.result = prerender(PrerenderRequest(str(cls.config_path), str(cls.shards), count=10, shard_size=4, seed=3))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_layout_resume_and_manifest(self):
        self.assertEqual(self.result["rendered"], 10)
        manifest = read_manifest(self.shards)
        self.assertEqual([s["count"] for s in manifest["shards"]], [4, 4, 2])
        again = prerender(PrerenderRequest(str(self.config_path), str(self.shards), count=10, shard_size=4, seed=3))
        self.assertEqual(again["rendered"], 0)
        self.assertEqual(again["skipped_shards"], 3)
        with self.assertRaises(ValueError):
            prerender(PrerenderRequest(str(self.config_path), str(self.shards), count=10, shard_size=5, seed=3))

    def test_items_match_live_format_and_collate(self):
        setup = TrainingSetup(str(self.config_path))
        data = PrerenderedDataset(self.shards, setup.spec.digest(), repeat=False, decode_interpretations=True)
        items = list(data)
        self.assertEqual(len(items), len(data))
        self.assertEqual(items[0]["features"].dtype, np.float32)
        self.assertEqual(items[0]["features"].shape[1], len(data.frequencies_hz()))
        np.testing.assert_allclose(data.frequencies_hz(), setup.renderer.grid.frequencies_hz)
        batch = Collator(setup.codec, data.frequencies_hz())(items[:3])
        self.assertEqual(tuple(batch["features"].shape[:2]), (3, len(setup.spec.grid.channels)))
        self.assertIsNotNone(batch["interpretations"][0])
        # Token sequences decode to the stored target interpretation.
        decoded = setup.codec.decode(items[0]["tokens"], items[0]["j_offsets"])
        self.assertEqual(len(decoded.components), len(items[0]["interpretation"].components))

    def test_workers_see_disjoint_complete_data(self):
        from torch.utils.data import DataLoader
        setup = TrainingSetup(str(self.config_path))
        for workers in (2, 5):
            data = as_torch_dataset(PrerenderedDataset(self.shards, repeat=False, shuffle=False))
            loader = DataLoader(data, batch_size=None, num_workers=workers, collate_fn=lambda item: item)
            scales = sorted(float(item["scale"]) for item in loader)
            reference = sorted(float(item["scale"]) for item in PrerenderedDataset(self.shards, repeat=False))
            self.assertEqual(scales, reference)

    def test_training_from_shards_and_digest_check(self):
        config = json.loads(self.config_path.read_text(encoding="utf-8"))
        config["data"] = {"prerendered": str(self.shards)}
        config["train"]["output_dir"] = str(Path(self.tmp.name) / "from_shards")
        setup = TrainingSetup(config)
        result = setup.build_trainer().fit()
        self.assertEqual(result["steps"], 4)
        with self.assertRaises(ValueError):
            PrerenderedDataset(self.shards, spec_digest="0" * 16)


if __name__ == "__main__":
    unittest.main()
