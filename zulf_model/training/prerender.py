"""Pre-rendered training shards for accelerator-rate training.

Rendering costs tens of milliseconds per spectrum on one CPU core, which is
too slow to keep a GPU (CUDA or Apple MPS) busy. `prerender` renders a fixed
number of spectra with the same sampler, renderer and codec as live training
and writes them as shards; `PrerenderedDataset` streams them back with
shuffling. Live rendering remains the default and the reference.

Layout of a shard directory:

    manifest.json           spec, spec digest, run config, stage, split, seeds, shard list
    frequency_hz.npy        model grid
    shard_00000.npz         features (float16), scale, set targets, padded tokens and J offsets
    shard_00000.jsonl       one JSON line per item: interpretation (targets) and family id

Seeds are (seed, shard index), so shards can be produced in parallel,
resumed, or extended without overlap.
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence

import numpy as np

SET_KEYS = ("component_mask", "isotope", "spin_mask", "group_id", "couplings", "log10_contribution",
            "group_class", "group_mask", "group_couplings")
MANIFEST = "manifest.json"


@dataclass(frozen=True)
class PrerenderRequest:
    run_config: str               # path to a run configuration (relative files resolve against it)
    output: str
    count: int
    shard_size: int = 1024
    split: str = "train"
    seed: int = 0
    stage: str = ""               # curriculum stage whose overrides apply ("" = base configuration)
    workers: int = 1
    feature_dtype: str = "float16"


def _stage_overrides(setup, stage: str) -> dict:
    if not stage:
        return {}
    for s in setup.train_config.curriculum:
        if s.name == stage:
            return s.overrides
    raise ValueError(f"Unknown curriculum stage '{stage}'.")


def _render_shard(request: PrerenderRequest, index: int, count: int) -> dict:
    from ..generator.sampler import build_default_sampler
    from .data import make_item
    from .setup import TrainingSetup
    setup = TrainingSetup(request.run_config)
    overrides = _stage_overrides(setup, request.stage)
    renderer = setup.make_renderer(overrides)
    sampler = build_default_sampler(setup.spec, dict(setup.generator_config,
                                                     mixture=dict(setup.generator_config.get("mixture", {}),
                                                                  **overrides.get("mixture", {}))))
    rng = np.random.default_rng([request.seed, index])
    renders = max(1, setup.train_config.renders_per_system)
    items, lines = [], []
    start = time.perf_counter()
    while len(items) < count:
        sample = sampler.draw(rng, request.split)
        for _ in range(min(renders, count - len(items))):
            item = make_item(renderer.render(sample.interpretation, rng), setup.codec, sample)
            items.append(item)
            lines.append(json.dumps({"interpretation": item["interpretation"].to_dict(),
                                     "family_id": item["family_id"]}))
    length = max(len(i["tokens"]) for i in items)
    tokens = np.zeros((count, length), np.int32)
    offsets = np.full((count, length), np.nan, np.float32)
    lengths = np.zeros(count, np.int32)
    for row, item in enumerate(items):
        n = len(item["tokens"])
        tokens[row, :n] = item["tokens"]
        offsets[row, :n] = np.asarray(item["j_offsets"], np.float32)
        lengths[row] = n
    arrays = {"features": np.stack([i["features"] for i in items]).astype(request.feature_dtype),
              "scale": np.asarray([i["scale"] for i in items], np.float32),
              "tokens": tokens, "token_length": lengths, "j_offsets": offsets}
    for key in SET_KEYS:
        arrays[key] = np.stack([i[key] for i in items])
    out = Path(request.output)
    np.savez(out / f"shard_{index:05d}.npz", **arrays)
    (out / f"shard_{index:05d}.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"index": index, "count": count, "seconds": time.perf_counter() - start,
            "spec_digest": setup.spec.digest()}


def prerender(request: PrerenderRequest) -> dict:
    """Render `count` items into shards; existing shards with the same manifest are kept (resume)."""
    from ..generator.sampler import GENERATOR_VERSION
    from .setup import TrainingSetup
    setup = TrainingSetup(request.run_config)
    _stage_overrides(setup, request.stage)
    out = Path(request.output)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"spec": setup.spec.to_dict(), "spec_digest": setup.spec.digest(), "run_config": setup.raw,
                "run_config_path": str(Path(request.run_config).resolve()), "stage": request.stage,
                "split": request.split, "seed": request.seed, "shard_size": request.shard_size,
                "feature_dtype": request.feature_dtype, "channels": list(setup.spec.grid.channels),
                "generator_version": GENERATOR_VERSION, "shards": []}
    manifest_path = out / MANIFEST
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        same = all(previous.get(k) == manifest[k] for k in ("spec_digest", "stage", "split", "seed", "shard_size"))
        if not same:
            raise ValueError(f"{out} holds shards made with different settings; use a new directory.")
        manifest["shards"] = previous.get("shards", [])
    np.save(out / "frequency_hz.npy", setup.renderer.grid.frequencies_hz)
    done = {s["index"] for s in manifest["shards"] if (out / f"shard_{s['index']:05d}.npz").exists()}
    plan = []
    remaining, index = int(request.count), 0
    while remaining > 0:
        size = min(request.shard_size, remaining)
        if index not in done:
            plan.append((index, size))
        remaining -= size
        index += 1
    start = time.perf_counter()
    results = []
    if request.workers > 1 and len(plan) > 1:
        with ProcessPoolExecutor(max_workers=request.workers) as pool:
            futures = [pool.submit(_render_shard, request, i, n) for i, n in plan]
            results = [f.result() for f in futures]
    else:
        results = [_render_shard(request, i, n) for i, n in plan]
    manifest["shards"] = sorted([s for s in manifest["shards"] if s["index"] in done] +
                                [{"index": r["index"], "count": r["count"]} for r in results],
                                key=lambda s: s["index"])
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    elapsed = time.perf_counter() - start
    rendered = sum(r["count"] for r in results)
    return {"output": str(out), "rendered": rendered, "skipped_shards": len(done),
            "total": sum(s["count"] for s in manifest["shards"]), "seconds": elapsed,
            "items_per_second": rendered / elapsed if elapsed > 0 else None}


def read_manifest(path) -> dict:
    return json.loads((Path(path) / MANIFEST).read_text(encoding="utf-8"))


class PrerenderedDataset:
    """Iterable over pre-rendered shards; picklable for DataLoader workers.

    Each worker reads a disjoint subset of shards; shard order and item order
    are reshuffled every pass from (seed, pass). Interpretations are decoded
    only when `decode_interpretations` is set (validation and metrics).
    """

    def __init__(self, path, spec_digest: Optional[str] = None, seed: int = 0, shuffle: bool = True,
                 repeat: bool = True, decode_interpretations: bool = False):
        self.path = Path(path)
        self.manifest = read_manifest(self.path)
        if spec_digest is not None and self.manifest["spec_digest"] != spec_digest:
            raise ValueError(f"Shards in {self.path} were rendered for spec {self.manifest['spec_digest']}, "
                             f"not {spec_digest}.")
        self.shards = [s["index"] for s in self.manifest["shards"]]
        if not self.shards:
            raise ValueError(f"No shards in {self.path}.")
        self.seed, self.shuffle, self.repeat = seed, shuffle, repeat
        self.decode = decode_interpretations

    def __len__(self) -> int:
        return sum(s["count"] for s in self.manifest["shards"])

    def frequencies_hz(self) -> np.ndarray:
        return np.load(self.path / "frequency_hz.npy")

    def _load(self, index: int):
        with np.load(self.path / f"shard_{index:05d}.npz") as blob:
            arrays = {k: blob[k] for k in blob.files}
        meta = None
        if self.decode:
            meta = [json.loads(line) for line in
                    (self.path / f"shard_{index:05d}.jsonl").read_text(encoding="utf-8").splitlines() if line]
        return arrays, meta

    def _item(self, arrays: dict, meta, row: int) -> dict:
        from ..spinsystem import Interpretation
        n = int(arrays["token_length"][row])
        item = {"features": arrays["features"][row].astype(np.float32), "scale": float(arrays["scale"][row]),
                "tokens": arrays["tokens"][row, :n].astype(np.int64).tolist(),
                "j_offsets": arrays["j_offsets"][row, :n].astype(float).tolist(),
                "interpretation": Interpretation.from_dict(meta[row]["interpretation"]) if meta else None,
                "family_id": meta[row]["family_id"] if meta else "", "render_ms": 0.0}
        for key in SET_KEYS:
            item[key] = arrays[key][row]
        return item

    def __iter__(self) -> Iterator[dict]:
        try:
            from torch.utils.data import get_worker_info
            info = get_worker_info()
        except ImportError:
            info = None
        workers, worker = (info.num_workers, info.id) if info is not None else (1, 0)
        epoch = 0
        while True:
            order = np.array(self.shards)
            if self.shuffle:
                np.random.default_rng([self.seed, epoch]).shuffle(order)
            # Partition shards across workers; with fewer shards than workers, partition items instead
            # so that no worker is left without data.
            by_shard = len(order) >= workers
            counter = 0
            for index in (order[worker::workers] if by_shard else order):
                arrays, meta = self._load(int(index))
                rows = np.arange(len(arrays["scale"]))
                if self.shuffle:
                    np.random.default_rng([self.seed, epoch, int(index)]).shuffle(rows)
                for row in rows:
                    counter += 1
                    if by_shard or (counter - 1) % workers == worker:
                        yield self._item(arrays, meta, int(row))
            if not self.repeat:
                return
            epoch += 1


def as_torch_dataset(dataset: PrerenderedDataset):
    """Wrap as a torch IterableDataset (torch imported lazily)."""
    from torch.utils.data import IterableDataset

    class _Wrapped(IterableDataset):
        def __init__(self, inner):
            self.inner = inner

        def __iter__(self):
            return iter(self.inner)

    return _Wrapped(dataset)
