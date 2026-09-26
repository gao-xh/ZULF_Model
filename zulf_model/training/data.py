"""Torch datasets and collation.

Spectra are rendered on the fly from ground-truth interpretations, either drawn
live from a sampler or read from shards. Every worker derives its RNG from
(seed, epoch, worker id), so runs are reproducible and workers never share
streams. Datasets are picklable for macOS spawn-based DataLoader workers.
"""
from __future__ import annotations

import math
from typing import Dict, Iterator, List, Optional, Sequence

import numpy as np
import torch
from torch.utils.data import Dataset, IterableDataset, get_worker_info

from ..codec import PAD, InterpretationCodec
from ..generator.sampler import MixtureSampler, Sample
from ..render.pipeline import SampleRenderer
from zulf_core.spinsystem import Interpretation


def make_item(rendered, codec: InterpretationCodec, sample: Optional[Sample] = None) -> dict:
    """Targets use the rendered weights (`rendered.target`) so contribution targets match the spectrum."""
    target = rendered.target if getattr(rendered, "target", None) is not None else rendered.interpretation
    enc = codec.encode(target)
    item = {"features": rendered.features, "scale": rendered.scale, "tokens": enc.tokens,
            "j_offsets": enc.j_offsets, "interpretation": target,
            "family_id": sample.family_id if sample else "", "render_ms": rendered.timings_ms.get("total", 0.0)}
    item.update(codec.encode_set(target))
    return item


class Collator:
    """Turns items into a batch; keeps interpretations as a Python list."""

    def __init__(self, codec: InterpretationCodec, frequency_hz: np.ndarray):
        self.codec = codec
        self.frequency_hz = torch.as_tensor(np.asarray(frequency_hz), dtype=torch.float32)

    def __call__(self, items: Sequence[dict]) -> dict:
        batch = {"features": torch.as_tensor(np.stack([i["features"] for i in items])),
                 "frequency_hz": self.frequency_hz,
                 "scale": torch.as_tensor([i["scale"] for i in items], dtype=torch.float32)}
        for key in ("component_mask", "isotope", "spin_mask", "group_id", "couplings", "log10_contribution",
            "group_class", "group_mask", "group_couplings"):
            batch[key] = torch.as_tensor(np.stack([i[key] for i in items]))
        length = max(len(i["tokens"]) for i in items)
        pad = self.codec.vocab[PAD]
        tokens = np.full((len(items), length), pad, np.int64)
        offsets = np.zeros((len(items), length), np.float32)
        mask = np.zeros((len(items), length), np.float32)
        for row, item in enumerate(items):
            n = len(item["tokens"])
            tokens[row, :n] = item["tokens"]
            o = np.asarray(item["j_offsets"], float)
            mask[row, :n] = ~np.isnan(o)
            offsets[row, :n] = np.nan_to_num(o)
        batch.update(tokens=torch.as_tensor(tokens), j_offsets=torch.as_tensor(offsets),
                     j_offset_mask=torch.as_tensor(mask))
        batch["interpretations"] = [i["interpretation"] for i in items]
        batch["family_ids"] = [i["family_id"] for i in items]
        batch["render_ms"] = float(np.mean([i["render_ms"] for i in items]))
        return batch


def to_device(batch: dict, device) -> dict:
    return {k: (v.to(device, non_blocking=True) if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}


def _worker_rng(seed: int, epoch: int) -> np.random.Generator:
    info = get_worker_info()
    worker = info.id if info is not None else 0
    return np.random.default_rng([seed, epoch, worker])


class OnTheFlyDataset(IterableDataset):
    """Infinite stream: draw a system from the sampler, render it `renders_per_system` times."""

    def __init__(self, sampler: MixtureSampler, renderer: SampleRenderer, codec: InterpretationCodec,
                 split: Optional[str] = "train", seed: int = 0, renders_per_system: int = 4):
        self.sampler, self.renderer, self.codec = sampler, renderer, codec
        self.split, self.seed, self.renders = split, seed, renders_per_system
        self.epoch = 0

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[dict]:
        rng = _worker_rng(self.seed, self.epoch)
        while True:
            sample = self.sampler.draw(rng, self.split)
            for _ in range(self.renders):
                yield make_item(self.renderer.render(sample.interpretation, rng), self.codec, sample)


class ShardDataset(IterableDataset):
    """Iterates stored samples (list or shard directory) with fresh renders each epoch."""

    def __init__(self, samples: Sequence[Sample], renderer: SampleRenderer, codec: InterpretationCodec,
                 seed: int = 0, renders_per_system: int = 1, shuffle: bool = True, repeat: bool = True):
        self.samples = list(samples)
        self.renderer, self.codec = renderer, codec
        self.seed, self.renders, self.shuffle, self.repeat = seed, renders_per_system, shuffle, repeat
        self.epoch = 0

    def set_epoch(self, epoch: int) -> None:
        self.epoch = epoch

    def __iter__(self) -> Iterator[dict]:
        info = get_worker_info()
        workers, worker = (info.num_workers, info.id) if info is not None else (1, 0)
        epoch = self.epoch
        while True:
            rng = np.random.default_rng([self.seed, epoch, worker])
            order = np.arange(len(self.samples))
            if self.shuffle:
                np.random.default_rng([self.seed, epoch]).shuffle(order)
            for index in order[worker::workers]:
                sample = self.samples[int(index)]
                for _ in range(self.renders):
                    yield make_item(self.renderer.render(sample.interpretation, rng), self.codec, sample)
            if not self.repeat:
                return
            epoch += 1


class FixedDataset(Dataset):
    """Pre-rendered, deterministic dataset for validation and benchmarks."""

    def __init__(self, samples: Sequence[Sample], renderer: SampleRenderer, codec: InterpretationCodec,
                 seed: int = 12345):
        self.items = [make_item(renderer.render(s.interpretation, np.random.default_rng([seed, k])), codec, s)
                      for k, s in enumerate(samples)]

    def __len__(self) -> int:
        return len(self.items)

    def __getitem__(self, index: int) -> dict:
        return self.items[index]
