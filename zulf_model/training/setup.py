"""Assemble a training run from one JSON configuration.

Keys (each value is an inline dict or a path to a JSON file):
problem, generator, processing, perturbation, protocol, model, train, validation.
Curriculum stages may override `perturbation`, `processing` and `mixture`.

Optional `data: {"prerendered": DIR}` trains from pre-rendered shards
(`training.prerender`) instead of live rendering; a curriculum stage may name
its own shard directory with the override key `prerendered`. Relative paths
resolve against the run configuration.
"""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Iterable, List, Optional, Union

import numpy as np
from torch.utils.data import DataLoader

from ..codec import InterpretationCodec
from ..generator.sampler import GENERATOR_VERSION, MixtureConfig, build_default_sampler
from ..models import ModelConfig, build_model
from zulf_core.physics.protocol import Protocol
from ..render.perturb import PerturbationConfig
from ..render.pipeline import ProcessingConfig, SampleRenderer
from ..spec import ProblemSpec
from zulf_core.timing import Timer
from .data import Collator, FixedDataset, OnTheFlyDataset
from .trainer import CurriculumStage, TrainConfig, Trainer


def _load(value, base: Optional[Path] = None) -> dict:
    if value is None:
        return {}
    if isinstance(value, dict):
        return value
    path = Path(value)
    if base is not None and not path.is_absolute() and not path.exists():
        path = base / path
    return json.loads(path.read_text(encoding="utf-8"))


def _merge(obj, overrides: dict):
    if not overrides:
        return obj
    data = dataclasses.asdict(obj)
    data.update(overrides)
    return type(obj).from_dict(data)


class TrainingSetup:
    def __init__(self, config: Union[dict, str, Path]):
        base = None
        if not isinstance(config, dict):
            base = Path(config).resolve().parent
            config = json.loads(Path(config).read_text(encoding="utf-8"))
        self.raw = config
        self.base = base
        self.spec = ProblemSpec.from_dict(_load(config.get("problem"), base)) if config.get("problem") else ProblemSpec()
        self.generator_config = dict(_load(config.get("generator"), base))
        couplings = self.generator_config.get("couplings_path")
        if couplings and not Path(couplings).is_absolute() and base is not None:
            generator_ref = config.get("generator")
            gen_dir = (base / generator_ref).parent if isinstance(generator_ref, str) else base
            if (gen_dir / couplings).exists():
                self.generator_config["couplings_path"] = str(gen_dir / couplings)
        self.processing = ProcessingConfig.from_dict(_load(config.get("processing"), base))
        self.perturbation = PerturbationConfig.from_dict(_load(config.get("perturbation"), base))
        protocol = _load(config.get("protocol"), base)
        self.protocol = Protocol.from_dict(protocol) if protocol else Protocol()
        self.model_config = ModelConfig.from_dict(_load(config.get("model"), base))
        self.train_config = TrainConfig.from_dict(_load(config.get("train"), base))
        self.validation = dict(count=128, split="val", seed=99, batch_size=32)
        self.validation.update(_load(config.get("validation"), base))
        self.codec = InterpretationCodec(self.spec)
        self.timer = Timer()
        self.sampler = build_default_sampler(self.spec, self.generator_config)
        self.renderer = self.make_renderer({})

    def make_renderer(self, overrides: dict) -> SampleRenderer:
        return SampleRenderer(self.spec, _merge(self.processing, overrides.get("processing", {})),
                              _merge(self.perturbation, overrides.get("perturbation", {})), self.protocol,
                              timer=self.timer)

    def collator(self) -> Collator:
        return Collator(self.codec, self.renderer.grid.frequencies_hz)

    def prerendered_path(self, stage: CurriculumStage) -> Optional[Path]:
        value = stage.overrides.get("prerendered") or (self.raw.get("data") or {}).get("prerendered")
        if not value:
            return None
        path = Path(value)
        if not path.is_absolute() and self.base is not None and not path.exists():
            path = self.base / path
        return path

    def loader_factory(self, stage: CurriculumStage, step: int) -> Iterable[dict]:
        path = self.prerendered_path(stage)
        if path is not None:
            from .prerender import PrerenderedDataset, as_torch_dataset
            dataset = as_torch_dataset(PrerenderedDataset(path, self.spec.digest(), seed=self.train_config.seed + step))
            return DataLoader(dataset, batch_size=self.train_config.batch_size, collate_fn=self.collator(),
                              num_workers=self.train_config.num_workers,
                              persistent_workers=self.train_config.num_workers > 0)
        renderer = self.make_renderer(stage.overrides)
        sampler = build_default_sampler(self.spec, dict(self.generator_config,
                                                        mixture=dict(self.generator_config.get("mixture", {}),
                                                                     **stage.overrides.get("mixture", {}))))
        dataset = OnTheFlyDataset(sampler, renderer, self.codec, "train", self.train_config.seed + step,
                                  self.train_config.renders_per_system)
        return DataLoader(dataset, batch_size=self.train_config.batch_size, collate_fn=self.collator(),
                          num_workers=self.train_config.num_workers,
                          persistent_workers=self.train_config.num_workers > 0)

    def validation_batches(self) -> List[dict]:
        v = self.validation
        if not v.get("count"):
            return []
        samples = list(self.sampler.generate(int(v["count"]), seed=int(v["seed"]), split=v.get("split")))
        data = FixedDataset(samples, self.renderer, self.codec, seed=int(v["seed"]))
        collate = self.collator()
        size = int(v["batch_size"])
        return [collate([data[i] for i in range(start, min(start + size, len(data)))])
                for start in range(0, len(data), size)]

    def build_trainer(self, model=None) -> Trainer:
        model = model or build_model(self.spec, self.model_config)
        metadata = {"generator_version": GENERATOR_VERSION, "spec_digest": self.spec.digest(), "run_config": self.raw}
        return Trainer(model, self.loader_factory, self.train_config, self.validation_batches(), self.timer, metadata)
