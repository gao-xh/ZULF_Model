"""Configurable trainer with device policy, AMP, schedules, checkpoints, curriculum and timing.

The trainer is model-agnostic: it needs a `CandidateModel`, a loader factory
that returns an iterable of collated batches for a curriculum stage, and an
optional validation dataset. Logs are JSON lines; checkpoints are
device-agnostic and resumable.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence

import numpy as np
import torch

from ..device import DevicePolicy, autocast, grad_scaler, seed_everything, select_device
from ..models import save_model
from ..models.base import CandidateModel
from zulf_core.timing import Timer
from .data import to_device
from .metrics import CoverageMetrics


@dataclass(frozen=True)
class CurriculumStage:
    until_step: int
    overrides: Dict[str, dict] = field(default_factory=dict)  # e.g. {"perturbation": {...}, "mixture": {...}}
    name: str = ""


@dataclass(frozen=True)
class TrainConfig:
    steps: int = 10000
    batch_size: int = 32
    learning_rate: float = 3e-4
    weight_decay: float = 0.01
    warmup_steps: int = 500
    min_lr_fraction: float = 0.05
    grad_clip: float = 1.0
    log_every: int = 50
    eval_every: int = 1000
    checkpoint_every: int = 1000
    eval_propose_k: int = 10
    eval_propose_limit: int = 64
    coverage_tolerance_hz: float = 1.0
    early_stopping_patience: int = 0          # evaluations without improvement; 0 disables
    seed: int = 0
    device: str = "auto"
    amp: Optional[str] = "auto"
    num_workers: int = 0
    renders_per_system: int = 4
    output_dir: str = "runs/train"
    curriculum: tuple = ()

    @classmethod
    def from_dict(cls, data: dict) -> "TrainConfig":
        data = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}
        if "curriculum" in data:
            data["curriculum"] = tuple(CurriculumStage(**s) for s in data["curriculum"])
        return cls(**data)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["curriculum"] = [asdict(s) for s in self.curriculum]
        return d


LoaderFactory = Callable[[CurriculumStage, int], Iterable[dict]]


class Trainer:
    def __init__(self, model: CandidateModel, loader_factory: LoaderFactory, config: TrainConfig,
                 val_batches: Optional[Sequence[dict]] = None, timer: Optional[Timer] = None,
                 extra_metadata: Optional[dict] = None):
        self.config = config
        self.policy: DevicePolicy = select_device(config.device, config.amp)
        seed_everything(config.seed)
        self.model = model.to(self.policy.device)
        self.loader_factory = loader_factory
        self.val_batches = list(val_batches or [])
        self.timer = timer or Timer()
        self.extra = extra_metadata or {}
        self.optimizer = torch.optim.AdamW(self.model.parameters(), lr=config.learning_rate,
                                           weight_decay=config.weight_decay)
        self.scheduler = torch.optim.lr_scheduler.LambdaLR(self.optimizer, self._lr_factor)
        self.scaler = grad_scaler(self.policy)
        self.step = 0
        self.best_metric = math.inf
        self.stale_evaluations = 0
        self.output = Path(config.output_dir)
        self.output.mkdir(parents=True, exist_ok=True)
        self.log_path = self.output / "log.jsonl"

    # -- schedule ------------------------------------------------------------------
    def _lr_factor(self, step: int) -> float:
        c = self.config
        if step < c.warmup_steps:
            return (step + 1) / max(1, c.warmup_steps)
        progress = (step - c.warmup_steps) / max(1, c.steps - c.warmup_steps)
        return c.min_lr_fraction + (1 - c.min_lr_fraction) * 0.5 * (1 + math.cos(math.pi * min(1.0, progress)))

    def stage_for(self, step: int) -> CurriculumStage:
        for stage in self.config.curriculum:
            if step < stage.until_step:
                return stage
        return CurriculumStage(until_step=self.config.steps, name="default")

    # -- logging and checkpoints -------------------------------------------------
    def log(self, record: dict) -> None:
        record = dict(record, step=self.step, time=time.time())
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")

    def save_checkpoint(self, name: str = "last.pt") -> Path:
        path = self.output / name
        state = {"model": {k: v.detach().cpu() for k, v in self.model.state_dict().items()},
                 "optimizer": self.optimizer.state_dict(), "scheduler": self.scheduler.state_dict(),
                 "step": self.step, "best_metric": self.best_metric, "train_config": self.config.to_dict(),
                 "spec": self.model.spec.to_dict(), "model_config": self.model.model_config.to_dict(),
                 "torch_rng": torch.get_rng_state(), "numpy_rng": np.random.get_state(), "extra": self.extra}
        torch.save(state, path)
        save_model(self.output / "model.pt", self.model, {"step": self.step, **self.extra})
        return path

    def resume(self, path) -> None:
        state = torch.load(path, map_location="cpu", weights_only=False)
        self.model.load_state_dict(state["model"])
        self.model.to(self.policy.device)
        self.optimizer.load_state_dict(state["optimizer"])
        self.scheduler.load_state_dict(state["scheduler"])
        self.step = state["step"]
        self.best_metric = state["best_metric"]
        torch.set_rng_state(state["torch_rng"])
        np.random.set_state(state["numpy_rng"])

    # -- evaluation --------------------------------------------------------------------
    @torch.no_grad()
    def evaluate(self) -> dict:
        if not self.val_batches:
            return {}
        self.model.eval()
        losses: Dict[str, List[float]] = {}
        metrics = CoverageMetrics(ks=sorted({1, 3, self.config.eval_propose_k}),
                                  tolerance_hz=self.config.coverage_tolerance_hz)
        proposed = 0
        with self.timer.section("eval"):
            for batch in self.val_batches:
                batch = to_device(batch, self.policy.device)
                with autocast(self.policy):
                    out = self.model(batch)
                    _, logs = self.model.loss(out, batch)
                for key, value in logs.items():
                    losses.setdefault(key, []).append(value)
                if proposed < self.config.eval_propose_limit:
                    take = min(len(batch["interpretations"]), self.config.eval_propose_limit - proposed)
                    candidates = self.model.propose(batch["features"][:take], batch["frequency_hz"],
                                                    k=self.config.eval_propose_k)
                    for truth, cands in zip(batch["interpretations"][:take], candidates):
                        metrics.update(truth, cands)
                    proposed += take
        self.model.train()
        summary = {f"val_{k}": float(np.mean(v)) for k, v in losses.items()}
        summary["coverage"] = metrics.summary()
        return summary

    # -- training loop -----------------------------------------------------------------
    def fit(self, max_steps: Optional[int] = None) -> dict:
        c = self.config
        target = min(c.steps, max_steps + self.step) if max_steps else c.steps
        self.model.train()
        stage = self.stage_for(self.step)
        loader = iter(self.loader_factory(stage, self.step))
        history = {"train": [], "eval": []}
        running: Dict[str, List[float]] = {}
        last_eval = {}
        data_start = time.perf_counter()
        while self.step < target:
            new_stage = self.stage_for(self.step)
            if new_stage != stage:
                stage = new_stage
                loader = iter(self.loader_factory(stage, self.step))
                self.log({"event": "curriculum_stage", "name": stage.name, "overrides": stage.overrides})
            batch = next(loader)
            self.timer.add("train.data_wait", time.perf_counter() - data_start)
            with self.timer.section("train.step"):
                batch = to_device(batch, self.policy.device)
                with autocast(self.policy):
                    out = self.model(batch)
                    loss, logs = self.model.loss(out, batch)
                self.optimizer.zero_grad(set_to_none=True)
                if self.scaler is not None:
                    self.scaler.scale(loss).backward()
                    self.scaler.unscale_(self.optimizer)
                else:
                    loss.backward()
                grad_norm = torch.nn.utils.clip_grad_norm_(self.model.parameters(), c.grad_clip)
                if self.scaler is not None:
                    self.scaler.step(self.optimizer)
                    self.scaler.update()
                else:
                    self.optimizer.step()
                self.scheduler.step()
            self.step += 1
            logs["grad_norm"] = float(grad_norm)
            logs["render_ms"] = batch.get("render_ms", 0.0)
            for key, value in logs.items():
                running.setdefault(key, []).append(value)
            if self.step % c.log_every == 0 or self.step == target:
                record = {k: float(np.mean(v)) for k, v in running.items()}
                record["lr"] = self.scheduler.get_last_lr()[0]
                record["stage"] = stage.name
                self.log({"event": "train", **record})
                history["train"].append(dict(record, step=self.step))
                running = {}
            if self.val_batches and (self.step % c.eval_every == 0 or self.step == target):
                last_eval = self.evaluate()
                self.log({"event": "eval", **last_eval})
                history["eval"].append(dict(last_eval, step=self.step))
                metric = last_eval.get("val_loss", math.inf)
                if metric < self.best_metric:
                    self.best_metric = metric
                    self.stale_evaluations = 0
                    save_model(self.output / "best_model.pt", self.model, {"step": self.step, **self.extra})
                else:
                    self.stale_evaluations += 1
                if c.early_stopping_patience and self.stale_evaluations >= c.early_stopping_patience:
                    self.log({"event": "early_stop"})
                    break
            if self.step % c.checkpoint_every == 0 or self.step == target:
                self.save_checkpoint()
            data_start = time.perf_counter()
        self.log({"event": "timing", "report": self.timer.report()})
        return {"steps": self.step, "history": history, "last_eval": last_eval, "timing": self.timer.report(),
                "device": self.policy.device, "amp": self.policy.amp_dtype}
