"""Device selection and precision policy for CUDA, Apple MPS and CPU.

Simulation and rendering always run in NumPy float64 on the CPU (MPS has no
float64). Networks train in float32; automatic mixed precision is optional and
enabled by default only on CUDA. Checkpoints are saved device-agnostic.
"""
from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class DevicePolicy:
    device: str
    amp_dtype: Optional[str]
    supports_float64: bool
    supports_compile: bool
    pin_memory: bool

    def torch_device(self):
        import torch
        return torch.device(self.device)


def available_devices() -> dict:
    import torch
    mps = getattr(torch.backends, "mps", None)
    return {"cuda": torch.cuda.is_available(),
            "cuda_count": torch.cuda.device_count() if torch.cuda.is_available() else 0,
            "mps": bool(mps and mps.is_available()),
            "torch": torch.__version__}


def select_device(preference: str = "auto", amp: Optional[str] = "auto") -> DevicePolicy:
    """Choose a device: 'auto' tries cuda, then mps, then cpu.

    `amp` is 'auto', None/'off', 'bf16' or 'fp16'. 'auto' gives bf16 on CUDA
    devices that support it (fp16 otherwise) and off on MPS and CPU.
    """
    import torch
    info = available_devices()
    if preference == "auto":
        device = "cuda" if info["cuda"] else "mps" if info["mps"] else "cpu"
    else:
        device = preference
        base = device.split(":")[0]
        if base == "cuda" and not info["cuda"]:
            raise RuntimeError("CUDA requested but not available.")
        if base == "mps" and not info["mps"]:
            raise RuntimeError("MPS requested but not available.")
        if base not in ("cuda", "mps", "cpu"):
            raise ValueError("Device must be auto, cpu, mps or cuda[:index].")
    base = device.split(":")[0]
    if base == "mps":
        os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    if amp == "auto":
        if base == "cuda":
            amp_dtype = "bf16" if torch.cuda.is_bf16_supported() else "fp16"
        else:
            amp_dtype = None
    elif amp in (None, "off", "none"):
        amp_dtype = None
    elif amp in ("bf16", "fp16"):
        amp_dtype = amp
    else:
        raise ValueError("amp must be auto, off, bf16 or fp16.")
    return DevicePolicy(device, amp_dtype, base != "mps", base == "cuda", base == "cuda")


def autocast(policy: DevicePolicy):
    """Autocast context for the forward pass; a no-op when AMP is off."""
    if policy.amp_dtype is None:
        return contextlib.nullcontext()
    import torch
    dtype = torch.bfloat16 if policy.amp_dtype == "bf16" else torch.float16
    return torch.autocast(device_type=policy.device.split(":")[0], dtype=dtype)


def grad_scaler(policy: DevicePolicy):
    """GradScaler for fp16 on CUDA; None otherwise."""
    import torch
    if policy.amp_dtype == "fp16" and policy.device.startswith("cuda"):
        return torch.amp.GradScaler("cuda")
    return None


def seed_everything(seed: int) -> None:
    import random
    import numpy as np
    import torch
    random.seed(seed)
    np.random.seed(seed % 2**32)
    torch.manual_seed(seed)
