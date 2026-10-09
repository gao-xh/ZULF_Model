"""Torch wrapper of PrerenderedDataset, in its own module so that DataLoader workers started with "spawn" (the
default on macOS) can pickle it; prerender.py stays importable without torch."""
from __future__ import annotations

from torch.utils.data import IterableDataset


class TorchPrerenderedDataset(IterableDataset):
    """A PrerenderedDataset as a torch IterableDataset (each worker reads its own shards)."""

    def __init__(self, inner):
        self.inner = inner

    def __iter__(self):
        return iter(self.inner)
