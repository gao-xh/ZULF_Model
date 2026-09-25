"""Sharded sample storage (gzip JSON lines) with a manifest.

Spectra are rendered on the fly during training, so shards hold only ground
truth, provenance and (optionally) transition lists. HDF5 can be added behind
the same `write_shards` / `iter_samples` interface.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Optional

from .sampler import GENERATOR_VERSION, Sample


def write_shards(samples: Iterable[Sample], directory, shard_size: int = 1000,
                 manifest_extra: Optional[dict] = None, prefix: str = "shard") -> dict:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    shards, count, index, handle, path = [], 0, 0, None, None
    digest = None

    def close():
        nonlocal handle
        if handle is not None:
            handle.close()
            shards.append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            handle = None

    for sample in samples:
        if handle is None or count % shard_size == 0 and count:
            close()
            path = directory / f"{prefix}_{index:05d}.jsonl.gz"
            handle = gzip.open(path, "wt", encoding="utf-8")
            index += 1
        handle.write(json.dumps(sample.to_dict(), separators=(",", ":")) + "\n")
        count += 1
    close()
    manifest = {"generator_version": GENERATOR_VERSION, "count": count, "shards": shards,
                "created_utc": datetime.now(timezone.utc).isoformat()}
    manifest.update(manifest_extra or {})
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def verify_shards(directory) -> dict:
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    for shard in manifest["shards"]:
        actual = hashlib.sha256((directory / shard["file"]).read_bytes()).hexdigest()
        if actual != shard["sha256"]:
            raise ValueError(f"Shard {shard['file']} changed after writing.")
    return manifest


def iter_samples(directory, verify: bool = True) -> Iterator[Sample]:
    directory = Path(directory)
    manifest = verify_shards(directory) if verify else json.loads((directory / "manifest.json").read_text())
    for shard in manifest["shards"]:
        with gzip.open(directory / shard["file"], "rt", encoding="utf-8") as handle:
            for line in handle:
                yield Sample.from_dict(json.loads(line))
