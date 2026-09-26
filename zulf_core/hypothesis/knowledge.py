"""Reference knowledge for hypotheses: measured couplings of known compounds.

A `KnowledgeBase` holds `ReferenceEntry` records: a compound, the fragment
template it was fitted with, named couplings (key names of that template,
e.g. 'J(Ca,Ha)'), conditions and a source. `nearest` ranks entries of the
same template by coupling distance. Every entry states its source; values
from literature must say so and are kept apart from measured ones
(`source_kind`), so blind analyses can exclude them.

The packaged file `data/reference_couplings.json` starts with the samples
confirmed in this project. Add entries with `add` and `save`, or load other
files with `KnowledgeBase.load(path)`; a database or literature source can
subclass `KnowledgeBase` and override `entries_for`.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

PACKAGED = Path(__file__).with_name("data") / "reference_couplings.json"


@dataclass
class ReferenceEntry:
    compound: str
    template: str                       # fragment template name (fragments.TEMPLATES)
    couplings: Dict[str, float]         # key name -> Hz
    source: str                         # where the numbers come from
    source_kind: str = "measured"       # measured | literature | computed
    conditions: str = ""                # solvent, pH, temperature, instrument
    uncertainty_hz: Dict[str, float] = field(default_factory=dict)
    notes: str = ""


class KnowledgeBase:
    def __init__(self, entries: Iterable[ReferenceEntry] = ()):
        self.entries: List[ReferenceEntry] = list(entries)

    @classmethod
    def load(cls, path=PACKAGED) -> "KnowledgeBase":
        data = json.loads(Path(path).read_text())
        return cls(ReferenceEntry(**e) for e in data["entries"])

    @classmethod
    def packaged(cls) -> "KnowledgeBase":
        return cls.load(PACKAGED)

    def save(self, path) -> None:
        Path(path).write_text(json.dumps({"entries": [asdict(e) for e in self.entries]}, indent=1) + "\n")

    def add(self, entry: ReferenceEntry) -> None:
        self.entries.append(entry)

    def entries_for(self, template: str, kinds: Optional[Sequence[str]] = None) -> List[ReferenceEntry]:
        return [e for e in self.entries if e.template == template and (kinds is None or e.source_kind in kinds)]

    def nearest(self, template: str, couplings: Mapping[str, float], keys: Optional[Sequence[str]] = None,
                kinds: Optional[Sequence[str]] = ("measured",), scale_hz: float = 1.0) -> List[Tuple[float, ReferenceEntry]]:
        """Entries of `template` ranked by RMS coupling difference over shared keys (divided by `scale_hz`)."""
        out = []
        for e in self.entries_for(template, kinds):
            shared = [k for k in (keys or couplings) if k in e.couplings and k in couplings]
            if not shared:
                continue
            d = math.sqrt(sum((couplings[k] - e.couplings[k]) ** 2 for k in shared) / len(shared)) / scale_hz
            out.append((d, e))
        return sorted(out, key=lambda t: t[0])
