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

    def entries_for(self, template: str, kinds: Optional[Sequence[str]] = None,
                    exclude: Sequence[str] = ()) -> List[ReferenceEntry]:
        """`exclude`: text fragments (e.g. sample ids); entries whose source mentions one are skipped, so a
        blind analysis never matches against a reference fitted on the sample under test."""
        return [e for e in self.entries if e.template == template and (kinds is None or e.source_kind in kinds)
                and not any(x and x in e.source for x in exclude)]

    def nearest(self, template: str, couplings: Mapping[str, float], keys: Optional[Sequence[str]] = None,
                kinds: Optional[Sequence[str]] = ("measured",), scale_hz: float = 1.0,
                exclude: Sequence[str] = ()) -> List[Tuple[float, ReferenceEntry]]:
        """Entries of `template` ranked by RMS coupling difference over shared keys (divided by `scale_hz`)."""
        out = []
        for e in self.entries_for(template, kinds, exclude):
            shared = [k for k in (keys or couplings) if k in e.couplings and k in couplings]
            if not shared:
                continue
            d = math.sqrt(sum((couplings[k] - e.couplings[k]) ** 2 for k in shared) / len(shared)) / scale_hz
            out.append((d, e))
        return sorted(out, key=lambda t: t[0])


def match_template(fragment, template_fragment) -> Optional[Dict[str, str]]:
    """Injective map from template labels onto fragment labels (element, proton-group sizes and bonds kept).

    The template may be a sub-fragment (the fragment can carry extra sites or protons). Returns None when no
    map exists; the first map in label order otherwise.
    """
    from itertools import permutations
    if fragment is None:
        return None
    t_sites = [s.label for s in template_fragment.sites
               if any(p.site == s.label for p in template_fragment.protons)]
    f_sites = [s.label for s in fragment.sites]
    t_bonds = {frozenset(b) for b in template_fragment.bonds}
    f_bonds = {frozenset(b) for b in fragment.bonds}
    protons_on = lambda frag, site: sorted(p.size for p in frag.protons if p.site == site)
    for image in permutations(f_sites, len(t_sites)):
        m = dict(zip(t_sites, image))
        if any(template_fragment.site(t).element != fragment.site(f).element for t, f in m.items()):
            continue
        if any(protons_on(template_fragment, t) != protons_on(fragment, f) for t, f in m.items()):
            continue
        if any(frozenset(m[x] for x in b) not in f_bonds for b in t_bonds if all(x in m for x in b)):
            continue
        for t, f in list(m.items()):         # proton groups follow their sites (sizes already equal)
            tp = [p for p in template_fragment.protons if p.site == t]
            fp = [p for p in fragment.protons if p.site == f]
            for a, b in zip(sorted(tp, key=lambda p: p.size), sorted(fp, key=lambda p: p.size)):
                m[a.label] = b.label
        return m
    return None


def translate_couplings(named: Mapping[str, float], mapping: Mapping[str, str], template_fragment) -> Dict[str, float]:
    """Refined couplings (fragment key names) renamed to the template's key names where both labels map."""
    from .fragment import pair
    inverse = {v: k for k, v in mapping.items()}
    orbits = template_fragment.coupling_orbits()
    out = {}
    for key, value in named.items():
        if not key.startswith("J(") or ":" in key:
            continue
        a, b = key[2:-1].split(",")
        if a in inverse and b in inverse:
            out[template_fragment.key_name(orbits[pair(inverse[a], inverse[b])])] = value
    return out


def knowledge_matches(model, named: Mapping[str, float], kb: "KnowledgeBase", kinds=("measured",),
                      top: int = 3, exclude: Sequence[str] = ()) -> List[dict]:
    """Nearest reference compounds over every template that maps onto the model's fragment."""
    from .fragments import TEMPLATES
    out = []
    for name, factory in TEMPLATES.items():
        template = factory()
        mapping = match_template(getattr(model, "fragment", None), template)
        if mapping is None:
            continue
        translated = translate_couplings(named, mapping, template)
        for dist, entry in kb.nearest(name, translated, kinds=kinds, exclude=exclude)[:top]:
            out.append({"template": name, "compound": entry.compound, "rms_hz": round(dist, 3),
                        "keys": sorted(k for k in translated if k in entry.couplings), "source": entry.source_kind})
    return sorted(out, key=lambda d: d["rms_hz"])[:top]
