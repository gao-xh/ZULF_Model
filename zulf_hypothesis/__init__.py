"""zulf_hypothesis: programmatic chemical hypotheses for ZULF spectra.

A separate package: it uses zulf_core (physics, the public solver API
`refine` / `RefineSettings`, spectra) and never changes it; neural models
(zulf_model) enter only through `adapters.model_hints` as hint providers.

Layers (each a registry or interface meant for extension):

* `fragment`  - label-based fragment descriptions with symmetry;
* `fragments` - templates used so far (`register_template` adds more);
* `labeling`  - natural abundance (default) or enriched labelling schemes;
* `builder`   - fragment -> isotopologue set with automatic ties, abundance
  ratios, minor isotopologues tied to parents, omission reasons
  (`HypothesisModel`);
* `checks`    - physical checks on refined results (`register_check`);
* `moves`     - one-step model extensions (`register_move`);
* `knowledge` - measured / literature reference couplings (`KnowledgeBase`).

* `inventory` - observed lines grouped into bands, with a support function;
* `groups`    - X-Hn group candidates per band (`register_pattern`);
* `hints`     - outside suggestions, e.g. neural proposers, as search hints
  (`HintProvider`, `ProposerHints`, `register_hint_provider`);
* `enumerate` - band covers, topologies and equivalence -> fragments, with
  coupling priors (`register_prior`);
* `pipeline`  - `propose_hypotheses`: spectrum -> built, ranked hypotheses;
* `scoring`   - one yardstick for all hypotheses (data cores, noise, BIC);
* `search`    - `search_hypotheses`: refine (free / fixed variants, parallel),
  check, extend by triggered moves, rank, match known compounds.
"""
from .builder import HypothesisModel, build_model, combine_models, min_ratio_for_snr, model_from_interpretation
from .checks import CHECKS, CheckContext, Finding, register_check, run_checks
from .fragment import DEFAULT_LABEL_ISOTOPES, Fragment, ProtonGroup, Site, pair
from .fragments import TEMPLATES, register_template, template
from .enumerate import PRIORS, CouplingPrior, FragmentProposal, GenericSp3Prior, enumerate_fragments, register_prior
from .groups import PATTERNS, GroupCandidate, GroupPattern, group_candidates, register_pattern
from .hints import (HINT_PROVIDERS, Hint, HintProvider, ProposerHints, apply_hints, group_hints_from_interpretation,
                    insight_report, register_hint_provider)
from .inventory import Band, Inventory, Line, band_inventory
from .knowledge import KnowledgeBase, ReferenceEntry
from .labeling import Labeling
from .pipeline import ProposalSet, propose_hypotheses
from .scoring import Yardstick, yardstick
from .search import Evaluated, SearchResult, SearchSettings, search_hypotheses
from .moves import MOVES, AddCoupledProton, ExtensionMove, propose_all, register_move

__all__ = ["HypothesisModel", "build_model", "combine_models", "min_ratio_for_snr", "model_from_interpretation",
           "Yardstick", "yardstick", "Evaluated", "SearchResult", "SearchSettings", "search_hypotheses", "PRIORS", "CouplingPrior", "FragmentProposal",
           "GenericSp3Prior", "enumerate_fragments", "register_prior", "PATTERNS", "GroupCandidate", "GroupPattern",
           "group_candidates", "register_pattern", "HINT_PROVIDERS", "Hint", "HintProvider", "ProposerHints",
           "apply_hints", "group_hints_from_interpretation", "insight_report", "register_hint_provider", "Band",
           "Inventory", "Line", "band_inventory", "ProposalSet", "propose_hypotheses", "Labeling", "CHECKS", "CheckContext", "Finding", "register_check", "run_checks",
           "DEFAULT_LABEL_ISOTOPES", "Fragment", "ProtonGroup", "Site", "pair", "TEMPLATES", "register_template",
           "template", "KnowledgeBase", "ReferenceEntry", "MOVES", "AddCoupledProton", "ExtensionMove",
           "propose_all", "register_move"]
