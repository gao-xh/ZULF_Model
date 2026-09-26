"""Programmatic chemical hypotheses for ZULF spectra.

Layers (each a registry or interface meant for extension):

* `fragment`  - label-based fragment descriptions with symmetry;
* `fragments` - templates used so far (`register_template` adds more);
* `builder`   - fragment -> natural-abundance isotopologue set with automatic
  ties, abundance ratios and omission reasons (`HypothesisModel`);
* `checks`    - physical checks on refined results (`register_check`);
* `moves`     - one-step model extensions (`register_move`);
* `knowledge` - measured / literature reference couplings (`KnowledgeBase`).

Planned (see docs/PLAN.md): band inventory and group candidates from a
spectrum, fragment enumeration, and a budgeted search loop over moves.
"""
from .builder import HypothesisModel, build_model
from .checks import CHECKS, CheckContext, Finding, register_check, run_checks
from .fragment import DEFAULT_LABEL_ISOTOPES, Fragment, ProtonGroup, Site, pair
from .fragments import TEMPLATES, register_template, template
from .knowledge import KnowledgeBase, ReferenceEntry
from .moves import MOVES, AddCoupledProton, ExtensionMove, propose_all, register_move

__all__ = ["HypothesisModel", "build_model", "CHECKS", "CheckContext", "Finding", "register_check", "run_checks",
           "DEFAULT_LABEL_ISOTOPES", "Fragment", "ProtonGroup", "Site", "pair", "TEMPLATES", "register_template",
           "template", "KnowledgeBase", "ReferenceEntry", "MOVES", "AddCoupledProton", "ExtensionMove",
           "propose_all", "register_move"]
