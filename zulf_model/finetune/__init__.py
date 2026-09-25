"""Closed-loop improvement: failure classification, focused resampling, fine-tuning rounds."""
from .failures import FailureRecord, classify_failure
from .focused import FocusedSource
from .loop import ActiveLearningLoop, LoopConfig

__all__ = ["FailureRecord", "classify_failure", "FocusedSource", "ActiveLearningLoop", "LoopConfig"]
