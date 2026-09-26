"""Adapters from other packages into zulf_hypothesis (imported lazily; optional dependencies)."""


def model_hints(*args, **kwargs):
    """See `model_hints.model_hints` (needs zulf_model and torch)."""
    from .model_hints import model_hints as _model_hints
    return _model_hints(*args, **kwargs)


__all__ = ["model_hints"]
