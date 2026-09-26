"""Training-data synthesis (perturbations, sample pipeline, features); re-exports zulf_core.render."""
from zulf_core.render.acquisition import Acquisition, evaluate_spectrum, process_record, spectrum_from_fid
from .features import spectrum_features, spectrum_scale
from zulf_core.render.grid import SpectrumGrid
from .perturb import (PerturbationConfig, RenderParams, render_observation, render_signal,
                      sample_render_params)
from .pipeline import ProcessingConfig, RenderedSample, SampleRenderer
from zulf_core.render.renderer import ContinuousRenderer, Renderer

__all__ = ["Acquisition", "process_record", "evaluate_spectrum", "spectrum_from_fid", "SpectrumGrid",
           "Renderer", "ContinuousRenderer", "PerturbationConfig", "RenderParams", "sample_render_params", "render_signal",
           "render_observation", "spectrum_features", "spectrum_scale", "ProcessingConfig",
           "RenderedSample", "SampleRenderer"]
