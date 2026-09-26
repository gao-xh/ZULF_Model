"""Core of the ZULF toolchain: spin systems, zero-field physics, rendering through the
experimental processing operator, refinement against spectra, and experimental data I/O.

Depends on NumPy and SciPy only. The learning side (`zulf_model`: generator, codec,
networks, training) depends on this package, never the other way round, so training
data and experimental fits share one forward model.
"""
__version__ = "0.1.0"
