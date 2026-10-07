"""Extensible registry of nuclei.

gamma values are gamma / (2 pi) in Hz per microtesla. The spin quantum number
is stored as a fraction so that quadrupolar nuclei can be added later without
changing the physics layer, which handles any spin quantum number.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Iterable, Tuple


@dataclass(frozen=True)
class Nucleus:
    symbol: str
    spin: Fraction
    gamma_hz_per_ut: float
    natural_abundance: float
    element: str

    @property
    def multiplicity(self) -> int:
        """Number of Zeeman states, 2 S + 1."""
        return int(2 * self.spin + 1)


class NucleusRegistry:
    """Mapping from symbol to `Nucleus`. Symbols are case-sensitive, e.g. '13C'."""

    def __init__(self, nuclei: Iterable[Nucleus] = ()):
        self._nuclei: Dict[str, Nucleus] = {}
        for nucleus in nuclei:
            self.register(nucleus)

    def register(self, nucleus: Nucleus, replace: bool = False) -> None:
        if nucleus.symbol in self._nuclei and not replace:
            raise ValueError(f"Nucleus {nucleus.symbol} already registered.")
        if nucleus.spin <= 0 or (2 * nucleus.spin).denominator != 1:
            raise ValueError("Spin quantum number must be a positive multiple of 1/2.")
        if not (abs(nucleus.gamma_hz_per_ut) > 0):
            raise ValueError("gamma must be finite and nonzero.")
        if not 0 <= nucleus.natural_abundance <= 1:
            raise ValueError("Natural abundance must lie in [0, 1].")
        self._nuclei[nucleus.symbol] = nucleus

    def __getitem__(self, symbol: str) -> Nucleus:
        try:
            return self._nuclei[symbol]
        except KeyError:
            raise KeyError(f"Unknown nucleus '{symbol}'. Registered: {sorted(self._nuclei)}") from None

    def __contains__(self, symbol: str) -> bool:
        return symbol in self._nuclei

    def symbols(self) -> Tuple[str, ...]:
        return tuple(self._nuclei)

    def gamma(self, symbol: str) -> float:
        return self[symbol].gamma_hz_per_ut

    def spin(self, symbol: str) -> Fraction:
        return self[symbol].spin

    def to_dict(self) -> dict:
        return {s: {"spin": str(n.spin), "gamma_hz_per_ut": n.gamma_hz_per_ut,
                    "natural_abundance": n.natural_abundance, "element": n.element}
                for s, n in self._nuclei.items()}


# gamma/(2 pi) values in Hz/uT (CODATA-derived, rounded) and natural abundances.
# 15N has a negative gamma and spin +1/2. 14N and 2H are registered for
# bookkeeping only; v1 generation does not place them in spin systems.
DEFAULT_NUCLEI = (
    Nucleus("1H", Fraction(1, 2), 42.577478, 0.999885, "H"),
    Nucleus("13C", Fraction(1, 2), 10.708395, 0.0107, "C"),
    Nucleus("15N", Fraction(1, 2), -4.316377, 0.00364, "N"),
    Nucleus("19F", Fraction(1, 2), 40.078, 1.0, "F"),
    Nucleus("31P", Fraction(1, 2), 17.235, 1.0, "P"),
    Nucleus("2H", Fraction(1), 6.536, 0.000115, "H"),
    Nucleus("14N", Fraction(1), 3.077, 0.99636, "N"),
)

REGISTRY = NucleusRegistry(DEFAULT_NUCLEI)


def get_registry() -> NucleusRegistry:
    """Return the process-wide default registry."""
    return REGISTRY


# ---- reduced coupling constants ---------------------------------------------------------------
# K_AB = 4 pi^2 J_AB / (h gamma_A gamma_B) removes the gyromagnetic ratios from a scalar coupling, leaving the
# electronic part: couplings of different nuclei (C-H, N-H, C-N) can be compared, and a coupling measured with one
# isotope converts to another (15N -> 14N, 1H -> 2H). With gamma / (2 pi) = g in Hz/T this is K = J / (h g_A g_B).
# Converted values neglect the primary isotope effect on K (of order 1 %): use them as starts or priors, then fit.

PLANCK_J_S = 6.62607015e-34              # exact (SI 2019)
REDUCED_COUPLING_UNIT = 1.0e19           # literature unit of K: 1e19 N A^-2 m^-3 (= 1e19 T^2 J^-1)


def _g_hz_per_t(symbol: str, registry: "NucleusRegistry" = None) -> float:
    return (registry or REGISTRY).gamma(symbol) * 1.0e6


def reduced_coupling(j_hz: float, nucleus_a: str, nucleus_b: str, registry: "NucleusRegistry" = None) -> float:
    """Reduced coupling constant K (N A^-2 m^-3) of a coupling J (Hz) between nuclei A and B."""
    return float(j_hz) / (PLANCK_J_S * _g_hz_per_t(nucleus_a, registry) * _g_hz_per_t(nucleus_b, registry))


def coupling_from_reduced(k: float, nucleus_a: str, nucleus_b: str, registry: "NucleusRegistry" = None) -> float:
    """Coupling J (Hz) between nuclei A and B for a reduced coupling K (N A^-2 m^-3)."""
    return float(k) * PLANCK_J_S * _g_hz_per_t(nucleus_a, registry) * _g_hz_per_t(nucleus_b, registry)


def convert_coupling(j_hz: float, pair_from: Tuple[str, str], pair_to: Tuple[str, str],
                     registry: "NucleusRegistry" = None) -> float:
    """J of the same bond path for other isotopes at equal K, e.g. J(15N,1H) -> J(14N,1H) (factor -0.713) or
    J(13C,1H) -> J(13C,2H) (factor 0.1535)."""
    return coupling_from_reduced(reduced_coupling(j_hz, *pair_from, registry=registry), *pair_to, registry=registry)


def nearest_nuclei(gamma_hz_per_ut: float, spin=None, registry: "NucleusRegistry" = None):
    """Registered nuclei ordered by |gamma - fitted| / |fitted| (optionally only those of one spin), as
    [(symbol, gamma, relative difference)]: names the nucleus a fitted gyromagnetic ratio points to (D53)."""
    reg = registry or REGISTRY
    rows = []
    for s in reg.symbols():
        if spin is not None and reg.spin(s) != Fraction(spin):
            continue
        g = reg.gamma(s)
        rows.append((s, g, abs(g - gamma_hz_per_ut) / max(abs(gamma_hz_per_ut), 1e-12)))
    return sorted(rows, key=lambda r: r[2])
