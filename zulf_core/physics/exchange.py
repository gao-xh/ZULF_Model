"""Zero-field spectra with chemical exchange of single spins with a solvent pool (Liouville space).

Model (high-temperature, linear in polarization; D45). The molecule's deviation density operator delta
(rho = (1 + delta) / N) evolves as

    d delta / dt = -2 pi i [H, delta] + sum_x k_x ( Tr_x(delta) (x) 1_x / d_x + eps_x I_z,x - delta )

for every exchanging spin x (rate k_x in 1/s): the spin is replaced by a solvent spin that carries no
correlation with the molecule and the single-spin deviation eps_x I_z (eps_x: its preparation weight,
by default the same as the molecule's spin, i.e. the solvent was prepolarized alike; its polarization is
static at zero field). The pool is large and not observed through the molecule. k = 0 is the static spin
system (physics.transitions); k -> infinity decouples x.

Only the part of operator space the signal can reach is used: zero-field H and the exchange map commute
with total rotations and conserve total M, and the preparation, the source term and the detection are
rank-1, q = 0 operators, so the dynamics are restricted to the rank-1 subspace of the q = 0 operators
(product basis |a><b| with M_a = M_b), found once per spin list as the eigenvalue-2 eigenspace of the
total-spin Casimir superoperator. With L the Liouvillian on that subspace (eigenvalues lambda_j =
-R_j + 2 pi i f_j) the signal is

    s(t) = norm Tr[D delta(t)] = sum_j w_j g_j exp(lambda_j t) + constant,
    g = V^-1 delta0 + V^-1 source / lambda,

and every mode with f_j > 0 becomes a line with amplitude 2 w_j g_j (its conjugate partner supplies the
mirror term) and intrinsic decay R_j (TransitionList.line_rates). Non-oscillating modes are not rendered
(their weight is reported in the metadata); only the static part enters `dc`.

Cost: the rank-1 q = 0 dimension is 1001 for 7 spin-1/2 and 3432 for 8; practical up to about 8 spins.

Linear-algebra backend (optional GPU): the dense projection of the Liouvillian, its eigen-decomposition and
the Casimir basis run in NumPy/SciPy by default, or in PyTorch (`set_backend("torch", "cuda")`, or the
environment variable ZULF_LINALG_DEVICE=cuda / cpu / numpy read at import). PyTorch is optional and only
imported when selected; results agree with the NumPy route to rounding (tests/test_exchange.py).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

from ..nuclei import get_registry
from ..spinsystem import SpinSystem
from .operators import angular_momentum
from .protocol import SUDDEN_DROP, Protocol
from .transitions import MERGE_TOLERANCE_HZ, TransitionList

Weight = Union[str, float]

_BACKEND = {"name": "numpy", "device": "cpu"}


def set_backend(name: str = "numpy", device: Optional[str] = None) -> dict:
    """Select the dense linear-algebra backend of exchange spectra: "numpy" or "torch" (device "cpu",
    "cuda", "cuda:1", ...). Raises if torch or the device is unavailable (no silent fallback)."""
    if name == "numpy":
        _BACKEND.update(name="numpy", device="cpu")
    elif name == "torch":
        import torch
        device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        if device.startswith("cuda") and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but torch.cuda.is_available() is False.")
        _BACKEND.update(name="torch", device=device)
    else:
        raise ValueError("backend must be 'numpy' or 'torch'.")
    _rank1_basis.cache_clear()
    return dict(_BACKEND)


def get_backend() -> dict:
    return dict(_BACKEND)


def backend_from_environment() -> None:
    """Apply ZULF_LINALG_DEVICE (numpy | cpu | cuda | cuda:N); cpu and cuda use torch."""
    import os
    value = os.environ.get("ZULF_LINALG_DEVICE", "").strip()
    if value and value != "numpy":
        set_backend("torch", value)


def _eigh(matrix: np.ndarray):
    if _BACKEND["name"] == "torch":
        import torch
        m = torch.as_tensor(matrix, device=_BACKEND["device"])
        e, v = torch.linalg.eigh(m)
        return e.cpu().numpy(), v.cpu().numpy()
    return np.linalg.eigh(matrix)



def _decompose(matrix: np.ndarray, inverse: bool):
    """Eigenvalues, right eigenvectors and (optionally) their inverse of a dense complex matrix."""
    if _BACKEND["name"] == "torch":
        import torch
        m = torch.as_tensor(matrix, dtype=torch.complex128, device=_BACKEND["device"])
        lam, vec = torch.linalg.eig(m)
        inv = torch.linalg.inv(vec) if inverse else None
        return lam.cpu().numpy(), vec.cpu().numpy(), None if inv is None else inv.cpu().numpy()
    lam, vec = sla.eig(matrix)
    return lam, vec, (np.linalg.inv(vec) if inverse else None)


def _matmul3(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """a @ b @ c on the selected backend (the costly step of the derivatives)."""
    if _BACKEND["name"] == "torch":
        import torch
        dev = _BACKEND["device"]
        t = [torch.as_tensor(x, dtype=torch.complex128, device=dev) for x in (a, b, c)]
        return (t[0] @ t[1] @ t[2]).cpu().numpy()
    return a @ b @ c


def _site_matrices(spins: Tuple[Fraction, ...]):
    """Sparse embedded (Ix, Iy, Iz) of every site and the site dimensions."""
    dims = [int(2 * s + 1) for s in spins]
    out = []
    for i, s in enumerate(spins):
        ops = angular_momentum(s)
        left = int(np.prod(dims[:i])) if i else 1
        right = int(np.prod(dims[i + 1:])) if i + 1 < len(dims) else 1
        out.append([sp.kron(sp.kron(sp.identity(left, format="csr"), sp.csr_matrix(o)),
                            sp.identity(right, format="csr"), format="csr") for o in ops])
    return out, dims


def _commutator_super(a: sp.spmatrix) -> sp.spmatrix:
    """Superoperator of rho -> [A, rho] on row-major vec(rho): A (x) 1 - 1 (x) A^T."""
    n = a.shape[0]
    eye = sp.identity(n, format="csr")
    return (sp.kron(a, eye, format="csr") - sp.kron(eye, a.T, format="csr")).tocsr()


def _replacement_super(site: int, dims) -> sp.spmatrix:
    """Superoperator of rho -> Tr_x(rho) (x) 1_x / d_x for site x: (1/d) sum_{m,m'} E_mm' (x) E_mm'."""
    d = dims[site]
    left = int(np.prod(dims[:site])) if site else 1
    right = int(np.prod(dims[site + 1:])) if site + 1 < len(dims) else 1
    total = None
    for m in range(d):
        for mp in range(d):
            e = sp.csr_matrix(([1.0], ([m], [mp])), shape=(d, d))
            full = sp.kron(sp.kron(sp.identity(left, format="csr"), e), sp.identity(right, format="csr"), format="csr")
            term = sp.kron(full, full, format="csr")
            total = term if total is None else total + term
    return (total / d).tocsr()



# projected pair superoperators are cached when the subspace is at most this large (memory m^2 per pair)
PAIR_CACHE_MAX_DIMENSION = 1600


class _Structure:
    """Spin-list data independent of couplings and rates: the rank-1 q = 0 basis B (real), projected
    single-spin I_z vectors, and projected superoperators (commutators with I_p . I_q, replacement maps)."""

    def __init__(self, spins: Tuple[Fraction, ...]):
        self.spins = spins
        self.site, self.dims = _site_matrices(spins)
        self.n = int(np.prod(self.dims))
        mz = np.real(sum(ops[2] for ops in self.site).diagonal())
        key = np.round(2 * mz).astype(int)
        a_idx, b_idx = np.nonzero(key[:, None] == key[None, :])
        self.q0 = a_idx * self.n + b_idx
        casimir = None
        for axis in range(3):
            c = _commutator_super(sp.csr_matrix(sum(ops[axis] for ops in self.site)))
            term = c @ c
            casimir = term if casimir is None else casimir + term
        cs = casimir.tocsr()[self.q0][:, self.q0].toarray()
        if np.abs(cs.imag).max(initial=0.0) < 1e-12:
            cs = cs.real
        cs = (cs + cs.conj().T) / 2
        e, v = _eigh(cs)
        self.basis = v[:, np.abs(e - 2.0) < 1e-6]
        self.m = self.basis.shape[1]
        self.iz = [self.project_vector(np.asarray(ops[2].toarray()).ravel()) for ops in self.site]
        self._pairs: Dict[Tuple[int, int], np.ndarray] = {}
        self._replace: Dict[int, np.ndarray] = {}

    def project_vector(self, vec: np.ndarray) -> np.ndarray:
        return self.basis.conj().T @ vec[self.q0]

    def project_super(self, sup: sp.spmatrix) -> np.ndarray:
        s = sup.tocsr()[self.q0][:, self.q0]
        out = self.basis.conj().T @ np.asarray(s @ self.basis)
        return out.real if np.abs(out.imag).max(initial=0.0) < 1e-12 else out

    def pair(self, p: int, q: int) -> np.ndarray:
        """B^H [I_p . I_q, .] B (without the -2 pi i factor)."""
        key = (min(p, q), max(p, q))
        if key in self._pairs:
            return self._pairs[key]
        op = sum(self.site[p][k] @ self.site[q][k] for k in range(3))
        out = self.project_super(_commutator_super(sp.csr_matrix(op)))
        if self.m <= PAIR_CACHE_MAX_DIMENSION:
            self._pairs[key] = out
        return out

    def replacement(self, x: int) -> np.ndarray:
        """B^H (Tr_x(.) (x) 1_x / d_x) B."""
        if x not in self._replace:
            self._replace[x] = self.project_super(_replacement_super(x, self.dims))
        return self._replace[x]


@lru_cache(maxsize=16)
def _structure(spins: Tuple[Fraction, ...]) -> _Structure:
    return _Structure(spins)


def _rank1_basis(spins: Tuple[Fraction, ...]):
    st = _structure(spins)
    return st.q0, st.basis, st.dims


_rank1_basis.cache_clear = _structure.cache_clear


def _weight(value: Weight, symbol: str, protocol: Protocol) -> float:
    if isinstance(value, str):
        if value == "preparation":
            return protocol.preparation_weight(symbol)
        raise ValueError("A solvent weight is 'preparation' or a number.")
    return float(value)


class _Modes:
    """Eigen-decomposition of the augmented Liouvillian [[L, source], [0, 0]] (state [delta; 1])."""

    def __init__(self, system: SpinSystem, rates_per_s: Mapping[int, float], protocol: Protocol,
                 solvent_weight, inverse: bool, spins: Optional[Tuple[Fraction, ...]] = None,
                 norm: Optional[float] = None):
        if protocol.pulses or protocol.has_field:
            raise NotImplementedError("Exchange spectra support the sudden-drop protocol without pulses or field.")
        registry = get_registry()
        # spins: site spin quantum numbers (default: the isotopes' spins; a reduced group site has its total spin)
        spins = spins or tuple(Fraction(registry.spin(s)) for s in system.isotopes)
        st = _structure(spins)
        self.st, self.system, self.protocol = st, system, protocol
        m = st.m
        j = system.couplings_hz
        liou = np.zeros((m, m), complex)
        for a in range(len(spins)):
            for b in range(a + 1, len(spins)):
                if j[a, b] != 0.0:
                    liou += (-2j * np.pi * j[a, b]) * st.pair(a, b)
        source = np.zeros(m, complex)
        self.rates, self.eps = {}, {}
        for x, k in rates_per_s.items():
            x, k = int(x), float(k)
            if k < 0 or not np.isfinite(k):
                raise ValueError("Exchange rates must be finite and nonnegative.")
            w = solvent_weight.get(x, "preparation") if isinstance(solvent_weight, Mapping) else solvent_weight
            self.eps[x] = _weight(w, system.isotopes[x], protocol)
            self.rates[x] = k
            if k == 0.0:
                continue
            liou += k * (st.replacement(x) - np.eye(m))
            source += k * self.eps[x] * st.iz[x]
        aug = np.zeros((m + 1, m + 1), complex)
        aug[:m, :m] = liou
        aug[:m, m] = source
        v0 = sum(protocol.preparation_weight(s) * st.iz[i] for i, s in enumerate(system.isotopes))
        d = sum(protocol.detection_weight(s) * st.iz[i] for i, s in enumerate(system.isotopes))
        self.lam, self.vec, self.inv = _decompose(aug, inverse)
        state = np.concatenate([v0, [1.0]]).astype(complex)
        self.g = self.inv @ state if inverse else np.linalg.solve(self.vec, state)
        self.w = np.concatenate([d, [0.0]]).conj() @ self.vec
        if norm is not None:
            self.norm = float(norm)
        else:
            self.norm = 1.0 / st.n if protocol.normalize_by_dimension else 1.0
        self.coef = self.norm * self.w * self.g          # signal = sum_j coef_j exp(lam_j t)

    def lines(self, tolerance_hz: float):
        f = self.lam.imag / (2 * np.pi)
        rate = np.maximum(-self.lam.real, 0.0)
        idx = np.flatnonzero(f > tolerance_hz)
        order = idx[np.lexsort((rate[idx], f[idx]))]
        groups, last = [], None
        for i in order:
            if last is not None and abs(f[i] - f[last]) <= tolerance_hz and \
                    abs(rate[i] - rate[last]) <= max(tolerance_hz, 1e-9 * rate[i]):
                groups[-1].append(i)
            else:
                groups.append([i])
            last = i
        return f, rate, groups


def _total_spin_multiplicities(n: int):
    """[(S, number of irreducible copies)] of n coupled spin-1/2 (S = n/2, n/2 - 1, ..., S >= 0)."""
    from math import comb
    out = []
    for k in range(n // 2 + 1):
        mult = comb(n, k) - (comb(n, k - 1) if k else 0)
        out.append((Fraction(n, 2) - k, mult))
    return out


class _Block:
    """One block of the group-reduced system: sites (each an original spin or the total spin of a
    non-exchanging equivalence group of spin-1/2), its SpinSystem over the sites, and its weight."""

    def __init__(self, system, sites, multiplicity):
        self.sites = sites                       # [(spin, members tuple)]
        self.multiplicity = multiplicity
        self.spins = tuple(sp_ for sp_, _ in sites)
        self.site_of = {m: i for i, (_, members) in enumerate(sites) for m in members}
        reps = [members[0] for _, members in sites]
        j = np.zeros((len(sites), len(sites)))
        for a in range(len(sites)):
            for b in range(a + 1, len(sites)):
                j[a, b] = j[b, a] = system.couplings_hz[reps[a], reps[b]]
        self.system = SpinSystem(tuple(system.isotopes[r] for r in reps), j, groups=[(i,) for i in range(len(sites))])


def _blocks(system: SpinSystem, exchanging) -> Optional[list]:
    """Group-reduced blocks, or None when no group can be reduced. A group is reduced when it has two or more
    spin-1/2 members and none of them exchanges: H, preparation and detection only see its total spin, so
    every total spin S (multiplicity from the coupling of n spin-1/2; S = 0 sites carry nothing) is one
    site. Exchanging spins stay individual (they exchange independently)."""
    registry = get_registry()
    options, reducible = [], False
    for g in system.groups:
        g = tuple(int(i) for i in g)
        half = all(Fraction(registry.spin(system.isotopes[i])) == Fraction(1, 2) for i in g)
        if len(g) >= 2 and half and not any(i in exchanging for i in g):
            reducible = True
            options.append([([(S, g)] if S > 0 else [], mult) for S, mult in _total_spin_multiplicities(len(g))])
        else:
            options.append([([(Fraction(registry.spin(system.isotopes[i])), (i,)) for i in g], 1)])
    if not reducible:
        return None
    import itertools
    blocks = []
    for combo in itertools.product(*options):
        sites = [site for part, _ in combo for site in part]
        mult = int(np.prod([m for _, m in combo]))
        if sites:
            blocks.append(_Block(system, sites, mult))
    return blocks


def _block_modes(system, rates_per_s, protocol, solvent_weight, inverse, reduce):
    """[(block or None, _Modes)] for the full system (one entry) or for every group-reduced block."""
    exchanging = {int(x) for x in rates_per_s}
    blocks = _blocks(system, exchanging) if reduce else None
    if blocks is None:
        return [(None, _Modes(system, rates_per_s, protocol, solvent_weight, inverse=inverse))]
    registry = get_registry()
    n_full = int(np.prod([int(2 * Fraction(registry.spin(s)) + 1) for s in system.isotopes]))
    out = []
    for b in blocks:
        rates = {b.site_of[int(x)]: k for x, k in rates_per_s.items()}
        weight = ({b.site_of[int(x)]: w for x, w in solvent_weight.items()} if isinstance(solvent_weight, Mapping)
                  else solvent_weight)
        norm = b.multiplicity / n_full if protocol.normalize_by_dimension else float(b.multiplicity)
        out.append((b, _Modes(b.system, rates, protocol, weight, inverse=inverse, spins=b.spins, norm=norm)))
    return out


def _concatenate(parts):
    """One TransitionList from the per-block lists (frequency order) and the permutation used."""
    ff = np.concatenate([t.frequencies_hz for t in parts])
    aa = np.concatenate([t.amplitudes for t in parts])
    rr = np.concatenate([t.rates_of_lines() for t in parts])
    order = np.argsort(ff, kind="stable")
    meta = dict(parts[0].metadata)
    meta.update({"method": "exchange (group-reduced)", "blocks": len(parts),
                 "subspace_dimension": [t.metadata["subspace_dimension"] for t in parts],
                 "nonoscillating_weight": float(sum(t.metadata["nonoscillating_weight"] for t in parts)),
                 "nonoscillating_rates": [r for t in parts for r in t.metadata["nonoscillating_rates"]],
                 "nonoscillating_amplitudes": [a for t in parts for a in t.metadata["nonoscillating_amplitudes"]]})
    return TransitionList(ff[order], aa[order], complex(sum(t.dc for t in parts)), metadata=meta,
                          line_rates=rr[order]), order


def exchange_transitions(system: SpinSystem, rates_per_s: Mapping[int, float], protocol: Protocol = SUDDEN_DROP,
                         solvent_weight: Union[Weight, Mapping[int, Weight]] = "preparation",
                         tolerance_hz: float = MERGE_TOLERANCE_HZ, relative_zero: float = 1e-12,
                         reduce: bool = True) -> TransitionList:
    """Transition list (with line_rates) of `system` whose spins `rates_per_s` keys exchange at those rates.

    rates_per_s: {spin index: k (1/s)}; spins not listed do not exchange. solvent_weight: the incoming spin's
    deviation weight eps ('preparation': the protocol's preparation weight of that nucleus; 0: unpolarized
    solvent), one value or one per exchanging spin. Longitudinal preparation and detection only (no pulses,
    no field). reduce: replace every non-exchanging equivalence group of spin-1/2 by its total spins (exact;
    lines of different blocks are not merged)."""
    modes = _block_modes(system, rates_per_s, protocol, solvent_weight, False, reduce)
    if modes[0][0] is None:
        return _transition_list(modes[0][1], tolerance_hz, relative_zero)[0]
    return _concatenate([_transition_list(md, tolerance_hz, relative_zero)[0] for _, md in modes])[0]


def _transition_list(md: _Modes, tolerance_hz: float, relative_zero: float):
    f, rate, groups = md.lines(tolerance_hz)
    ff = np.array([f[g[0]] for g in groups])
    rr = np.array([rate[g[0]] for g in groups])
    aa = np.array([2 * md.coef[g].sum() for g in groups], complex)
    big = np.abs(aa) > relative_zero * np.abs(aa).max() if len(aa) else np.zeros(0, bool)
    still = np.abs(f) <= tolerance_hz
    moving = rate > 1e-12 * max(1.0, float(np.abs(md.lam).max(initial=0.0)))
    dc = complex(md.coef[still & ~moving].sum())
    meta = {"method": "exchange", "protocol": md.protocol.name, "n_spins": md.system.n_spins,
            "exchange_rates_per_s": {int(x): float(k) for x, k in md.rates.items()},
            "subspace_dimension": int(md.st.m), "backend": dict(_BACKEND),
            "nonoscillating_weight": float(np.abs(md.coef[still & moving]).sum()),
            # decaying modes at zero frequency: signal term amplitude * exp(-rate t) (not rendered)
            "nonoscillating_rates": rate[still & moving].tolist(),
            "nonoscillating_amplitudes": np.real(md.coef[still & moving]).tolist()}
    tl = TransitionList(ff[big], aa[big], dc, metadata=meta, line_rates=rr[big])
    return tl, [g for g, b in zip(groups, big) if b]


class ExchangeDerivatives:
    """Line derivatives of an exchange transition list along given parameter directions.

    For direction k and line l (the list's order): d_amplitudes[k, l] is dA, t_terms[k, l] the coefficient of
    t exp(lambda t) (A d lambda for an isolated mode) and phase_terms[k, l] = A df (for the phase delay), so
    that d s / d theta = sum_l Re((dA + T t) exp(lambda_l t)) for the list convention."""

    def __init__(self, transitions, d_amplitudes, t_terms, phase_terms):
        self.transitions = transitions
        self.d_amplitudes = d_amplitudes
        self.t_terms = t_terms
        self.phase_terms = phase_terms


def exchange_transition_derivatives(system: SpinSystem, rates_per_s: Mapping[int, float],
                                    coupling_pairs: Sequence[Sequence[Tuple[int, int]]],
                                    exchange_spins: Sequence[Sequence[int]] = (),
                                    protocol: Protocol = SUDDEN_DROP, solvent_weight="preparation",
                                    tolerance_hz: float = MERGE_TOLERANCE_HZ,
                                    relative_zero: float = 1e-12, reduce: bool = True) -> ExchangeDerivatives:
    """Transition list and its derivatives with respect to (a) couplings: each direction is a list of spin
    pairs whose coupling moves together by one Hz (a group coupling), (b) log exchange rates: each direction is
    a list of spins whose log k moves together.

    With the augmented generator La = V diag(lambda) V^-1, X = V^-1 dLa V and M_ij = w_i X_ij g_j, the signal
    derivative is sum_ij M_ij F_ij(t) with F the divided difference of exp(lambda t) (Daleckii-Krein): the
    coefficient of exp(lambda_k t) moves by sum_{j not in cluster(k)} (M_kj + M_jk) / (lambda_k - lambda_j), and
    pairs inside one cluster of (numerically) equal eigenvalues give t exp(lambda t) terms.

    With reduce (group-reduced blocks, see exchange_transitions) a coupling direction must move every member
    pair of the groups it touches (a group coupling); pairs inside one reduced group have no effect."""
    modes = _block_modes(system, rates_per_s, protocol, solvent_weight, True, reduce)
    parts = []
    for block, md in modes:
        if block is None:
            parts.append(_mode_derivatives(md, coupling_pairs, exchange_spins, tolerance_hz, relative_zero))
            continue
        pairs_b = []
        for pairs in coupling_pairs:
            count: Dict[Tuple[int, int], int] = {}
            for p, q in pairs:
                a, b = block.site_of.get(int(p)), block.site_of.get(int(q))
                if a is None or b is None or a == b:
                    continue                       # an S = 0 site, or inside one reduced group
                key = (min(a, b), max(a, b))
                count[key] = count.get(key, 0) + 1
            for (a, b), c in count.items():
                need = len(block.sites[a][1]) * len(block.sites[b][1])
                if c != need:
                    raise ValueError("A coupling direction must move every member pair of a reduced group "
                                     "(reduce=False for partial directions).")
            pairs_b.append(list(count))
        spins_b = [[block.site_of[int(x)] for x in spins_] for spins_ in exchange_spins]
        parts.append(_mode_derivatives(md, pairs_b, spins_b, tolerance_hz, relative_zero))
    if len(parts) == 1:
        return parts[0]
    tl, order = _concatenate([p.transitions for p in parts])
    cat = lambda name: np.concatenate([getattr(p, name) for p in parts], axis=1)[:, order]
    return ExchangeDerivatives(tl, cat("d_amplitudes"), cat("t_terms"), cat("phase_terms"))


def _mode_derivatives(md: "_Modes", coupling_pairs, exchange_spins, tolerance_hz, relative_zero):
    tl, groups = _transition_list(md, tolerance_hz, relative_zero)
    st, m = md.st, md.st.m
    directions = []
    for pairs in coupling_pairs:
        dl = np.zeros((m + 1, m + 1), complex)
        for p, q in pairs:
            dl[:m, :m] += -2j * np.pi * st.pair(p, q)
        directions.append(dl)
    for spins_ in exchange_spins:
        dl = np.zeros((m + 1, m + 1), complex)
        for x in spins_:
            k = md.rates.get(int(x), 0.0)
            dl[:m, :m] += k * (st.replacement(int(x)) - np.eye(m))
            dl[:m, m] += k * md.eps[int(x)] * st.iz[int(x)]
        directions.append(dl)
    lam = md.lam
    scale = max(1.0, float(np.abs(lam).max(initial=0.0)))
    diff = lam[:, None] - lam[None, :]
    same = np.abs(diff) <= max(2 * np.pi * tolerance_hz, 1e-10 * scale)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_diff = np.where(same, 0.0, 1.0 / np.where(same, 1.0, diff))
    n_dir = len(directions)
    da = np.zeros((n_dir, len(tl)), complex)
    tt = np.zeros((n_dir, len(tl)), complex)
    pt = np.zeros((n_dir, len(tl)), complex)
    for k, dl in enumerate(directions):
        x = _matmul3(md.inv, dl, md.vec)
        mm = md.norm * (md.w[:, None] * x * md.g[None, :])
        dcoef = ((mm + mm.T) * inv_diff).sum(axis=1)
        tcoef = (mm * same).sum(axis=1)
        dfreq = np.diag(x).imag / (2 * np.pi)
        for l, g in enumerate(groups):
            da[k, l] = 2 * dcoef[g].sum()
            tt[k, l] = 2 * tcoef[g].sum()
            pt[k, l] = 2 * (md.coef[g] * dfreq[g]).sum()
    return ExchangeDerivatives(tl, da, tt, pt)


class ExchangeCache:
    """Exact-key cache of exchange transition lists (system content, rates, protocol, solvent weight)."""

    def __init__(self, max_entries: int = 64):
        self.max_entries = max_entries
        self._store: Dict[tuple, TransitionList] = {}

    def get(self, system: SpinSystem, rates_per_s: Mapping[int, float], protocol: Protocol = SUDDEN_DROP,
            solvent_weight: Weight = "preparation") -> TransitionList:
        key = (system.isotopes, system.couplings_hz.tobytes(), tuple(sorted((int(k), float(v))
                                                                          for k, v in rates_per_s.items())),
               repr(protocol.to_dict()), repr(solvent_weight))
        if key not in self._store:
            if len(self._store) >= self.max_entries:
                self._store.pop(next(iter(self._store)))
            self._store[key] = exchange_transitions(system, rates_per_s, protocol, solvent_weight)
        return self._store[key]


def reference_exchange_signal(system: SpinSystem, rates_per_s: Mapping[int, float], times_s: np.ndarray,
                              protocol: Protocol = SUDDEN_DROP, solvent_weight: float = None) -> np.ndarray:
    """Brute-force reference (tests only): the full Liouville space (all q, all ranks), dense superoperators
    built from explicit partial traces, affine propagation by one augmented matrix exponential per time."""
    registry = get_registry()
    spins = [Fraction(registry.spin(s)) for s in system.isotopes]
    dims = [int(2 * s + 1) for s in spins]
    n = int(np.prod(dims))
    ops = []
    for i, s in enumerate(spins):
        mats = angular_momentum(s)
        row = []
        for m in mats:
            full = np.array([[1.0]], complex)
            for k, dk in enumerate(dims):
                full = np.kron(full, m if k == i else np.eye(dk))
            row.append(full)
        ops.append(row)
    h = np.zeros((n, n), complex)
    for a in range(len(spins)):
        for b in range(a + 1, len(spins)):
            h += system.couplings_hz[a, b] * sum(ops[a][k] @ ops[b][k] for k in range(3))
    rho0 = sum(protocol.preparation_weight(s) * o[2] for s, o in zip(system.isotopes, ops))
    det = sum(protocol.detection_weight(s) * o[2] for s, o in zip(system.isotopes, ops))
    basis_ops = [np.eye(1, n * n, k).reshape(n, n) for k in range(n * n)]

    def apply(rho):
        out = -2j * np.pi * (h @ rho - rho @ h)
        for x, k in rates_per_s.items():
            left = int(np.prod(dims[:x])) if x else 1
            right = int(np.prod(dims[x + 1:])) if x + 1 < len(dims) else 1
            t = rho.reshape(left, dims[x], right, left, dims[x], right)
            reduced = np.einsum("iajkal->ijkl", t)                     # partial trace over spin x
            full = np.einsum("ijkl,ab->iajkbl", reduced, np.eye(dims[x]) / dims[x]).reshape(n, n)
            eps = protocol.preparation_weight(system.isotopes[x]) if solvent_weight is None else solvent_weight
            out = out + k * (full + eps * ops[x][2] - rho)
        return out

    # linear part as a dense matrix on vec(rho); the constant source separately
    zero = np.zeros((n, n), complex)
    src = apply(zero).reshape(-1)
    lin = np.column_stack([(apply(b) - apply(zero)).reshape(-1) for b in basis_ops])
    aug = np.zeros((n * n + 1, n * n + 1), complex)
    aug[:n * n, :n * n] = lin
    aug[:n * n, n * n] = src
    state = np.concatenate([rho0.reshape(-1), [1.0]])
    norm = 1.0 / n if protocol.normalize_by_dimension else 1.0
    out = np.empty(len(times_s))
    for i, t in enumerate(times_s):
        y = sla.expm(aug * t) @ state
        out[i] = norm * np.real(np.trace(y[:n * n].reshape(n, n) @ det))
    return out


backend_from_environment()
