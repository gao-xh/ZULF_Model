"""Usage: python scripts/field_switch_check.py C2|C4

Brute-force simulation of the Ajoy-lab field sequence (PNAS Nexus 2025, pgaf187) for pyridine 13C isotopomers:
state along x at 80 uT -> GF1 ramped linearly to 0 over T_ramp with GF2 = 1 uT along -y on -> hold at 1 uT ->
sudden switch-off -> zero-field evolution, detection along y. Compare line amplitudes with the ideal sudden drop
(rho = sum gamma I_y) after removing an overall complex scale and a time shift (both absorbed by the fit)."""
import numpy as np, sys
GH, GC = 42.577, 10.7084   # Hz/uT
sx = np.array([[0, .5], [.5, 0]], complex); sy = np.array([[0, -.5j], [.5j, 0]]); sz = np.diag([.5, -.5]).astype(complex)
def ops(n):
    out = []
    for i in range(n):
        row = []
        for s in (sx, sy, sz):
            m = np.array([[1]], complex)
            for k in range(n):
                m = np.kron(m, s if k == i else np.eye(2))
            row.append(m)
        out.append(row)
    return out
def hj(J, O):
    n = len(O); H = 0
    for i in range(n):
        for k in range(i + 1, n):
            if J[i, k]:
                H = H + J[i, k] * sum(O[i][a] @ O[k][a] for a in range(3))
    return H
def hz(g, O, B):  # Hz; B in uT (Bx, By, Bz)
    return -sum(g[i] * sum(B[a] * O[i][a] for a in range(3)) for i in range(len(O)))
def lines(H, rho, D, fmin=100):
    w, V = np.linalg.eigh(H); r = V.conj().T @ rho @ V; d = V.conj().T @ D @ V
    f, a = [], []
    for i in range(len(w)):
        for k in range(len(w)):
            fr = w[k] - w[i]
            if fr > fmin:
                f.append(fr); a.append(2 * r[i, k] * d[k, i])
    f = np.array(f); a = np.array(a); o = np.argsort(f)
    f, a = f[o], a[o]
    # merge degenerate lines
    F, A = [], []
    for fi, ai in zip(f, a):
        if F and abs(fi - F[-1]) < 1e-6: A[-1] += ai
        else: F.append(fi); A.append(ai)
    return np.array(F), np.array(A)
def compare(fr, ar, a):
    taus = np.linspace(-60e-3, 60e-3, 24001)
    E = np.exp(2j * np.pi * np.outer(taus, fr)) * ar[None, :]
    c = (E.conj() @ a) / np.vdot(ar, ar)
    res = np.linalg.norm(a[None, :] - c[:, None] * E, axis=1) / np.linalg.norm(a)
    i = int(np.argmin(res))
    return res[i], taus[i], c[i]
def run(label, J, iso, T_ramp_list, holds):
    n = len(iso); O = ops(n); g = [GH if s == "H" else GC for s in iso]
    HJ = hj(J, O); Dy = sum(g[i] * O[i][1] for i in range(n))
    fr, ar = lines(HJ, sum(g[i] * O[i][1] for i in range(n)), Dy)
    keep = np.abs(ar) > 1e-3 * np.abs(ar).max(); fr, ar = fr[keep], ar[keep]
    # state at 80 uT along x: dephased Zeeman state (diagonal part in the 80 uT eigenbasis)
    H80 = HJ + hz(g, O, (80.0, -1.0, 0.0)); w, V = np.linalg.eigh(H80)
    r0 = V.conj().T @ sum(g[i] * O[i][0] for i in range(n)) @ V
    rho_start = V @ np.diag(np.diag(r0)) @ V.conj().T
    for T in T_ramp_list:
        dt = 5e-6; steps = int(round(T / dt)); rho = rho_start.copy()
        for s in range(steps):
            b1 = 80.0 * (1 - (s + 0.5) / steps)
            w, V = np.linalg.eigh(HJ + hz(g, O, (b1, -1.0, 0.0)))
            U = V @ np.diag(np.exp(-2j * np.pi * w * dt)) @ V.conj().T
            rho = U @ rho @ U.conj().T
        w1, V1 = np.linalg.eigh(HJ + hz(g, O, (0.0, -1.0, 0.0)))
        for th in holds:
            U = V1 @ np.diag(np.exp(-2j * np.pi * w1 * th)) @ V1.conj().T
            r = U @ rho @ U.conj().T
            f, a = lines(HJ, r, Dy); a = np.array([a[np.argmin(abs(f - x))] for x in fr])
            res, tau, c = compare(fr, ar, a)
            # the part of the signal that is magnetisation along y
            my = np.trace(r @ sum(g[i] * O[i][1] for i in range(n))).real / np.trace(rho_start @ sum(g[i] * O[i][0] for i in range(n))).real
            print(f"{label} ramp {T*1e3:5.1f} ms hold {th*1e3:6.1f} ms: residual vs sudden drop {res:.3f} "
                  f"(tau {tau*1e3:+.3f} ms, |scale| {abs(c):.3f}); lines {len(fr)}; M_y/M_x0 {my:.3f}", flush=True)
# pyridine ring: indices C?, H2..H6
Jl = {("H2", "H3"): 4.9, ("H2", "H4"): 1.8, ("H2", "H5"): 0.9, ("H2", "H6"): -0.1, ("H3", "H4"): 7.7,
      ("H3", "H5"): 1.4, ("H3", "H6"): 0.9, ("H4", "H5"): 7.7, ("H4", "H6"): 1.8, ("H5", "H6"): 4.9}
CH = {"C2": {"H2": 178.0, "H3": 3.0, "H4": 6.8, "H5": -1.5, "H6": 11.0},
      "C4": {"H2": 6.8, "H3": 0.9, "H4": 162.0, "H5": 0.9, "H6": 6.8}}
Hs = ["H2", "H3", "H4", "H5", "H6"]
which = sys.argv[1] if len(sys.argv) > 1 else "C2"
J = np.zeros((6, 6))
for (a, b), v in Jl.items():
    i, k = Hs.index(a) + 1, Hs.index(b) + 1; J[i, k] = J[k, i] = v
for h, v in CH[which].items():
    k = Hs.index(h) + 1; J[0, k] = J[k, 0] = v

def dephased(label, J, iso, T):
    n = len(iso); O = ops(n); g = [GH if s == "H" else GC for s in iso]
    HJ = hj(J, O); Dy = sum(g[i] * O[i][1] for i in range(n))
    fr, ar = lines(HJ, Dy, Dy)
    keep = np.abs(ar) > 1e-3 * np.abs(ar).max(); fr, ar = fr[keep], ar[keep]
    H80 = HJ + hz(g, O, (80.0, -1.0, 0.0)); w, V = np.linalg.eigh(H80)
    r0 = V.conj().T @ sum(g[i] * O[i][0] for i in range(n)) @ V
    rho = V @ np.diag(np.diag(r0)) @ V.conj().T
    dt = 5e-6; steps = int(round(T / dt))
    for s in range(steps):
        b1 = 80.0 * (1 - (s + 0.5) / steps)
        w, V = np.linalg.eigh(HJ + hz(g, O, (b1, -1.0, 0.0)))
        U = V @ np.diag(np.exp(-2j * np.pi * w * dt)) @ V.conj().T
        rho = U @ rho @ U.conj().T
    w1, V1 = np.linalg.eigh(HJ + hz(g, O, (0.0, -1.0, 0.0)))
    r1 = V1.conj().T @ rho @ V1
    # dephased in the 1 uT eigenbasis (long hold with field inhomogeneity): keep near-degenerate blocks
    mask = np.abs(w1[:, None] - w1[None, :]) < 1e-3
    r = V1 @ (r1 * mask) @ V1.conj().T
    f, a = lines(HJ, r, Dy); a = np.array([a[np.argmin(abs(f - x))] for x in fr])
    res, tau, c = compare(fr, ar, a)
    print(f"{label} ramp {T*1e3:.0f} ms, long hold (dephased at 1 uT): residual vs sudden drop {res:.3f} (tau {tau*1e3:+.3f} ms)", flush=True)

if __name__ == "__main__":
    run(f"pyridine-{which}", J, ["C"] + ["H"] * 5, [0.03], [0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3])
    dephased(f"pyridine-{which}", J, ["C"] + ["H"] * 5, 0.03)
