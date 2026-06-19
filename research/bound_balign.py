"""Isolate the effect of b's alignment on the OCS support |S|.

Holds the relationship matrix G (hence its full spectrum) AND the coancestry cap k
fixed, and varies only how many of G's dominant eigendirections the objective b
excites. Question: does |S| track that number (b's spectral spread), explaining why
two panels of very different effective rank can show a similar support?

Construction: G = U diag(lambda) U^T with a fixed decaying spectrum and a fixed
random orthonormal basis U. For each d, b lies exactly in span(u_1..u_d) (the top-d
eigendirections). Everything else identical across d.

Run:  python3 research/bound_balign.py    (self-contained; numpy/scipy)
"""
import numpy as np
from scipy.linalg import null_space

# --- dense-G simplex support-first (matvec = G c, GSS = G[S,S]) -------------
def closed_form(G, b, k, S):
    GS = G[np.ix_(S, S)]; bS = b[S]; nS = len(S)
    if nS == 1: return (np.array([1.]), bS[0], 0.) if GS[0, 0] <= k else None
    c0 = np.ones(nS) / nS; N = null_space(np.ones((1, nS)))
    Gt = N.T @ GS @ N; bt = N.T @ bS; gt = N.T @ (GS @ c0); q0 = c0 @ GS @ c0
    gi = np.linalg.solve(Gt, gt); bi = np.linalg.solve(Gt, bt); rho2 = gt @ gi - (q0 - k)
    if rho2 <= 0: return None
    cS = c0 + N @ (-gi + np.sqrt(rho2) * bi / np.sqrt(bt @ bi))
    (mu, lam), *_ = np.linalg.lstsq(np.column_stack([np.ones(nS), 2 * (GS @ cS)]), bS, rcond=None)
    return cS, mu, lam

def support_first(G, b, k, tol=1e-7, max_iter=8000):
    n = len(b); S = [int(np.argmax(b))]; c_cur = np.zeros(n); c_cur[S[0]] = 1.
    for _ in range(max_iter):
        sol = closed_form(G, b, k, S)
        if sol is None:
            Gc = G @ c_cur
            for j in np.argsort(Gc):
                if j not in S: S.append(int(j)); break
            continue
        cS, mu, lam = sol
        if cS.min() < -tol: S = [S[i] for i in range(len(S)) if cS[i] > tol]; continue
        c = np.zeros(n); c[S] = cS; c_cur = c
        r = b - mu - 2 * lam * (G @ c); r[S] = -np.inf; j = int(np.argmax(r))
        if r[j] <= tol: return c, sorted(S)
        S.append(j)
    return c_cur, sorted(S)

n = 800
rng = np.random.default_rng(0)
# fixed orthonormal eigenbasis and a fixed decaying spectrum (+ floor -> PD)
U, _ = np.linalg.qr(rng.standard_normal((n, n)))
lam = 1.0 / np.arange(1, n + 1) + 1e-3              # harmonic decay, fixed for all d
G = (U * lam) @ U.T
G = 0.5 * (G + G.T)

# fixed binding cap in the SMALL-support regime: a fraction of the way from the
# uniform-plan minimum coancestry up to a typical single-candidate coancestry.
one = np.ones(n); kmin = float(one @ G @ one) / n ** 2
meand = float(np.diag(G).mean())
k = kmin + 0.4 * (meand - kmin)
pr_spectrum = float(lam.sum() ** 2 / (lam ** 2).sum())
print(f"n={n}  fixed spectrum: participation ratio {pr_spectrum:.1f}, lambda_1={lam[0]:.3f}")
print(f"fixed cap k = {k:.5f}  (kmin={kmin:.5f}, mean diag={meand:.5f})\n")
print(f"  {'d (top eig-dirs)':>16} {'b spread':>9} {'mean|S|':>8} {'range':>9}  (10 random b in span(top-d))")
for d in [1, 2, 3, 5, 8, 13, 21, 34, 55]:
    sizes, spreads = [], []
    for _ in range(10):
        a = rng.standard_normal(d)
        b = U[:, :d] @ a                            # b lies in span(u_1..u_d)
        b = (b - b.mean()) / b.std()
        coeff = U.T @ b
        spreads.append((coeff ** 2).sum() ** 2 / (coeff ** 4).sum())
        c, S = support_first(G, b, k)
        assert float(c @ G @ c) <= k * (1 + 1e-4), "infeasible (max_iter)"
        sizes.append(len(S))
    print(f"  {d:>16} {np.mean(spreads):>9.1f} {np.mean(sizes):>8.1f} {f'{min(sizes)}-{max(sizes)}':>9}")
print("\n(G and k are identical across rows; only b's eigen-spread changes.")
print(" If |S| tracks the spread, b's alignment is the factor the spectrum alone misses.)")
