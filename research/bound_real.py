"""Real-GRM check of the support / spectrum relationship — the missing real-data
datapoint for the support-bound question.

Loads the exported real kinship matrices (wheat n=599, mouse n=1814) and asks how
the realised support |S| relates to the spectrum of K: is |S| of the order of the
effective rank (participation ratio), or far below it?

Prerequisite (writes /tmp/bench_K_*.csv, _bc_*.csv, _ub_*.txt):
    Rscript research/repro/wheat_export.R
    Rscript research/repro/mouse_export.R
Then:  python3 research/bound_real.py
Self-contained (dense-K simplex solver mirroring support_first.py); numpy/scipy.
"""
import numpy as np
from scipy.linalg import null_space

# --- dense-K simplex support-first (matvec = K c, GSS = K[S,S]) -------------
def closed_form(K, b, k, S):
    GS = K[np.ix_(S, S)]; bS = b[S]; nS = len(S)
    if nS == 1: return (np.array([1.]), bS[0], 0.) if GS[0, 0] <= k else None
    c0 = np.ones(nS) / nS; N = null_space(np.ones((1, nS)))
    Gt = N.T @ GS @ N; bt = N.T @ bS; gt = N.T @ (GS @ c0); q0 = c0 @ GS @ c0
    gi = np.linalg.solve(Gt, gt); bi = np.linalg.solve(Gt, bt); rho2 = gt @ gi - (q0 - k)
    if rho2 <= 0: return None
    cS = c0 + N @ (-gi + np.sqrt(rho2) * bi / np.sqrt(bt @ bi))
    (mu, lam), *_ = np.linalg.lstsq(np.column_stack([np.ones(nS), 2 * (GS @ cS)]), bS, rcond=None)
    return cS, mu, lam

def support_first(K, b, k, tol=1e-7, max_iter=8000):
    n = len(b); S = [int(np.argmax(b))]; c_cur = np.zeros(n); c_cur[S[0]] = 1.
    for _ in range(max_iter):
        sol = closed_form(K, b, k, S)
        if sol is None:
            Gc = K @ c_cur
            for j in np.argsort(Gc):
                if j not in S: S.append(int(j)); break
            continue
        cS, mu, lam = sol
        if cS.min() < -tol: S = [S[i] for i in range(len(S)) if cS[i] > tol]; continue
        c = np.zeros(n); c[S] = cS; c_cur = c
        r = b - mu - 2 * lam * (K @ c); r[S] = -np.inf; j = int(np.argmax(r))
        if r[j] <= tol: return c, sorted(S)
        S.append(j)
    return c_cur, sorted(S)

def spectrum(K):
    ev = np.linalg.eigvalsh(K); ev = np.clip(ev, 0, None)[::-1]
    full = int((ev > 1e-9 * ev.max()).sum())
    pr = float(ev.sum() ** 2 / (ev ** 2).sum())        # participation ratio (effective #eigs)
    stable = float((ev ** 2).sum() / ev.max() ** 2)    # stable rank
    return ev, full, pr, stable

for name, n, opti_support in [("wheat", 599, 26), ("mouse", 1814, 26)]:
    K = np.loadtxt(f"/tmp/bench_K_{n}.csv", delimiter=",")
    b = np.genfromtxt(f"/tmp/bench_bc_{n}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{n}.txt").read().strip())
    ev, full, pr, stable = spectrum(K)
    print(f"\n=== {name} (real GRM, n={n}) ===")
    print(f"  spectrum: full rank {full} | participation ratio {pr:.1f} | stable rank {stable:.1f}")
    print(f"  top eigenvalues: {', '.join(f'{x:.2f}' for x in ev[:6])} ...")
    print(f"  optiSel (sexed) support at the working cap (ub={ub:.4f}): {opti_support}")
    print(f"  {'cap k':>12} {'|S| simplex':>12} {'cKc':>10}")
    for fr in [0.5, 1.0, 2.0, 5.0]:
        k = fr * ub
        c, S = support_first(K, b, k)
        print(f"  {k:>12.5f} {len(S):>12} {float(c @ K @ c):>10.5f}")
print("\n(Question: is |S| of the order of the participation ratio, or far below it?")
print(" Either way, this is the real-data datapoint the bound theory must reproduce.)")
