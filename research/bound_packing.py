"""c>=0-aware attack: is the optimal support a low-relatedness (packing) set?

The nonnegativity constraint is the irreducible core (bound_kernel.py). A genetically
natural mechanism: under the diversity cap the optimum selects candidates that are
mutually LITTLE related — a low-relatedness set — so |S| is bounded by a packing
number of the population (how many mutually-near-unrelated individuals it admits at
the cap level), n-independent when structure limits that, and = n for G = eps*I
(everyone unrelated -> packing = n -> Theorem 2).

Test: compare the within-support relatedness of the optimum to random sets of the
same size, on the real wheat/mouse kinship matrices.

Result (negative). The optimum is mildly diversity-favoring — within-support relatedness sits
below the gain-greedy (top-b) set (wheat -0.007 vs +0.015; mouse -0.005 vs +0.008) — but it is
NOT a low-relatedness independent set: it tolerates related high-merit pairs (max within-support
kinship 0.54 / 0.13) and |S| does not track a packing number (vacuous: 575 / 113 >> 25 / 17). So
the packing characterisation fails too. Across the natural candidates — spectral, relaxation,
packing — none gives a clean structural shortcut; the support cardinality is the genuine
combinatorial core (cardinality of a nonnegativity-constrained QP optimum).

Run:  python3 research/bound_packing.py   (needs the /tmp exports)
"""
import numpy as np
from scipy.linalg import null_space

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

def support_first(G, b, k, tol=1e-7, max_iter=12000):
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

def offdiag_mean(K, idx):
    sub = K[np.ix_(idx, idx)]
    m = len(idx)
    return (sub.sum() - np.trace(sub)) / (m * (m - 1))

def greedy_packing(K, tau, order):
    """Largest set, taken greedily in `order`, with all pairwise kinship <= tau."""
    chosen = []
    for j in order:
        if all(K[j, i] <= tau for i in chosen):
            chosen.append(j)
    return len(chosen)

rng = np.random.default_rng(0)
for name, nn in [("wheat", 599), ("mouse", 1814)]:
    K = np.loadtxt(f"/tmp/bench_K_{nn}.csv", delimiter=",")
    b = np.genfromtxt(f"/tmp/bench_bc_{nn}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{nn}.txt").read().strip())
    c, S = support_first(K, b, ub)
    S = np.array(S)
    within = offdiag_mean(K, S)
    maxin = (K[np.ix_(S, S)] - np.diag(np.diag(K[np.ix_(S, S)]))).max()
    # random same-size subsets
    rnd = np.array([offdiag_mean(K, rng.choice(nn, len(S), replace=False)) for _ in range(300)])
    # top-|S|-by-b subset (gain-greedy, ignores relatedness) for contrast
    topb = np.argsort(b)[::-1][:len(S)]
    within_topb = offdiag_mean(K, topb)
    # packing number at the support's own max within-relatedness, taken in b-order
    pack = greedy_packing(K, maxin, np.argsort(b)[::-1])
    print(f"\n=== {name} (n={nn}) ===  |S|={len(S)}")
    print(f"  mean off-diag kinship WITHIN support : {within:+.4f}   (max within: {maxin:+.4f})")
    print(f"  mean off-diag over random |S|-subsets: {rnd.mean():+.4f} ± {rnd.std():.4f}")
    print(f"  mean off-diag within top-|S|-by-b set: {within_topb:+.4f}  (gain-greedy, ignores kinship)")
    print(f"  greedy packing #(pairwise <= max within) in b-order: {pack}   (vs |S|={len(S)})")
print("\n(If WITHIN-support relatedness << random and << gain-greedy, the optimum picks a")
print(" low-relatedness set; |S| then tracks a packing number — n-independent under structure,")
print(" = n for G=εI where every pair is unrelated. That is the c>=0-aware mechanism.)")
