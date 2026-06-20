"""The kernel identity: effective support PR(c*) = R(c*)/k.

At an optimum with the kinship cap active (c'Gc* = k) and simplex budget (sum c* = 1),
the participation ratio PR(c*) = 1/||c*||^2 (the effective number of contributors)
equals R(c*)/k, where R(c*) = c*'Gc*/c*'c* is the optimum's Rayleigh quotient — because
R*||c*||^2 = c*'Gc* = k. Hence PR <= lambda_1/k (R <= lambda_1), and the conditional
support bound reduces to bounding R(c*), small exactly when c* avoids the expensive
top directions (assumption A2).

This checks the identity (PR == R/k), the bound (PR <= lambda_1/k), and compares the
effective support PR to the raw cardinality |S|, on the real panels and synthetic
cases (G = eps*I; a gapped G with b benign vs b on the top eigenvector).

Run:  python3 research/bound_kernel.py   (needs the /tmp exports for the real panels)
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

def report(name, G, b, k, lam1):
    c, S = support_first(G, b, k)
    nrm2 = float(c @ c); pr = 1.0 / nrm2
    cGc = float(c @ G @ c); R = cGc / nrm2
    print(f"  {name:>22} |S|={len(S):>4}  PR={pr:>7.1f}  R/k={R/k:>7.1f}  R={R:>8.4f}  "
          f"λ1/k={lam1/k:>8.1f}  cGc/k={cGc/k:.3f}")

print("=== PR(c*) = R(c*)/k  (effective contributors);  PR <= λ1/k ===")
print("  (PR and R/k must match; PR <= λ1/k must hold; compare PR to raw |S|)\n")

# synthetic: G = eps*I (no structure) at a tight cap -> PR ~ n
n = 400; I = np.eye(n); rng = np.random.default_rng(0)
b = rng.standard_normal(n)
kmin = 1.0 / n
report("G=I, tight cap", I, b, 1.5 * kmin, 1.0)

# synthetic: gapped G, b benign (random) vs b on top eigenvector
U, _ = np.linalg.qr(rng.standard_normal((n, n)))
lam = np.arange(1, n + 1).astype(float) ** (-1.0) + 1e-3
G = (U * lam) @ U.T; G = 0.5 * (G + G.T)
one = np.ones(n); kmn = float(one @ G @ one) / n ** 2; meand = float(np.diag(G).mean())
k = kmn + 0.4 * (meand - kmn); lam1 = float(lam.max())
b_benign = U[:, 5:40] @ rng.standard_normal(35); b_benign = (b_benign - b_benign.mean()) / b_benign.std()
b_top = U[:, 0].copy(); b_top = (b_top - b_top.mean()) / b_top.std()
report("gapped, b benign", G, b_benign, k, lam1)
report("gapped, b=top eig", G, b_top, k, lam1)

# real panels
for nm, nn in [("wheat", 599), ("mouse", 1814)]:
    K = np.loadtxt(f"/tmp/bench_K_{nn}.csv", delimiter=",")
    bb = np.genfromtxt(f"/tmp/bench_bc_{nn}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{nn}.txt").read().strip())
    l1 = float(np.linalg.eigvalsh(K)[-1])
    report(f"real {nm}", K, bb, ub, l1)

print("\n(PR is the effective support; |S| the raw cardinality. The conditional bound")
print(" reduces to bounding R(c*): small under A2 (c* in cheap directions), ~λ1 when b")
print(" is on the top eigendirection, = the only eigenvalue for G=I.)")
