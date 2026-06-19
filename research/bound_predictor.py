"""Search for a scalar that predicts the OCS support |S| ACROSS regimes.

The b-alignment finding (bound_balign.py) says |S| is driven by where b sits in G's
spectrum, not by G's spectrum alone. This tests candidate predictors that couple b
and the spectrum, over a sweep that varies BOTH the spectrum and b's alignment, plus
the real wheat/mouse panels. The winner (if any) is the candidate law; the spectrum-
only proxy PR(G) is included as the control known to fail.

Candidates for each instance (G, b) at a fixed binding cap:
  rayleigh   = b'Gb / b'b              directional coancestry cost of b
  spread_b   = PR(U'b)                 how many eigendirections b excites
  pr_ginvb   = PR of (G^{-1} b)_+      support of the unconstrained optimal direction
  pr_G       = PR(spectrum of G)       spectrum-only control

Run:  python3 research/bound_predictor.py   (numpy/scipy)
"""
import numpy as np
from scipy.linalg import null_space
from scipy.stats import spearmanr, pearsonr

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

def pr(x):  # participation ratio of a nonneg weight vector
    x = np.asarray(x, float); s2 = (x ** 2).sum()
    return float(s2 ** 2 / (x ** 4).sum()) if s2 > 0 else 0.0

def predictors(G, b, lam=None, U=None):
    if lam is None:
        lam, U = np.linalg.eigh(G); lam = np.clip(lam, 0, None)
    beta = U.T @ b
    g = np.linalg.solve(G, b)
    return dict(rayleigh=float(b @ G @ b / (b @ b)),
                spread_b=pr(beta),
                pr_ginvb=pr(np.clip(g, 0, None)),
                pr_G=pr(lam))

def cap(G, frac=0.4):
    n = G.shape[0]; one = np.ones(n)
    kmin = float(one @ G @ one) / n ** 2
    return kmin + frac * (float(np.diag(G).mean()) - kmin)

def synth_G(n, decay, seed):
    rng = np.random.default_rng(seed)
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    lam = np.arange(1, n + 1).astype(float) ** (-decay) + 1e-3
    G = (U * lam) @ U.T
    return 0.5 * (G + G.T), lam, U

rows = []  # (label, |S|, rayleigh, spread_b, pr_ginvb, pr_G)
print("=== synthetic: 3 spectra x 5 b-alignments, fixed binding cap ===")
print(f"  {'spectrum':>10} {'d':>3} {'|S|':>4} {'rayleigh':>9} {'spread_b':>9} {'pr_ginvb':>9} {'pr_G':>7}")
n = 600
for decay, tag in [(0.5, "flat"), (1.0, "med"), (2.0, "steep")]:
    G, lam, U = synth_G(n, decay, 7)
    k = cap(G)
    for d in [1, 3, 8, 21, 55]:
        rng = np.random.default_rng(100 + d)
        b = U[:, :d] @ rng.standard_normal(d); b = (b - b.mean()) / b.std()
        c, S = support_first(G, b, k)
        P = predictors(G, b, lam, U)
        rows.append((f"{tag}/d{d}", len(S), P["rayleigh"], P["spread_b"], P["pr_ginvb"], P["pr_G"]))
        print(f"  {tag:>10} {d:>3} {len(S):>4} {P['rayleigh']:>9.3f} {P['spread_b']:>9.1f} {P['pr_ginvb']:>9.1f} {P['pr_G']:>7.1f}")

print("\n=== real panels (working cap from the exports) ===")
for name, nn in [("wheat", 599), ("mouse", 1814)]:
    K = np.loadtxt(f"/tmp/bench_K_{nn}.csv", delimiter=",")
    b = np.genfromtxt(f"/tmp/bench_bc_{nn}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{nn}.txt").read().strip())
    c, S = support_first(K, b, ub)
    P = predictors(K, b)
    rows.append((name, len(S), P["rayleigh"], P["spread_b"], P["pr_ginvb"], P["pr_G"]))
    print(f"  {name:>10} |S|={len(S):>3}  rayleigh={P['rayleigh']:.3f}  spread_b={P['spread_b']:.1f}"
          f"  pr_ginvb={P['pr_ginvb']:.1f}  pr_G={P['pr_G']:.1f}")

S = np.array([r[1] for r in rows], float)
print("\n=== correlation with |S| (all instances, synthetic + real) ===")
for j, nm in [(2, "rayleigh"), (3, "spread_b"), (4, "pr_ginvb"), (5, "pr_G")]:
    x = np.array([r[j] for r in rows], float)
    print(f"  {nm:>10}: Pearson r={pearsonr(x, S)[0]:+.2f}  Spearman rho={spearmanr(x, S)[0]:+.2f}")
print("\n(A predictor that tracks |S| ACROSS regimes — high |rho| over both synthetic")
print(" spectra/alignments AND the real panels — is the candidate law. pr_G is the control.)")
