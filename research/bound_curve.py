"""Characterise the support-vs-cap curve |S|(k), per instance.

A single |S| value resists prediction; the CURVE may be cleaner. Hand-checking the
real-panel data hinted |S|*k ~ const, i.e. |S| ~ k^{-1}. This fits |S| ~ k^{-alpha}
per instance over a k-sweep, across synthetic spectra x b-alignments and the real
panels, and asks: is the exponent alpha near-universal (the curve SHAPE is clean),
leaving only the prefactor C = |S|*k^alpha to carry the joint (spectrum, b) scale?

Run:  python3 research/bound_curve.py   (numpy/scipy; needs the /tmp exports for real)
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

def synth_G(n, decay, seed):
    rng = np.random.default_rng(seed)
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    lam = np.arange(1, n + 1).astype(float) ** (-decay) + 1e-3
    G = (U * lam) @ U.T
    return 0.5 * (G + G.T), U

def fit_curve(G, b, ks):
    xs, ys = [], []
    for k in ks:
        c, S = support_first(G, b, k)
        if float(c @ G @ c) <= k * (1 + 1e-4) and len(S) >= 2:
            xs.append(k); ys.append(len(S))
    lx, ly = np.log(xs), np.log(ys)
    slope, intc = np.polyfit(lx, ly, 1)
    pred = slope * lx + intc
    r2 = 1 - ((ly - pred) ** 2).sum() / max(((ly - ly.mean()) ** 2).sum(), 1e-12)
    return -slope, r2, np.exp(intc), list(zip(xs, ys))   # alpha, R^2, prefactor, points

n = 600
print("=== |S| ~ k^(-alpha) per instance ===")
print(f"  {'instance':>16} {'alpha':>6} {'R^2':>6} {'C=|S|.k^a':>10}  |S| across the k-sweep")
for decay, tag in [(0.5, "flat"), (1.0, "med"), (2.0, "steep")]:
    G, U = synth_G(n, decay, 7)
    one = np.ones(n); kmin = float(one @ G @ one) / n ** 2; meand = float(np.diag(G).mean())
    ks = [kmin + f * (meand - kmin) for f in [0.2, 0.3, 0.45, 0.6, 0.8]]
    for d in [1, 8, 55]:
        rng = np.random.default_rng(7 * d)
        b = U[:, :d] @ rng.standard_normal(d); b = (b - b.mean()) / b.std()
        a, r2, C, pts = fit_curve(G, b, ks)
        sizes = " ".join(f"{s:>3}" for _, s in pts)
        print(f"  {tag + '/d' + str(d):>16} {a:>6.2f} {r2:>6.3f} {C:>10.3f}  {sizes}")

print("  --- real panels ---")
for name, nn in [("wheat", 599), ("mouse", 1814)]:
    K = np.loadtxt(f"/tmp/bench_K_{nn}.csv", delimiter=",")
    b = np.genfromtxt(f"/tmp/bench_bc_{nn}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{nn}.txt").read().strip())
    ks = [f * ub for f in [0.5, 0.8, 1.0, 1.5, 2.5, 4.0]]
    a, r2, C, pts = fit_curve(K, b, ks)
    sizes = " ".join(f"{s:>3}" for _, s in pts)
    print(f"  {name:>16} {a:>6.2f} {r2:>6.3f} {C:>10.3f}  {sizes}")
print("\n(If alpha clusters (near 1) with high R^2, the curve SHAPE |S| ~ k^-alpha is")
print(" universal; the prefactor C still carries the joint (spectrum, b) scale.)")
