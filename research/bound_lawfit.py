"""Fit |S| to a dimensionless law over the full (b-alignment x spectrum x cap) space.

No single raw scalar predicts |S| (bound_predictor.py). But the mechanism — the
support dilutes a costly objective until coancestry meets the cap — suggests a
dimensionless ratio. A concentrated plan in b's direction costs ~rayleigh = b'Gb/b'b;
diluted over |S| candidates it costs ~rayleigh/|S|; setting that to k gives the guess
|S| ~ rayleigh/k. Because the cap k itself scales with the spectrum, such a ratio may
absorb the spectral gap that defeated the raw scalars.

Tests four dimensionless ratios against |S| over a sweep that varies the spectrum,
b's alignment AND the cap, plus the real panels, and reports which collapses |S|
onto a single power law (log-log fit). Run:  python3 research/bound_lawfit.py
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

def synth_G(n, decay, seed):
    rng = np.random.default_rng(seed)
    U, _ = np.linalg.qr(rng.standard_normal((n, n)))
    lam = np.arange(1, n + 1).astype(float) ** (-decay) + 1e-3
    G = (U * lam) @ U.T
    return 0.5 * (G + G.T), U

def instance_stats(G, b, k):
    n = len(b); kmin = float(np.ones(n) @ G @ np.ones(n)) / n ** 2
    jb = int(np.argmax(b)); kg = float(G[jb, jb])                 # gain-greedy single-candidate cost
    rayleigh = float(b @ G @ b / (b @ b))
    return kmin, kg, rayleigh

records = []  # (label, |S|, rayleigh, kg, kmin, k)
n = 600
print("=== sweep: 3 spectra x 5 alignments x 3 caps (mean |S| over 3 b) ===")
for decay in [0.5, 1.0, 2.0]:
    G, U = synth_G(n, decay, 7)
    one = np.ones(n); kmin = float(one @ G @ one) / n ** 2; meand = float(np.diag(G).mean())
    for d in [1, 3, 8, 21, 55]:
        for frac in [0.2, 0.4, 0.8]:
            k = kmin + frac * (meand - kmin)
            sizes = []
            for sd in range(3):
                rng = np.random.default_rng(1000 * sd + 7 * d)
                b = U[:, :d] @ rng.standard_normal(d); b = (b - b.mean()) / b.std()
                c, S = support_first(G, b, k)
                if float(c @ G @ c) <= k * (1 + 1e-4): sizes.append(len(S))
            if sizes:
                _, kg, ray = instance_stats(G, b, k)
                records.append((f"dec{decay}/d{d}/f{frac}", float(np.mean(sizes)), ray, kg, kmin, k))

for name, nn in [("wheat", 599), ("mouse", 1814)]:
    K = np.loadtxt(f"/tmp/bench_K_{nn}.csv", delimiter=",")
    b = np.genfromtxt(f"/tmp/bench_bc_{nn}.csv", delimiter=",", skip_header=1, usecols=0).astype(float)
    ub = float(open(f"/tmp/bench_ub_{nn}.txt").read().strip())
    c, S = support_first(K, b, ub)
    kmin, kg, ray = instance_stats(K, b, ub)
    records.append((name, float(len(S)), ray, kg, kmin, ub))
    print(f"  real {name}: |S|={len(S)}  rayleigh={ray:.4f}  kg={kg:.4f}  kmin={kmin:.5f}  k={ub:.5f}")

S = np.array([r[1] for r in records])
ray = np.array([r[2] for r in records]); kg = np.array([r[3] for r in records])
kmin = np.array([r[4] for r in records]); k = np.array([r[5] for r in records])
ratios = {
    "rayleigh/k":            ray / k,
    "(rayleigh-kmin)/(k-kmin)": (ray - kmin) / np.maximum(k - kmin, 1e-12),
    "kg/k":                  kg / k,
    "(kg-kmin)/(k-kmin)":    (kg - kmin) / np.maximum(k - kmin, 1e-12),
}
print("\n=== log-log fit of |S| vs each dimensionless ratio (all instances) ===")
print(f"  {'ratio':>26} {'exponent':>9} {'R^2':>6} {'Pearson(logS,logr)':>20}")
ls = np.log(S)
for nm, x in ratios.items():
    m = (x > 0) & np.isfinite(x)
    lx = np.log(x[m]); lsi = ls[m]
    a, b0 = np.polyfit(lx, lsi, 1)
    pred = a * lx + b0; ss_res = ((lsi - pred) ** 2).sum(); ss_tot = ((lsi - lsi.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot; pear = np.corrcoef(lx, lsi)[0, 1]
    print(f"  {nm:>26} {a:>9.2f} {r2:>6.2f} {pear:>20.2f}")
print("\n(High R^2 = |S| ~ A * ratio^exponent collapses across spectra, alignments, caps")
print(" AND the two real panels. That ratio is the candidate law.)")
