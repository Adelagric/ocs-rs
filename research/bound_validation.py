"""
Validation of the support-size bound for genomic OCS.

THEOREM (epsilon = 0, G0 = ZZ'/s of rank r, q equality rows): there exists an
optimal contribution vector with  |support| <= q + r + 1,  independent of n.
Proof: by strong duality the optimum maximises the Lagrangian
  L(c) = b'c - lambda* c'G0 c
over {Ac=d, c>=0}; L depends on c only through (Z'c, b'c) in R^{r+1}. The slice
  P = {c>=0 : Ac=d, Z'c=Z'c*, b'c=b'c*}
is therefore entirely optimal (same L, feasible, same objective), and is an LP
polytope with q+r+1 equality rows, so it has a vertex with <= q+r+1 nonzeros.
An active-set solver (support-first) returns such a vertex.

This script (simplex form, q=1; self-contained, mirrors support_first.py):
  (1) sanity-checks the solver against scipy SLSQP,
  (2) verifies |S| <= q + r + 1 over an EXACT low-rank sweep (noise-free, rank=r),
  (3) confirms n-independence at fixed (r, m),
  (4) tests whether |S| tracks the EFFECTIVE rank (full-rank G0, few dominant
      eigenvalues) rather than the proven worst case r.

Run:  python3 research/bound_validation.py
"""
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize

RIDGE = 1e-5
Q = 1  # simplex form: one equality row 1'c = 1

# --- solver core (verbatim from support_first.py) ---------------------------
def matvec(Z, s, c): return RIDGE * c + (Z @ (Z.T @ c)) / s
def GSS(Z, s, S): ZS = Z[S]; return ZS @ ZS.T / s + RIDGE * np.eye(len(S))

def closed_form(Z, s, b, k, S):
    GS = GSS(Z, s, S); bS = b[S]; nS = len(S)
    if nS == 1: return (np.array([1.]), bS[0], 0.) if GS[0, 0] <= k else None
    c0 = np.ones(nS) / nS; N = null_space(np.ones((1, nS)))
    Gt = N.T @ GS @ N; bt = N.T @ bS; gt = N.T @ (GS @ c0); q0 = c0 @ GS @ c0
    gi = np.linalg.solve(Gt, gt); bi = np.linalg.solve(Gt, bt); rho2 = gt @ gi - (q0 - k)
    if rho2 <= 0: return None
    cS = c0 + N @ (-gi + np.sqrt(rho2) * bi / np.sqrt(bt @ bi))
    (mu, lam), *_ = np.linalg.lstsq(np.column_stack([np.ones(nS), 2 * (GS @ cS)]), bS, rcond=None)
    return cS, mu, lam

def support_first(Z, s, b, k, tol=1e-7, max_iter=6000):
    n = len(b); nprod = 0; S = [int(np.argmax(b))]; c_cur = np.zeros(n); c_cur[S[0]] = 1.
    for it in range(max_iter):
        sol = closed_form(Z, s, b, k, S)
        if sol is None:
            Gc = matvec(Z, s, c_cur); nprod += 1
            for j in np.argsort(Gc):
                if j not in S: S.append(int(j)); break
            continue
        cS, mu, lam = sol
        if cS.min() < -tol: S = [S[i] for i in range(len(S)) if cS[i] > tol]; continue
        c = np.zeros(n); c[S] = cS; c_cur = c; Gc = matvec(Z, s, c); nprod += 1
        r = b - mu - 2 * lam * Gc; r[S] = -np.inf; j = int(np.argmax(r))
        if r[j] <= tol: return c, sorted(S), it + 1, nprod
        S.append(j)
    return c_cur, sorted(S), max_iter, nprod

# --- low-rank-controlled genotype clouds ------------------------------------
def latent_pop(n, m, n_factors, seed, noise=0.0):
    """A genotype cloud driven by `n_factors` latent factors. With noise=0,
    rank(G0) = n_factors exactly (exact low rank). With noise>0, G0 is full rank
    but has only ~n_factors dominant eigenvalues (low EFFECTIVE rank)."""
    rng = np.random.default_rng(seed)
    F = rng.standard_normal((n, n_factors))          # individual factor scores
    L = rng.standard_normal((n_factors, m))          # marker loadings
    raw = F @ L
    if noise: raw = raw + noise * rng.standard_normal((n, m))
    Z = raw - raw.mean(0)                            # centre (VanRaden-style)
    s = (Z * Z).sum() / n                            # any positive scale; k tracks it
    return Z, s

def decay_pop(n, m, K, alpha, seed, floor=0.3):
    """K latent factors with strengths sigma_k = k^{-alpha}: alpha tunes the
    spectral decay (0 = flat = high effective rank; large = steep = low effective
    rank), a small noise floor keeps G0 full rank."""
    rng = np.random.default_rng(seed)
    sig = (np.arange(1, K + 1).astype(float)) ** (-alpha)
    F = rng.standard_normal((n, K)) * sig            # scale each factor
    raw = F @ rng.standard_normal((K, m)) + floor * rng.standard_normal((n, m))
    Z = raw - raw.mean(0); s = (Z * Z).sum() / n
    return Z, s

def ranks(Z, s):
    n, m = Z.shape
    Gram = (Z @ Z.T) / s if n <= m else (Z.T @ Z) / s   # smaller matrix, same nonzero spectrum
    ev = np.linalg.eigvalsh(Gram)
    ev = np.clip(ev, 0, None)[::-1]
    full = int((ev > 1e-9 * ev.max()).sum())                 # numerical rank
    stable = float((ev ** 2).sum() / ev.max() ** 2)          # stable rank ||.||_F^2/||.||_2^2
    pr = float(ev.sum() ** 2 / (ev ** 2).sum())              # participation ratio (eff. #eigs)
    return full, stable, pr

def krange(Z, s):
    n = Z.shape[0]; d = RIDGE + (Z * Z).sum(1) / s
    one = np.ones(n); kmin = (one @ matvec(Z, s, one)) / n ** 2
    return kmin, d.mean()

def ebv(Z, seed):
    rng = np.random.default_rng(seed + 999)
    g = Z @ rng.standard_normal(Z.shape[1]); return (g - g.mean()) / g.std()

FRAC = 0.25  # binding cap (quadratic active) -> theorem regime

# === (1) sanity: solver vs scipy SLSQP =====================================
print("=== (1) sanity: support-first vs scipy (n=300, exact rank 8) ===")
Z, s = latent_pop(300, 1500, 8, 1); b = ebv(Z, 1); G = Z @ Z.T / s + RIDGE * np.eye(300)
kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
c, S, it, npd = support_first(Z, s, b, k)
cons = [{"type": "eq", "fun": lambda c: c.sum() - 1}, {"type": "ineq", "fun": lambda c: k - c @ G @ c}]
rr = minimize(lambda c: -(b @ c), np.ones(300) / 300, bounds=[(0, None)] * 300,
              constraints=cons, method="SLSQP", options={"ftol": 1e-11, "maxiter": 3000})
print(f"  |S|={len(S)}  gain_sf={b@c:.6f}  gain_scipy={-rr.fun:.6f}  Δgain={abs(b@c+rr.fun):.1e}")

# === (2) the proven bound, exact low rank ===================================
print("\n=== (2) |S| <= q + rank(G0) + 1  (q=1, noise-free exact rank, n=800, m=3000) ===")
print(f"  {'r set':>6} {'rank(G0)':>9} {'q+r+1':>7} {'|S|':>5} {'holds':>6} {'gain':>9}")
for r_set in [2, 5, 10, 20, 40, 80]:
    Z, s = latent_pop(800, 3000, r_set, 7)
    b = ebv(Z, 7); kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
    c, S, it, npd = support_first(Z, s, b, k)
    full, _, _ = ranks(Z, s); bound = Q + full + 1
    ok = "yes" if len(S) <= bound else "NO!"
    conv = "" if it < 6000 else " (maxiter)"
    print(f"  {r_set:>6} {full:>9} {bound:>7} {len(S):>5} {ok:>6} {b@c:>9.4f}{conv}")

# === (3) n-independence at fixed (rank, m) ==================================
print("\n=== (3) n-independence: exact rank r=10, m=2000, vary n ===")
print(f"  {'n':>6} {'rank(G0)':>9} {'q+r+1':>7} {'|S|':>5} {'products':>9}")
for n in [500, 1000, 2000, 5000, 10000]:
    Z, s = latent_pop(n, 2000, 10, 5)
    b = ebv(Z, 5); kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
    c, S, it, npd = support_first(Z, s, b, k)
    full, _, _ = ranks(Z, s)
    print(f"  {n:>6} {full:>9} {Q+full+1:>7} {len(S):>5} {npd:>9}")

# === (4) effective-rank conjecture: full rank, few dominant eigenvalues =====
print("\n=== (4) |S| vs EFFECTIVE rank (n=800, m=3000, full-rank G0 + noise) ===")
print(f"  {'factors':>8} {'rank(G0)':>9} {'stable_rk':>10} {'particip':>9} {'q+r+1':>7} {'|S|':>5}")
for nf in [3, 6, 12, 25, 50]:
    Z, s = latent_pop(800, 3000, nf, 11, noise=0.5)
    b = ebv(Z, 11); kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
    c, S, it, npd = support_first(Z, s, b, k)
    full, stable, pr = ranks(Z, s)
    print(f"  {nf:>8} {full:>9} {stable:>10.1f} {pr:>9.1f} {Q+full+1:>7} {len(S):>5}")
print("  (full rank ~ min(m,n): the proven bound is loose; |S| tracks the few")
print("   dominant eigenvalues = effective rank, not the full rank.)")

# === (5) ridge selects the DENSE optimum: |S| vs ridge at exact low rank ====
# The eps=0 theorem guarantees a sparse optimum (<= q+r+1 = 12 here). Does the
# ridged solver approach it as eps->0, or does the ridge select a dense one?
print("\n=== (5) ridge-driven spreading: exact rank r=10, n=1500, vary ridge eps ===")
print(f"  {'eps':>8} {'|S|':>5} {'gain':>9}   (eps=0 theorem allows |S| <= q+r+1 = 12)")
for eps in [1e-2, 1e-4, 1e-6, 1e-8]:
    RIDGE = eps   # module global; matvec/GSS/closed_form read it at call time
    Z, s = latent_pop(1500, 2000, 10, 5)
    b = ebv(Z, 5); kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
    c, S, it, npd = support_first(Z, s, b, k)
    print(f"  {eps:>8.0e} {len(S):>5} {b@c:>9.4f}")
print("  (if |S| stays large as eps->0, the ridge selects the dense min-norm")
print("   optimum: the eps=0 sparse optimum is not the eps->0 limit.)")

# === (6) the law: |S| vs effective rank across spectral-decay profiles ======
print("\n=== (6) |S| vs effective rank, sweeping spectral decay (n=800, m=3000, K=60) ===")
print(f"  {'alpha':>6} {'particip':>9} {'stable_rk':>10} {'full_rk':>8} {'|S|':>5}")
RIDGE = 1e-5  # reset after the eps sweep above
for alpha in [0.0, 0.5, 1.0, 2.0, 4.0]:
    Z, s = decay_pop(800, 3000, 60, alpha, 23)
    b = ebv(Z, 23); kmin, kmax = krange(Z, s); k = kmin + FRAC * (kmax - kmin)
    c, S, it, npd = support_first(Z, s, b, k)
    full, stable, pr = ranks(Z, s)
    print(f"  {alpha:>6.1f} {pr:>9.1f} {stable:>10.1f} {full:>8} {len(S):>5}")
print("  (alpha 0=flat .. 4=steep. |S| rises and falls with the effective rank")
print("   (participation ratio), NOT with full rank ~min(m,n)=800 — the law a")
print("   theorem for the ridged problem must reproduce.)")
