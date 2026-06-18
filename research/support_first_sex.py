"""
Sex-constrained support-first OCS solver — reference Python prototype.

This is the "true OCS" variant: contributions are split by sex, so the simplex
constraint Sum(c)=1 is replaced by two equality constraints
    Sum_{i in males}   c_i = 1/2
    Sum_{i in females} c_i = 1/2
on top of c >= 0 and the kinship cap  c^T K c <= ub.

Method (exact, not a heuristic):
  - Active-set column generation on the c >= 0 face. Start from the best male and
    best female, grow the support S by the most-promising reduced-cost candidate,
    drop indices that go negative.
  - For a fixed support S the problem is "maximise a linear objective over an
    ellipsoid intersected with the two affine sex constraints". closed_form(S)
    solves it in closed form: project onto the null space of the 2x|S| sex matrix
    A_S, reduce to "max linear over a ball", and read off the multipliers
    (mu_male, mu_female, lambda) by least squares. No inner iterative solver.
  - Feasibility phase adds least-related candidates in a chunk (doubling the
    support) rather than one at a time; this was the dominant cost otherwise.

I/O: reads the benchmark export written by the R scripts
  /tmp/bench_K_{n}.csv   kinship matrix sKin (n x n)
  /tmp/bench_bc_{n}.csv  columns: bv, oc (optiSel solution), sex
  /tmp/bench_ub_{n}.txt  kinship upper bound
and prints the support-first optimum next to optiSel's for the same instance.

Usage:  python support_first_sex.py <n>
"""
import numpy as np, csv, time, sys
from scipy.linalg import null_space

n = int(sys.argv[1])
K = np.loadtxt(f"/tmp/bench_K_{n}.csv", delimiter=",")
rows = list(csv.reader(open(f"/tmp/bench_bc_{n}.csv")))[1:]
bv = np.array([float(r[0]) for r in rows])
oc_opti = np.array([float(r[1]) for r in rows])
sex = np.array([r[2].strip('"') for r in rows])
ub = float(open(f"/tmp/bench_ub_{n}.txt").read())
male = (sex == "male")
A = np.vstack([male.astype(float), (~male).astype(float)])
rhs = np.array([0.5, 0.5])

def closed_form(S):
    S = np.array(S); KS = K[np.ix_(S, S)]; bS = bv[S]; AS = A[:, S]
    if AS[0].sum() == 0 or AS[1].sum() == 0:
        return None
    c0 = AS.T @ np.linalg.solve(AS @ AS.T, rhs)
    N = null_space(AS)
    if N.shape[1] == 0:
        cS = c0
        if cS @ KS @ cS > ub + 1e-9 or cS.min() < -1e-9:
            return None
        Mx = np.column_stack([AS[0], AS[1], 2*(KS@cS)])
        (muM, muF, lam), *_ = np.linalg.lstsq(Mx, bS, rcond=None)
        return cS, muM, muF, max(lam, 0.0)
    Gt = N.T@KS@N; bt = N.T@bS; gt = N.T@(KS@c0); q0 = c0@KS@c0
    gi = np.linalg.solve(Gt, gt); bi = np.linalg.solve(Gt, bt)
    rho2 = gt@gi - (q0 - ub)
    if rho2 <= 0:
        return None
    y = -gi + np.sqrt(rho2)*bi/np.sqrt(bt@bi)
    cS = c0 + N@y
    Mx = np.column_stack([AS[0], AS[1], 2*(KS@cS)])
    (muM, muF, lam), *_ = np.linalg.lstsq(Mx, bS, rcond=None)
    return cS, muM, muF, lam

def support_first(tol=1e-7, max_iter=8000):
    bestM = int(np.where(male)[0][np.argmax(bv[male])])
    bestF = int(np.where(~male)[0][np.argmax(bv[~male])])
    S = [bestM, bestF]; c_cur = np.zeros(n); c_cur[S] = 0.5
    dropped = np.zeros(n, bool)
    for it in range(max_iter):
        sol = closed_form(S)
        if sol is None:
            # Support infeasible: add least-related candidates in a CHUNK
            # (double the support) rather than one at a time — the feasibility
            # phase was the cost (675/768 iters at n=1000). Chunking it cuts
            # iterations ~30x with the identical optimum.
            Kc = K @ c_cur
            order = [int(x) for x in np.argsort(Kc) if x not in S and not dropped[x]]
            if not order: break
            S.extend(order[:max(1, len(S))]); continue
        cS, muM, muF, lam = sol
        if cS.min() < -tol:
            keep = [S[i] for i in range(len(S)) if cS[i] > tol]
            for i in range(len(S)):
                if cS[i] <= tol: dropped[S[i]] = True
            if not any(male[i] for i in keep): keep.append(bestM)
            if not any(~male[i] for i in keep): keep.append(bestF)
            S = keep; continue
        c = np.zeros(n)
        for idx, i in enumerate(S): c[i] = cS[idx]
        c_cur = c
        Kc = K @ c
        r = bv - np.where(male, muM, muF) - 2*lam*Kc
        r[S] = -np.inf
        j = int(np.argmax(r))
        if r[j] <= tol:
            return c, sorted(S), it+1
        dropped[:] = False
        S.append(j)
    return c_cur, sorted(S), max_iter

t0 = time.perf_counter()
c, S, iters = support_first()
t_sf = time.perf_counter() - t0
gain_sf, gain_opti = bv @ c, bv @ oc_opti
print(f"  support-first n={n}: temps={t_sf:.3f}s gain={gain_sf:.5f} support={len(S)} "
      f"cKc={c@K@c:.5f}<=ub={ub:.5f}? {c@K@c<=ub+1e-6} iters={iters}")
print(f"  vs optiSel gain={gain_opti:.5f}  Δgain={gain_sf-gain_opti:+.2e}  "
      f"(sf {'≥' if gain_sf>=gain_opti-1e-9 else '<'} optiSel)")
