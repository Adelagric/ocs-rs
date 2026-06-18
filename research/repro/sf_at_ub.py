"""Sexed support-first at an arbitrary kinship bound ub (passed on argv).

Same exact solver as research/support_first_sex.py, but ub is a parameter so we
can evaluate the exact OCS frontier at exactly the coancestry AlphaMate reaches,
for an equal-diversity head-to-head. Reports gain, achieved c'Kc, support, time.

Usage:  python3 sf_at_ub.py <n> <ub> [<ub> ...]
"""
import numpy as np, csv, time, sys
from scipy.linalg import null_space

n = int(sys.argv[1])
K = np.loadtxt(f"/tmp/bench_K_{n}.csv", delimiter=",")
rows = list(csv.reader(open(f"/tmp/bench_bc_{n}.csv")))[1:]
bv = np.array([float(r[0]) for r in rows])
sex = np.array([r[2].strip('"') for r in rows])
male = (sex == "male")
A = np.vstack([male.astype(float), (~male).astype(float)])
rhs = np.array([0.5, 0.5])


def closed_form(S, ub):
    S = np.array(S); KS = K[np.ix_(S, S)]; bS = bv[S]; AS = A[:, S]
    if AS[0].sum() == 0 or AS[1].sum() == 0:
        return None
    c0 = AS.T @ np.linalg.solve(AS @ AS.T, rhs)
    N = null_space(AS)
    if N.shape[1] == 0:
        cS = c0
        if cS @ KS @ cS > ub + 1e-9 or cS.min() < -1e-9:
            return None
        Mx = np.column_stack([AS[0], AS[1], 2 * (KS @ cS)])
        (muM, muF, lam), *_ = np.linalg.lstsq(Mx, bS, rcond=None)
        return cS, muM, muF, max(lam, 0.0)
    Gt = N.T @ KS @ N; bt = N.T @ bS; gt = N.T @ (KS @ c0); q0 = c0 @ KS @ c0
    gi = np.linalg.solve(Gt, gt); bi = np.linalg.solve(Gt, bt)
    rho2 = gt @ gi - (q0 - ub)
    if rho2 <= 0:
        return None
    y = -gi + np.sqrt(rho2) * bi / np.sqrt(bt @ bi)
    cS = c0 + N @ y
    Mx = np.column_stack([AS[0], AS[1], 2 * (KS @ cS)])
    (muM, muF, lam), *_ = np.linalg.lstsq(Mx, bS, rcond=None)
    return cS, muM, muF, lam


def support_first(ub, tol=1e-7, max_iter=8000):
    bestM = int(np.where(male)[0][np.argmax(bv[male])])
    bestF = int(np.where(~male)[0][np.argmax(bv[~male])])
    S = [bestM, bestF]; c_cur = np.zeros(n); c_cur[S] = 0.5
    dropped = np.zeros(n, bool)
    for it in range(max_iter):
        sol = closed_form(S, ub)
        if sol is None:
            Kc = K @ c_cur
            order = [int(x) for x in np.argsort(Kc) if x not in S and not dropped[x]]
            if not order:
                break
            S.extend(order[:max(1, len(S))]); continue
        cS, muM, muF, lam = sol
        if cS.min() < -tol:
            keep = [S[i] for i in range(len(S)) if cS[i] > tol]
            for i in range(len(S)):
                if cS[i] <= tol:
                    dropped[S[i]] = True
            if not any(male[i] for i in keep):
                keep.append(bestM)
            if not any(~male[i] for i in keep):
                keep.append(bestF)
            S = keep; continue
        c = np.zeros(n)
        for idx, i in enumerate(S):
            c[i] = cS[idx]
        c_cur = c
        Kc = K @ c
        r = bv - np.where(male, muM, muF) - 2 * lam * Kc
        r[S] = -np.inf
        j = int(np.argmax(r))
        if r[j] <= tol:
            return c, sorted(S), it + 1
        dropped[:] = False
        S.append(j)
    return c_cur, sorted(S), max_iter


for ub in [float(x) for x in sys.argv[2:]]:
    t0 = time.perf_counter()
    c, S, iters = support_first(ub)
    dt = time.perf_counter() - t0
    print(f"  support-first @ ub={ub:.6f}: gain={bv @ c:+.5f}  c'Kc={c @ K @ c:.6f}  "
          f"support={len(S)}  iters={iters}  temps={dt:.3f}s")
