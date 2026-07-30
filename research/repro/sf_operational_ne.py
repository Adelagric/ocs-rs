"""Given-G active-set (support-first handed dense K) at operational Ne, generic over
the panel size n. Given K, pricing is K@c = O(n^2), independent of the marker count,
so this isolates the algorithmic advantage from the matrix-free O(n*m) penalty.

Needs /tmp/bench_K_{n}.csv and /tmp/bench_bc_{n}.csv (export_K_panels.R / the pig
export). Usage: python3 research/repro/sf_operational_ne.py <n> [<n> ...]
"""
import csv
import sys
import time

import numpy as np
from scipy.linalg import null_space


def run(n):
    K = np.loadtxt(f"/tmp/bench_K_{n}.csv", delimiter=",")
    rows = list(csv.reader(open(f"/tmp/bench_bc_{n}.csv")))[1:]
    bv = np.array([float(r[0]) for r in rows])
    sex = np.array([r[2].strip('"') for r in rows])
    male = sex == "male"
    A = np.vstack([male.astype(float), (~male).astype(float)])
    rhs = np.array([0.5, 0.5])

    def closed_form(S, ub):
        S = np.array(S)
        KS = K[np.ix_(S, S)]
        bS = bv[S]
        AS = A[:, S]
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
        Gt = N.T @ KS @ N
        bt = N.T @ bS
        gt = N.T @ (KS @ c0)
        q0 = c0 @ KS @ c0
        gi = np.linalg.solve(Gt, gt)
        bi = np.linalg.solve(Gt, bt)
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
        S = [bestM, bestF]
        c_cur = np.zeros(n)
        c_cur[S] = 0.5
        dropped = np.zeros(n, bool)
        for _ in range(max_iter):
            sol = closed_form(S, ub)
            if sol is None:
                Kc = K @ c_cur
                order = [int(x) for x in np.argsort(Kc) if x not in S and not dropped[x]]
                if not order:
                    break
                S.extend(order[: max(1, len(S))])
                continue
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
                S = keep
                continue
            c = np.zeros(n)
            for idx, i in enumerate(S):
                c[i] = cS[idx]
            c_cur = c
            Kc = K @ c
            r = bv - np.where(male, muM, muF) - 2 * lam * Kc
            r[S] = -np.inf
            j = int(np.argmax(r))
            if r[j] <= tol:
                return c, sorted(S), _ + 1
            dropped[:] = False
            S.append(j)
        return c_cur, sorted(S), max_iter

    def ne_of(c):
        return 1.0 / np.sum(c**2)

    nm = male.sum()
    nf = (~male).sum()
    cg = np.zeros(n)
    cg[int(np.where(male)[0][np.argmax(bv[male])])] = 0.5
    cg[int(np.where(~male)[0][np.argmax(bv[~male])])] = 0.5
    ub_loose = cg @ K @ cg
    cu = np.where(male, 0.5 / nm, 0.5 / nf)
    ub_tight = cu @ K @ cu

    print(f"n={n}  K={K.shape[0]}x{K.shape[1]}")
    for t_ne in (25.0, 50.0, 100.0):
        lo, hi, chosen = ub_tight, ub_loose, 0.5 * (ub_tight + ub_loose)
        for _ in range(30):
            ub = 0.5 * (lo + hi)
            chosen = ub
            c, S, iters = support_first(ub)
            if ne_of(c) > t_ne:
                lo = ub
            else:
                hi = ub
            if abs(ne_of(c) - t_ne) / t_ne < 0.01:
                break
        ub = chosen
        ts = []
        for _ in range(3):
            t0 = time.perf_counter()
            c, S, iters = support_first(ub)
            ts.append(time.perf_counter() - t0)
        ne = ne_of(c)
        males = sum(1 for i in S if male[i])
        s_arr = np.array(S)
        kappa = np.linalg.cond(K[np.ix_(s_arr, s_arr)])  # conditioning the closed form faces
        print(
            f"  Ne={ne:5.1f} dF={50/ne:4.2f}% |S|={len(S):4d} ({males}m/{len(S)-males}f)"
            f" iters={iters:4d} algo={np.median(ts):.4f}s  cond(G_S)={kappa:.2e}"
        )


for n in (int(a) for a in sys.argv[1:]):
    run(n)
