"""Dense sexed prototype with an explicit primal/KKT acceptance gate."""
import csv
import hashlib
import os
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, null_space


class ConvergenceError(RuntimeError):
    pass


def load_panel(panel, legacy=False):
    if legacy:
        base = Path(os.environ.get("REPRO_LEGACY_DIR", "/tmp"))
        kp, bp = base / f"bench_K_{panel}.csv", base / f"bench_bc_{panel}.csv"
    else:
        base = Path(os.environ.get("REPRO_INPUT_DIR", "artifacts/repro_inputs"))
        kp, bp = base / f"{panel}_K.csv", base / f"{panel}_bc.csv"
    K = np.loadtxt(kp, delimiter=",")
    with bp.open() as stream:
        rows = list(csv.DictReader(stream))
    bv = np.array([float(r["bv"]) for r in rows])
    sex = np.array([r["sex"] for r in rows])
    if K.shape != (len(bv), len(bv)) or not np.isfinite(K).all() or not np.isfinite(bv).all():
        raise ValueError("nonfinite or misaligned dense panel")
    if not np.all(np.isin(sex, ["male", "female"])) or len(set(sex)) != 2:
        raise ValueError("both sex groups and aligned sex labels are required")
    if not np.allclose(K, K.T, atol=1e-10, rtol=1e-10):
        raise ValueError("kinship matrix must be symmetric")
    # Positive ridge makes automated genomic exports strictly positive definite.
    # Check once outside timed solves; symmetric indefinite K invalidates convex
    # KKT sufficiency and must not be reported as a certified global optimum.
    try:
        cholesky(K, lower=True, check_finite=False)
    except np.linalg.LinAlgError as error:
        raise ValueError("kinship input must be positive definite; use a positive ridge") from error
    hashes = {"K_md5": hashlib.md5(kp.read_bytes()).hexdigest(),
              "bc_md5": hashlib.md5(bp.read_bytes()).hexdigest()}
    return K, bv, sex == "male", hashes


def verify_solution(K, bv, male, c, ub, mu, lam, tol=1e-7):
    if not np.isfinite(ub) or ub <= 0:
        raise ConvergenceError("a finite positive kinship cap is required")
    ktol = 10 * tol * max(1.0, float(np.max(np.abs(bv))))
    if not np.isfinite(c).all() or not np.isfinite(lam) or not np.isfinite(mu).all():
        raise ConvergenceError("nonfinite solution or multipliers")
    q = float(c @ K @ c)
    if c.min() < -tol or abs(c[male].sum() - .5) > tol or abs(c[~male].sum() - .5) > tol:
        raise ConvergenceError("primal sex/nonnegativity constraints failed")
    if not np.isfinite(q):
        raise ConvergenceError("quadratic value must be finite")
    if q > ub + tol * max(1.0, abs(ub)) or lam < -tol:
        raise ConvergenceError("primal kinship or dual feasibility failed")
    residual = bv - np.where(male, mu[0], mu[1]) - 2 * lam * (K @ c)
    if not np.isfinite(residual).all() or not np.isfinite(lam * (q - ub)):
        raise ConvergenceError("KKT residual/complementarity must be finite")
    active = c > tol
    if np.any(np.abs(residual[active]) > ktol) or np.any(residual[~active] > ktol):
        raise ConvergenceError("stationarity/reduced-cost KKT check failed")
    if abs(lam * (q - ub)) > ktol * max(1.0, abs(ub)):
        raise ConvergenceError("quadratic complementarity failed")


def support_first(K, bv, male, ub, tol=1e-7, max_iter=8000):
    if not np.isfinite(ub) or ub <= 0:
        raise ConvergenceError("a finite positive kinship cap is required")
    n = len(bv)
    A = np.vstack([male.astype(float), (~male).astype(float)])
    rhs = np.array([.5, .5])
    bestM = int(np.where(male)[0][np.argmax(bv[male])])
    bestF = int(np.where(~male)[0][np.argmax(bv[~male])])
    S = [bestM, bestF]
    c_cur = np.zeros(n)
    c_cur[S] = .5
    if float(c_cur @ K @ c_cur) <= ub + tol:
        mu = bv[[bestM, bestF]]
        verify_solution(K, bv, male, c_cur, ub, mu, 0., tol)
        return c_cur, sorted(S), 0
    dropped = np.zeros(n, bool)

    def closed_form():
        idx = np.array(S)
        KS, bS, AS = K[np.ix_(idx, idx)], bv[idx], A[:, idx]
        c0 = AS.T @ np.linalg.solve(AS @ AS.T, rhs)
        N = null_space(AS)
        if N.shape[1] == 0:
            return None
        Gt, bt, gt = N.T @ KS @ N, N.T @ bS, N.T @ (KS @ c0)
        gi, bi = np.linalg.solve(Gt, gt), np.linalg.solve(Gt, bt)
        rho2 = float(gt @ gi - (c0 @ KS @ c0 - ub))
        if rho2 <= 0:
            return None
        bnorm = float(bt @ bi)
        if bnorm <= 1e-24:
            cS, lam = c0 - N @ gi, 0.
        else:
            cS = c0 + N @ (-gi + np.sqrt(rho2 / bnorm) * bi)
            lam = .5 * np.sqrt(bnorm / rho2)
        mu = np.linalg.lstsq(AS.T, bS - 2 * lam * (KS @ cS), rcond=None)[0]
        return cS, mu, lam

    for it in range(max_iter):
        sol = closed_form()
        if sol is None:
            order = [int(x) for x in np.argsort(K @ c_cur) if x not in S and not dropped[x]]
            if not order:
                raise ConvergenceError("no feasible support found; cap may be infeasible")
            S.extend(order[:max(1, len(S))])
            continue
        cS, mu, lam = sol
        if cS.min() < -tol:
            keep = [S[i] for i in range(len(S)) if cS[i] > tol]
            dropped[np.array(S)[cS <= tol]] = True
            if not any(male[i] for i in keep):
                keep.append(bestM)
            if not any(not male[i] for i in keep):
                keep.append(bestF)
            S = keep
            continue
        c = np.zeros(n)
        c[S] = cS
        c_cur = c
        residual = bv - np.where(male, mu[0], mu[1]) - 2 * lam * (K @ c)
        residual[S] = -np.inf
        j = int(np.argmax(residual))
        if residual[j] <= tol:
            verify_solution(K, bv, male, c, ub, mu, lam, tol)
            return c, np.flatnonzero(c > tol).tolist(), it + 1
        dropped[:] = False
        S.append(j)
    raise ConvergenceError(f"iteration limit ({max_iter}) reached before KKT convergence")
