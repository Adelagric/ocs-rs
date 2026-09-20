"""Route-A OCS as a CVXPY family: max bᵀc  s.t. ‖Lᵀc‖ ≤ r, Σc = 1, c ≥ 0.
Parameters: L (n×n Cholesky factor of G+εI), b (n), r (= √k). Generated with cvxgenrust."""
import sys, time, os
import cvxpy as cp, cvxgenrust as cgr

n = int(sys.argv[1])
L = cp.Parameter((n, n), name="L")
b = cp.Parameter(n, name="b")
r = cp.Parameter(nonneg=True, name="r")
c = cp.Variable(n, name="c")
problem = cp.Problem(cp.Maximize(b @ c), [cp.norm(L.T @ c, 2) <= r, cp.sum(c) == 1, c >= 0])
assert problem.is_dcp(dpp=True), "not DPP"

t0 = time.perf_counter()
proj = cgr.generate_code(problem, code_dir=f"ocs_n{n}_cgr", module_name=f"ocs_n{n}", wrapper=False)
t1 = time.perf_counter()
size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(f"ocs_n{n}_cgr/src") for f in fs)
print(f"[gen] n={n}: generate_code {t1-t0:.1f}s, generated src/ = {size/1e6:.1f} MB")
