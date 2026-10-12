"""Run with Python after installing ocs-rs and NumPy."""

import numpy as np
import ocs_rs

rng = np.random.default_rng(11)
n, m = 12, 20
dosages = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m))
p = dosages.mean(axis=0) / 2
Z = dosages - 2 * p
s = ocs_rs.vanraden_scale(p)
b = rng.standard_normal(n)
male = np.arange(n) % 2 == 0
caps = np.full(n, 0.15)
ridge = 1e-5  # G = Z @ Z.T / s + ridge * I

# A feasible plan contributes half from each sex, uniformly within that sex.
c0 = np.where(male, 0.5 / male.sum(), 0.5 / (~male).sum())
assert np.all(c0 <= caps)
q0 = float(np.sum((Z.T @ c0) ** 2) / s + ridge * (c0 @ c0))
k = q0 + 0.05  # This cap admits c0; it is not a universal breeding threshold.

res = ocs_rs.solve(Z, b, k=k, s=s, ridge=ridge, male=male, caps=caps)
if res.status != "Solved":
    raise RuntimeError(f"No solved optimum: {res.status}")
print(res)
print(res.support)  # 0-based NumPy indices
print(res.c[res.support])
