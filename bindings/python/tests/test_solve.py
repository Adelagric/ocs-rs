"""Verify the Python bindings end to end.

The optimum is checked against an *independent* optimiser (scipy SLSQP on the dense
formulation), not just against the Rust solver itself — so a marshalling bug (row
order, sex mask, caps) cannot pass silently. The sexed and capped paths are checked
on their constraints, which are what those code paths exist to enforce.
"""

import numpy as np
import pytest

import ocs_rs

RIDGE = 1e-5


def make_data(n=200, m=400, seed=0):
    """Synthetic centred genotypes with a breeding value correlated to them."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.05, 0.5, m)
    genotypes = rng.binomial(2, p, size=(n, m)).astype(float)
    Z = genotypes - 2 * p
    s = float(2 * np.sum(p * (1 - p)))
    g = Z @ rng.standard_normal(m)
    b = (g - g.mean()) / g.std()
    return Z, s, b


def dense_G(Z, s, ridge=RIDGE):
    return Z @ Z.T / s + ridge * np.eye(Z.shape[0])


def binding_cap(G, frac=0.3):
    """A cap between the uniform-plan minimum and a typical single-candidate cost."""
    n = G.shape[0]
    one = np.ones(n)
    kmin = float(one @ G @ one) / n**2
    return kmin + frac * (float(np.diag(G).mean()) - kmin)


def test_simplex_matches_independent_optimiser():
    scipy_opt = pytest.importorskip("scipy.optimize")
    Z, s, b = make_data()
    G = dense_G(Z, s)
    n = len(b)
    k = binding_cap(G)

    res = ocs_rs.solve(Z, b, k=k, s=s, ridge=RIDGE)

    assert res.status == "Solved"
    assert abs(res.c.sum() - 1.0) < 1e-9
    assert res.c.min() > -1e-9
    assert res.quad <= k * (1 + 1e-9)

    cons = [
        {"type": "eq", "fun": lambda c: c.sum() - 1.0},
        {"type": "ineq", "fun": lambda c: k - c @ G @ c},
    ]
    ref = scipy_opt.minimize(
        lambda c: -(b @ c),
        np.ones(n) / n,
        bounds=[(0, None)] * n,
        constraints=cons,
        method="SLSQP",
        options={"ftol": 1e-12, "maxiter": 800},
    )
    assert abs(res.gain - (-ref.fun)) < 1e-6, (res.gain, -ref.fun)
    # the support is a handful of candidates, not the whole panel
    assert len(res.support) < n // 4


def test_sexed_splits_the_budget_by_sex():
    Z, s, b = make_data(seed=1)
    n = len(b)
    male = np.zeros(n, dtype=bool)
    male[::2] = True
    k = binding_cap(dense_G(Z, s))

    res = ocs_rs.solve(Z, b, k=k, s=s, ridge=RIDGE, male=male)

    assert res.status == "Solved"
    assert abs(res.c[male].sum() - 0.5) < 1e-9
    assert abs(res.c[~male].sum() - 0.5) < 1e-9
    assert res.quad <= k * (1 + 1e-9)


def test_caps_are_respected():
    Z, s, b = make_data(seed=2)
    n = len(b)
    k = binding_cap(dense_G(Z, s))
    caps = np.full(n, 0.05)

    res = ocs_rs.solve(Z, b, k=k, s=s, ridge=RIDGE, caps=caps)

    assert res.status == "Solved"
    assert res.c.max() <= 0.05 + 1e-9
    assert abs(res.c.sum() - 1.0) < 1e-9
    assert res.quad <= k * (1 + 1e-9)


def test_shape_errors_are_reported():
    Z, s, b = make_data(n=50, m=80, seed=3)
    with pytest.raises(ValueError):
        ocs_rs.solve(Z, b[:-1], k=0.1, s=s)          # b too short
    with pytest.raises(ValueError):
        ocs_rs.solve(Z[0], b, k=0.1, s=s)            # Z not 2-D
