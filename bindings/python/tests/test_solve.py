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
    """A cap between a feasible uniform plan and a typical single-candidate cost."""
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
    assert ref.success, ref.message
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


def test_sexed_individual_caps_match_analytic_optimum():
    # Zero centred marker columns and ridge=1 give G=I, independently of layout.
    # The two high-value candidates saturate their distinct caps. Remaining male
    # and female budgets are .3 and .25. Write their paired contributions as
    # (.15+d_m, .15-d_m) and (.125+d_f, .125-d_f):
    # q = .17875 + 2*(d_m**2+d_f**2), gain varies as 2*d_m+d_f.
    # Cauchy-Schwarz therefore gives d_m=2*sqrt(.001), d_f=sqrt(.001) at k=.18875.
    # Their positive upper-bound multipliers confirm that saturating the first
    # two caps is optimal; this is an analytic reference, not another solver call.
    Z = np.zeros((6, 2))
    b = np.array([10.0, 9.0, 2.0, 1.0, 0.0, 0.0])
    male = np.array([True, False, True, False, True, False])
    caps = np.array([0.2, 0.25, 0.4, 0.4, 0.4, 0.4])
    k = 0.18875
    delta = np.sqrt(0.001)
    expected = np.array([0.2, 0.25, 0.15 + 2 * delta, 0.125 + delta,
                         0.15 - 2 * delta, 0.125 - delta])

    res = ocs_rs.solve(Z, b, k=k, s=1.0, ridge=1.0, male=male, caps=caps)

    assert res.status == "Solved"
    assert np.all(res.c >= -1e-10)
    assert np.all(res.c <= caps + 1e-10)
    assert abs(res.c[male].sum() - 0.5) < 1e-10
    assert abs(res.c[~male].sum() - 0.5) < 1e-10
    assert res.c @ res.c <= k + 1e-10
    np.testing.assert_allclose(res.c, expected, atol=1e-8, rtol=0)
    assert abs(res.gain - b @ expected) < 1e-8
    assert abs(res.quad - res.c @ res.c) < 1e-10
    np.testing.assert_array_equal(np.sort(res.support), np.flatnonzero(res.c > 1e-12))


def test_shape_errors_are_reported():
    Z, s, b = make_data(n=50, m=80, seed=3)
    with pytest.raises(ValueError):
        ocs_rs.solve(Z, b[:-1], k=0.1, s=s)          # b too short
    with pytest.raises(ValueError):
        ocs_rs.solve(Z[0], b, k=0.1, s=s)            # Z not 2-D
