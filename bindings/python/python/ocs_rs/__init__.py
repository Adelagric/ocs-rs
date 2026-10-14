"""Exact, matrix-free optimum contribution selection (OCS) at genomic scale.

OCS maximises genetic gain ``bᵀc`` subject to a quadratic relationship bound
``cᵀGc ≤ k``, with non-negative contributions summing to one (or,
in the sexed form, to one half per sex). This package wraps the Rust `support-first`
solver: it finds the active support by column generation, solves each fixed
support in closed form, and never forms the dense n×n relationship matrix — the
relationship products come straight from the genotype matrix ``Z``. Small support
is empirical; individual caps can require large supports. With the additive
relationship convention ``G = Z Zᵀ / s + ridge·I``, group coancestry is
``cᵀ(G/2)c = quad/2`` and its bound is ``k/2``. This is not a realised offspring
inbreeding rate.

The solver targets the optimum for continuous contributions under a single
coancestry constraint, subject to its numerical termination tolerances. A complete
runnable example is provided in the package README and examples/quickstart.py in
the source checkout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _ocs_rs as _rs

__all__ = ["solve", "read_plink", "OcsResult", "Panel", "vanraden_scale"]
__version__ = "0.4.2"


@dataclass
class OcsResult:
    """Outcome of an OCS solve.

    ``c`` and ``gain`` are a usable optimum only when ``status == "Solved"``.
    ``support`` contains 0-based NumPy indices.
    """

    c: np.ndarray
    support: np.ndarray
    gain: float
    quad: float
    iterations: int
    products: int
    status: str

    @property
    def solved(self) -> bool:
        return self.status == "Solved"

    def __repr__(self) -> str:  # concise, the vectors are usually long
        return (
            f"OcsResult(status={self.status!r}, |support|={len(self.support)}, "
            f"gain={self.gain:.6g}, quad={self.quad:.6g}, products={self.products})"
        )


@dataclass
class Panel:
    """Dense marker panel; ``ids`` contains IIDs in ``.fam`` row order.

    FID and sex metadata are not returned. Align all external vectors to these
    rows before handing ``Z`` to :func:`solve`.
    """

    Z: np.ndarray
    p: np.ndarray
    s: float
    ids: list

    @property
    def n(self) -> int:
        return int(self.Z.shape[0])

    @property
    def m(self) -> int:
        return int(self.Z.shape[1])

    def __repr__(self) -> str:  # concise, Z is usually large
        return f"Panel(n={self.n}, m={self.m}, s={self.s:.6g})"


def read_plink(prefix) -> Panel:
    """Read a PLINK 1 trio ``<prefix>.bed`` / ``.bim`` / ``.fam``.

    Genotypes are centred by twice the allele frequency estimated from the
    non-missing calls, and missing calls are imputed to the marker mean — so
    ``panel.Z`` and ``panel.s`` can be used by :func:`solve`:

    >>> panel = ocs_rs.read_plink("cattle")
    >>> panel.ids  # Align external breeding values, sex masks and caps to these IIDs.

    Only SNP-major ``.bed`` files are read (PLINK ≥ 1.9 writes those by default).
    ``ids`` contains the IID column only; FID and ``.fam`` sex are not returned.
    Preserve them separately and check IID uniqueness before matching. File
    dimension checks do not verify identity or order agreement among the trio.
    This reader returns a dense matrix; packed solving is Rust-only.
    """
    out = _rs.read_plink(str(prefix))
    # `out["z"]` is raw row-major float64 bytes; frombuffer views them without a copy,
    # reshape lays out the (n, m) matrix. A copy makes it writable for the caller.
    z = np.frombuffer(out["z"], dtype=np.float64).reshape(int(out["n"]), int(out["m"])).copy()
    return Panel(
        Z=z,
        p=np.asarray(out["p"], dtype=np.float64),
        s=float(out["s"]),
        ids=list(out["ids"]),
    )


def vanraden_scale(freqs) -> float:
    """VanRaden scaling ``s = 2 Σ p(1-p)`` from allele frequencies ``p``."""
    p = np.asarray(freqs, dtype=np.float64)
    return float(2.0 * np.sum(p * (1.0 - p)))


def solve(
    Z,
    b,
    k: float,
    *,
    s: float,
    ridge: float = 1e-5,
    male=None,
    caps=None,
    max_iter: int = 10_000,
    tol: float = 1e-9,
) -> OcsResult:
    """Solve OCS on the centred genotype matrix ``Z``.

    Parameters
    ----------
    Z : array (n, m)
        Genotypes centred by twice the allele frequencies. The relationship matrix
        is ``G = Z Zᵀ / s + ridge·I``; it is never formed.
    b : array (n,)
        Breeding values (the selection criterion), in exactly the row order of Z.
    k : float
        Quadratic relationship bound: ``cᵀGc ≤ k``. Group coancestry uses ``G/2``
        and therefore has bound ``k/2``.
    s : float
        VanRaden scale ``2 Σ p(1-p)`` — see :func:`vanraden_scale`.
    ridge : float, optional
        Positive ridge making ``G`` positive definite (default 1e-5). Its term
        ``ridge * sum(c**2)`` is included in the quadratic relationship constraint.
    male : array (n,) of bool, optional
        Sex indicator. When given, the sexed form is solved
        (``Σ_males = Σ_females = ½``); otherwise the simplex form (``Σc = 1``).
        Must be in the row order of Z; :func:`read_plink` does not return sex.
    caps : array (n,), optional
        Per-candidate upper bounds ``c ≤ caps``, in the row order of Z.
        Omit for no caps.
    max_iter, tol : optional
        Active-set iteration cap and reduced-cost tolerance.

    Returns
    -------
    OcsResult
        Check ``status == "Solved"`` before using the contributions as an optimum.
        ``support`` contains 0-based indices; ``MaxIter`` has no optimality guarantee.

    Notes
    -----
    This function accepts Z, not G. Passing G as Z would use
    ``G @ G.T / s + ridge * I``. Arbitrary or modified relationship matrices
    require an adapter/factorization that reproduces the intended operator.
    Python makes dense copies of Z (contiguous float64, bytes, native layout).
    Packed 2-bit solving is available only through Rust's ``read_packed`` and
    ``solve_sexed_packed`` APIs.
    """
    Z = np.ascontiguousarray(Z, dtype=np.float64)
    if Z.ndim != 2:
        raise ValueError(f"Z must be 2-D (n, m), got shape {Z.shape}")
    n, m = Z.shape
    b = np.ascontiguousarray(b, dtype=np.float64).ravel()
    if b.shape[0] != n:
        raise ValueError(f"b has length {b.shape[0]}, expected n = {n}")
    # Raw row-major float64 bytes, not a Python-level sequence: the native side
    # reads them in one pass, avoiding an element-by-element copy through the
    # interpreter that dominates the call once the solve is milliseconds. Z is
    # already C-contiguous float64 above, so tobytes() is a single memcpy.
    zb = Z.tobytes()

    if male is None and caps is None:
        out = _rs.solve(zb, n, m, s, ridge, b, k, max_iter, tol)
    elif male is None:
        c_ = np.ascontiguousarray(caps, dtype=np.float64).ravel()
        out = _rs.solve_capped(zb, n, m, s, ridge, b, c_, k, max_iter, tol)
    elif caps is None:
        mk = np.ascontiguousarray(male, dtype=bool).ravel().tolist()
        out = _rs.solve_sexed(zb, n, m, s, ridge, b, mk, k, max_iter, tol)
    else:
        mk = np.ascontiguousarray(male, dtype=bool).ravel().tolist()
        c_ = np.ascontiguousarray(caps, dtype=np.float64).ravel()
        out = _rs.solve_sexed_capped(zb, n, m, s, ridge, b, mk, c_, k, max_iter, tol)

    return OcsResult(
        c=np.asarray(out["c"], dtype=np.float64),
        support=np.asarray(out["support"], dtype=np.int64),
        gain=float(out["gain"]),
        quad=float(out["quad"]),
        iterations=int(out["iterations"]),
        products=int(out["products"]),
        status=str(out["status"]),
    )
