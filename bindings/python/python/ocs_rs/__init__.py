"""Exact, matrix-free optimum contribution selection (OCS) at genomic scale.

OCS maximises genetic gain ``bᵀc`` subject to a cap on the mean coancestry of the
next generation, ``cᵀGc ≤ k``, with non-negative contributions summing to one (or,
in the sexed form, to one half per sex). This package wraps the Rust `support-first`
solver: it finds the tiny active support by column generation, solves each fixed
support in closed form, and never forms the dense n×n relationship matrix — the
kinship products come straight from the genotype matrix ``Z``.

The result is the exact optimum (KKT-certified, not a heuristic) for continuous
contributions under a single coancestry constraint.

Example
-------
>>> import numpy as np, ocs_rs
>>> res = ocs_rs.solve(Z, b, k=0.03, s=s, male=male)   # Z: (n, m) centred genotypes
>>> res.support, res.gain
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import _ocs_rs as _rs

__all__ = ["solve", "read_plink", "OcsResult", "Panel", "vanraden_scale"]
__version__ = "0.2.0"


@dataclass
class OcsResult:
    """Outcome of an OCS solve.

    ``c`` and ``gain`` are a usable optimum only when ``status == "Solved"``.
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
    """A marker panel read from PLINK files, ready to hand to :func:`solve`."""

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
    ``panel.Z`` and ``panel.s`` go straight into :func:`solve`:

    >>> panel = ocs_rs.read_plink("cattle")
    >>> res = ocs_rs.solve(panel.Z, b, k=0.03, s=panel.s)

    Only SNP-major ``.bed`` files are read (PLINK ≥ 1.9 writes those by default).
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
        Breeding values (the selection criterion).
    k : float
        Coancestry cap: the constraint is ``cᵀGc ≤ k``.
    s : float
        VanRaden scale ``2 Σ p(1-p)`` — see :func:`vanraden_scale`.
    ridge : float, optional
        Small ridge making ``G`` positive definite (default 1e-5).
    male : array (n,) of bool, optional
        Sex indicator. When given, the sexed form is solved
        (``Σ_males = Σ_females = ½``); otherwise the simplex form (``Σc = 1``).
    caps : array (n,), optional
        Per-candidate upper bounds ``c ≤ caps``. Omit for no caps.
    max_iter, tol : optional
        Active-set iteration cap and reduced-cost tolerance.

    Returns
    -------
    OcsResult
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
