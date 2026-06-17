//! Casting OCS into Clarabel's standard conic form.
//!
//! Clarabel solves `min ½xᵀPx + qᵀx  s.t.  Ax + s = b, s ∈ K`. OCS maps to it
//! with `P = 0`, `q = -b` (maximising genetic gain), and the constraint stack
//!
//! | block        | rows  | A          | b        | cone                 |
//! |--------------|-------|------------|----------|----------------------|
//! | `Σcᵢ = 1`    | 1     | `1ᵀ`       | `1`      | `ZeroConeT(1)`       |
//! | `c ≥ 0`      | n     | `-I`       | `0`      | `NonnegativeConeT(n)`|
//! | `‖Fᵀc‖ ≤ r`  | d+1   | `[0ᵀ;-Fᵀ]` | `[r,0…]` | `SecondOrderConeT`   |
//! | `c ≤ u`(opt) | n     | `I`        | `u`      | `NonnegativeConeT(n)`|
//!
//! The SOC slack is `s = b − Ax = (r, Fᵀc)`, so `‖Fᵀc‖ ≤ r`. With the Cholesky
//! factor `F = L` (`d = n`, `r = √k`) the kinship form is `‖Lᵀc‖² = cᵀ(G+εI)c`;
//! with the raw factor `F = Z` (`d = m`, `r = √(k·s)`) it is `‖Zᵀc‖² = s·cᵀGc`.
//! Both enforce `cᵀGc ≤ k` (Route A conservatively by `ε‖c‖²`).
//!
//! `A` is built directly in CSC: per candidate column the entries are emitted in
//! strictly increasing row order, so no sort or triplet buffer is needed even at
//! Route-B scale where the conic block alone has `n·m` nonzeros.

use clarabel::algebra::CscMatrix;
use clarabel::solver::{NonnegativeConeT, SecondOrderConeT, SupportedConeT, ZeroConeT};
use faer::Mat;

/// Which factorisation supplies the second-order cone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Route {
    /// Cholesky factor `L` of `G+εI`; SOC dimension `n+1`.
    Cholesky,
    /// Raw centred genotypes `Z`; SOC dimension `m+1`, no Cholesky.
    Raw,
}

impl Route {
    /// Short tag for CSV/report output.
    pub fn tag(self) -> &'static str {
        match self {
            Route::Cholesky => "A_chol",
            Route::Raw => "B_raw",
        }
    }
}

/// The factor matrix backing the cone, borrowed from the caller.
pub enum Factor<'a> {
    /// Lower-triangular Cholesky factor `L`, shape `n×n`.
    Cholesky(&'a Mat<f64>),
    /// Centred genotype matrix `Z`, shape `n×m`.
    Raw(&'a Mat<f64>),
}

impl Factor<'_> {
    /// Cone tail dimension `d` (`n` for Cholesky, `m` for raw).
    pub fn d(&self) -> usize {
        match self {
            Factor::Cholesky(l) => l.nrows(),
            Factor::Raw(z) => z.ncols(),
        }
    }

    /// The route this factor implements.
    pub fn route(&self) -> Route {
        match self {
            Factor::Cholesky(_) => Route::Cholesky,
            Factor::Raw(_) => Route::Raw,
        }
    }
}

/// A fully-assembled conic problem ready for [`crate::solve`].
pub struct ConicProblem {
    /// Zero objective Hessian, `n×n`.
    pub p: CscMatrix<f64>,
    /// Linear objective `-b`.
    pub q: Vec<f64>,
    /// Stacked constraint matrix in CSC.
    pub a: CscMatrix<f64>,
    /// Right-hand side.
    pub b: Vec<f64>,
    /// Cone composition, in row order.
    pub cones: Vec<SupportedConeT<f64>>,
    /// Number of decision variables (candidates).
    pub n: usize,
    /// SOC tail dimension actually used.
    pub d: usize,
    /// Route used.
    pub route: Route,
    /// Number of structural nonzeros in `A` (a proxy for conic-block density).
    pub a_nnz: usize,
}

/// Assemble the OCS conic program.
///
/// `b_gebv` are the genetic-gain coefficients (length `n`). `k` is the kinship
/// bound, `s` the VanRaden scale (only used by [`Route::Raw`]). `cap` is the
/// optional per-candidate upper bound `u` (length `n`); `None` leaves it off.
pub fn build(
    factor: Factor<'_>,
    b_gebv: &[f64],
    k: f64,
    s: f64,
    cap: Option<&[f64]>,
) -> ConicProblem {
    let route = factor.route();
    let d = factor.d();
    let n = b_gebv.len();
    let r = match route {
        Route::Cholesky => k.sqrt(),
        Route::Raw => (k * s).sqrt(),
    };

    let has_cap = cap.is_some();
    let n_rows = 1 + n + (d + 1) + if has_cap { n } else { 0 };

    // Row offsets within the stacked constraint matrix.
    let soc_r_row = 1 + n; // the radius component r
    let soc_tail0 = soc_r_row + 1; // first row of Fᵀc
    let cap_row0 = soc_tail0 + d; // first row of the c ≤ u block

    // Upper bound on nonzeros: eq(1) + nonneg(1) [+cap(1)] per column, plus the
    // conic tail (≤ j+1 for Cholesky, m for raw).
    let per_col_fixed = 1 + 1 + usize::from(has_cap);
    let soc_nnz = match factor {
        Factor::Cholesky(_) => n * (n + 1) / 2,
        Factor::Raw(_) => n * d,
    };
    let nnz_cap = n * per_col_fixed + soc_nnz;

    let mut colptr = Vec::with_capacity(n + 1);
    let mut rowval = Vec::with_capacity(nnz_cap);
    let mut nzval = Vec::with_capacity(nnz_cap);
    colptr.push(0);

    for j in 0..n {
        // Σc = 1
        rowval.push(0);
        nzval.push(1.0);
        // c ≥ 0  (row 1+j of -I)
        rowval.push(1 + j);
        nzval.push(-1.0);
        // SOC tail: column j of -Fᵀ
        match factor {
            Factor::Cholesky(l) => {
                // -Lᵀ column j = -(L row j), nonzero for i = 0..=j.
                for i in 0..=j {
                    rowval.push(soc_tail0 + i);
                    nzval.push(-l[(j, i)]);
                }
            }
            Factor::Raw(z) => {
                // -Zᵀ column j = -(Z row j), dense over i = 0..m.
                for i in 0..d {
                    rowval.push(soc_tail0 + i);
                    nzval.push(-z[(j, i)]);
                }
            }
        }
        // c ≤ u  (row cap_row0+j of I)
        if has_cap {
            rowval.push(cap_row0 + j);
            nzval.push(1.0);
        }
        colptr.push(rowval.len());
    }

    let a_nnz = rowval.len();
    let a = CscMatrix::new(n_rows, n, colptr, rowval, nzval);
    let p = CscMatrix::<f64>::zeros((n, n));
    let q: Vec<f64> = b_gebv.iter().map(|&v| -v).collect();

    // b: [1 | 0…0 | r 0…0 | u]
    let mut b = Vec::with_capacity(n_rows);
    b.push(1.0);
    b.extend(std::iter::repeat_n(0.0, n));
    b.push(r);
    b.extend(std::iter::repeat_n(0.0, d));
    if let Some(u) = cap {
        b.extend_from_slice(u);
    }
    debug_assert_eq!(b.len(), n_rows);

    let mut cones: Vec<SupportedConeT<f64>> = Vec::with_capacity(4);
    cones.push(ZeroConeT(1));
    cones.push(NonnegativeConeT(n));
    cones.push(SecondOrderConeT(d + 1));
    if has_cap {
        cones.push(NonnegativeConeT(n));
    }

    ConicProblem {
        p,
        q,
        a,
        b,
        cones,
        n,
        d,
        route,
        a_nnz,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clarabel::algebra::CscMatrix;

    // Densify a CscMatrix into row-major Vec<Vec<f64>> for comparison.
    // `col` indexes the column-pointer pair `colptr[col]..colptr[col+1]`, the
    // canonical CSC traversal, so an `enumerate` rewrite does not apply.
    #[allow(clippy::needless_range_loop)]
    fn densify(a: &CscMatrix<f64>) -> Vec<Vec<f64>> {
        let mut dense = vec![vec![0.0; a.n]; a.m];
        for col in 0..a.n {
            for idx in a.colptr[col]..a.colptr[col + 1] {
                dense[a.rowval[idx]][col] = a.nzval[idx];
            }
        }
        dense
    }

    #[test]
    fn direct_csc_matches_triplets_cholesky() {
        // Build the same A via the direct column path and via new_from_triplets,
        // and assert they densify identically. This guards the assembler.
        let n = 5;
        let l = Mat::<f64>::from_fn(n, n, |i, j| if j <= i { (i + j + 1) as f64 } else { 0.0 });
        let b = vec![0.3, -0.1, 0.7, 0.2, 0.5];
        let prob = build(Factor::Cholesky(&l), &b, 0.4, 1.0, None);

        let d = n;
        let mut ri = Vec::new();
        let mut ci = Vec::new();
        let mut vi = Vec::new();
        for j in 0..n {
            ri.push(0);
            ci.push(j);
            vi.push(1.0);
            ri.push(1 + j);
            ci.push(j);
            vi.push(-1.0);
            for i in 0..=j {
                ri.push(2 + n + i);
                ci.push(j);
                vi.push(-l[(j, i)]);
            }
        }
        let a_ref = CscMatrix::new_from_triplets(1 + n + (d + 1), n, ri, ci, vi);
        assert_eq!(densify(&prob.a), densify(&a_ref));
    }

    #[test]
    fn raw_route_dimensions_and_density() {
        let n = 4;
        let m = 7;
        let z = Mat::<f64>::from_fn(n, m, |i, j| (i as f64) - (j as f64) * 0.5);
        let b = vec![1.0, 2.0, 3.0, 4.0];
        let prob = build(Factor::Raw(&z), &b, 0.5, 2.0, None);
        assert_eq!(prob.d, m);
        assert_eq!(prob.route, Route::Raw);
        // SOC tail is fully dense: n*m entries.
        assert_eq!(prob.a_nnz, n * 2 + n * m);
        // r = sqrt(k*s)
        assert!((prob.b[1 + n] - (0.5_f64 * 2.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn soc_slack_reproduces_quadratic_form() {
        // For an explicit c, the SOC tail rows of (b - A c) must equal Fᵀc, so
        // that ‖tail‖² is the kinship quadratic form. Check with Cholesky.
        let n = 3;
        let l = Mat::<f64>::from_fn(n, n, |i, j| if j <= i { 1.0 + (i * j) as f64 } else { 0.0 });
        let b = vec![0.5, 0.5, 0.5];
        let k = 0.9;
        let prob = build(Factor::Cholesky(&l), &b, k, 1.0, None);
        let c = [0.2, 0.3, 0.5];

        // s = b - A c, computed from the CSC (canonical column traversal).
        let mut ax = vec![0.0; prob.a.m];
        #[allow(clippy::needless_range_loop)]
        for col in 0..prob.a.n {
            for idx in prob.a.colptr[col]..prob.a.colptr[col + 1] {
                ax[prob.a.rowval[idx]] += prob.a.nzval[idx] * c[col];
            }
        }
        let slack: Vec<f64> = prob.b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();

        // SOC block starts at row 1+n: (r, (Lᵀc)_0, …).
        let r = slack[1 + n];
        assert!((r - k.sqrt()).abs() < 1e-12);
        let tail = &slack[1 + n + 1..1 + n + 1 + n];
        // ‖tail‖² == cᵀ L Lᵀ c.
        let tail_norm_sq: f64 = tail.iter().map(|x| x * x).sum();
        let mut llt_quad = 0.0;
        for i in 0..n {
            for jj in 0..n {
                let mut g_ij = 0.0;
                for t in 0..n {
                    g_ij += l[(i, t)] * l[(jj, t)];
                }
                llt_quad += c[i] * g_ij * c[jj];
            }
        }
        assert!((tail_norm_sq - llt_quad).abs() < 1e-10);
    }

    #[test]
    #[allow(clippy::needless_range_loop)] // dense[row][col] diagonal check
    fn cap_block_is_identity_with_u_rhs() {
        // The optional c ≤ u block must be +I in A (last n rows) and carry u in
        // b. This path is never exercised by the experiments, so it is pinned
        // here directly against the densified matrix.
        let n = 3;
        let l = Mat::<f64>::identity(n, n);
        let b = vec![1.0, 2.0, 3.0];
        let u = vec![0.4, 0.5, 0.6];
        let prob = build(Factor::Cholesky(&l), &b, 0.5, 1.0, Some(&u));

        // Rows: 1 (eq) + n (nonneg) + (n+1) (soc) + n (cap).
        assert_eq!(prob.a.m, 1 + n + (n + 1) + n);
        assert_eq!(prob.cones.len(), 4);

        let dense = densify(&prob.a);
        let cap0 = 1 + n + (n + 1); // first cap row
        for i in 0..n {
            for j in 0..n {
                let expect = if i == j { 1.0 } else { 0.0 };
                assert_eq!(dense[cap0 + i][j], expect, "cap block not +I at ({i},{j})");
            }
            // b carries u on the cap rows.
            assert_eq!(prob.b[cap0 + i], u[i]);
        }
    }
}
