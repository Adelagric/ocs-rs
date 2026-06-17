//! VanRaden genomic relationship matrix and its Cholesky factor (Route A).
//!
//! `G = Z Zᵀ / s + ε I`. The matrix kept in memory is the *ridged* one, because
//! that is what Route A factors and what Clarabel's cone actually enforces. The
//! true (unridged) quadratic form is recovered analytically when needed
//! ([`quad_form`]), avoiding a second `n×n` allocation at scale.

use crate::error::OcsError;
use faer::{Mat, Side};

/// A ridged genomic relationship matrix.
pub struct Grm {
    /// `G + εI`, symmetric positive definite, shape `n×n`.
    pub g: Mat<f64>,
    /// The ridge `ε` baked into `g`.
    pub ridge: f64,
    /// Dimension.
    pub n: usize,
}

impl Grm {
    /// Build `G = Z Zᵀ / s + ridge·I` from centred genotypes.
    ///
    /// Computed in place: the `Z Zᵀ` GEMM result is scaled by `1/s` and the
    /// ridge is added to its diagonal, so no extra `n×n` buffer is allocated
    /// beyond the GEMM output itself.
    pub fn build(z: &Mat<f64>, s: f64, ridge: f64) -> Grm {
        let n = z.nrows();
        let mut g = z * z.transpose(); // n×n, symmetric PSD
        let inv_s = 1.0 / s;
        for j in 0..n {
            for i in 0..n {
                g[(i, j)] *= inv_s;
            }
            g[(j, j)] += ridge;
        }
        Grm { g, ridge, n }
    }

    /// Lower Cholesky factor `L` of the ridged `G` (so `L Lᵀ = G + εI`).
    ///
    /// Returns [`OcsError::NonPositivePivot`] if the ridge is too small for the
    /// conditioning of this instance.
    pub fn cholesky_lower(&self) -> Result<Mat<f64>, OcsError> {
        match self.g.llt(Side::Lower) {
            Ok(llt) => Ok(llt.L().to_owned()),
            Err(_) => {
                // faer reports `NonPositivePivot { index }`; surface the index.
                let index = self.first_nonpositive_pivot_hint();
                Err(OcsError::NonPositivePivot {
                    index,
                    ridge: self.ridge,
                })
            }
        }
    }

    // faer's LltError carries the failing index, but its fields are not all
    // public across versions; this re-derives a best-effort hint cheaply for
    // reporting only (not used in any numerical path).
    fn first_nonpositive_pivot_hint(&self) -> usize {
        for i in 0..self.n {
            if self.g[(i, i)] <= 0.0 {
                return i;
            }
        }
        0
    }
}

/// Build `G` and factor it, escalating the ridge by ×10 (up to `max_tries`)
/// whenever Cholesky reports a non-positive pivot. Returns the factored `Grm`
/// (carrying the ridge that finally worked), its lower factor `L`, and the
/// number of escalations performed (0 = first ridge sufficed).
pub fn build_and_factor(
    z: &Mat<f64>,
    s: f64,
    initial_ridge: f64,
    max_tries: u32,
) -> Result<(Grm, Mat<f64>, u32), OcsError> {
    let mut ridge = initial_ridge;
    for attempt in 0..max_tries {
        let grm = Grm::build(z, s, ridge);
        match grm.cholesky_lower() {
            Ok(l) => return Ok((grm, l, attempt)),
            Err(OcsError::NonPositivePivot { .. }) => {
                ridge *= 10.0;
            }
            Err(e) => return Err(e),
        }
    }
    Err(OcsError::RidgeExhausted {
        ridge: ridge / 10.0,
    })
}

/// True quadratic form `cᵀ G_true c`, where `G_true = G_ridged − εI`.
///
/// This is what the kinship constraint bounds; the ridge is only a Cholesky
/// stabiliser. `cᵀ(G−εI)c = cᵀG c − ε‖c‖²`.
pub fn quad_form(grm: &Grm, c: &[f64]) -> f64 {
    let n = grm.n;
    debug_assert_eq!(c.len(), n);
    let mut acc = 0.0;
    let mut norm_sq = 0.0;
    for (i, &ci) in c.iter().enumerate() {
        norm_sq += ci * ci;
        let mut row = 0.0;
        for (j, &cj) in c.iter().enumerate() {
            row += grm.g[(i, j)] * cj;
        }
        acc += ci * row;
    }
    acc - grm.ridge * norm_sq
}

/// True quadratic form via the raw factor: `cᵀG c = ‖Zᵀc‖² / s` (Route B), with
/// no `G` ever formed. `y_i = (Zᵀc)_i = Σ_k Z[k,i] c[k]`.
pub fn quad_form_z(z: &Mat<f64>, s: f64, c: &[f64]) -> f64 {
    let m = z.ncols();
    debug_assert_eq!(c.len(), z.nrows());
    let mut acc = 0.0;
    for i in 0..m {
        let mut yi = 0.0;
        for (k, &ck) in c.iter().enumerate() {
            yi += z[(k, i)] * ck;
        }
        acc += yi * yi;
    }
    acc / s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::datagen;

    #[test]
    fn grm_is_symmetric() {
        let d = datagen::generate(40, 200, 11);
        let grm = Grm::build(&d.z, d.s, 1e-5);
        for i in 0..40 {
            for j in 0..40 {
                assert!((grm.g[(i, j)] - grm.g[(j, i)]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn grm_diagonal_near_one() {
        // VanRaden self-relationships scatter around 1.
        let d = datagen::generate(50, 5000, 3);
        let grm = Grm::build(&d.z, d.s, 1e-5);
        let mean_diag: f64 = (0..50).map(|i| grm.g[(i, i)]).sum::<f64>() / 50.0;
        assert!(
            (mean_diag - 1.0).abs() < 0.15,
            "mean diagonal {mean_diag} far from 1"
        );
    }

    #[test]
    fn cholesky_reconstructs_g() {
        let d = datagen::generate(30, 300, 21);
        let grm = Grm::build(&d.z, d.s, 1e-5);
        let l = grm.cholesky_lower().unwrap();
        // L Lᵀ must equal the ridged G on the lower triangle.
        for i in 0..30 {
            for j in 0..=i {
                let mut acc = 0.0;
                for t in 0..=j {
                    acc += l[(i, t)] * l[(j, t)];
                }
                assert!((acc - grm.g[(i, j)]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn quad_form_routes_agree() {
        // cᵀG c computed via the ridged G (minus ridge) and via Z must match.
        let d = datagen::generate(60, 400, 8);
        let grm = Grm::build(&d.z, d.s, 1e-5);
        let c: Vec<f64> = (0..60).map(|i| ((i % 7) as f64 + 1.0) / 100.0).collect();
        let via_g = quad_form(&grm, &c);
        let via_z = quad_form_z(&d.z, d.s, &c);
        assert!(
            (via_g - via_z).abs() < 1e-7,
            "G route {via_g} vs Z route {via_z}"
        );
        assert!(via_g > 0.0);
    }

    #[test]
    fn rank_deficient_needs_ridge() {
        // m < n makes Z Zᵀ singular (rank ≤ m); a zero ridge must fail, and a
        // positive ridge must rescue it. (A small ridge can already suffice,
        // because shifting zero eigenvalues by ε > roundoff makes the matrix
        // numerically PD — so the guarantee tested is *recovery*, not a fixed
        // escalation count.)
        let d = datagen::generate(80, 20, 4);
        let zero_ridge = Grm::build(&d.z, d.s, 0.0);
        assert!(
            zero_ridge.cholesky_lower().is_err(),
            "rank-deficient G must fail Cholesky at ridge 0"
        );

        let (grm, l, tries) = build_and_factor(&d.z, d.s, 1e-6, 12).unwrap();
        assert!(grm.ridge > 0.0);
        // Bookkeeping: zero escalations ⟺ the initial ridge was kept.
        if tries == 0 {
            assert!((grm.ridge - 1e-6).abs() < 1e-18);
        } else {
            assert!(grm.ridge > 1e-6);
        }
        // The returned factor must reconstruct the ridged matrix it factored.
        for i in 0..80 {
            for j in 0..=i {
                let mut acc = 0.0;
                for t in 0..=j {
                    acc += l[(i, t)] * l[(j, t)];
                }
                assert!((acc - grm.g[(i, j)]).abs() < 1e-7);
            }
        }
    }
}
