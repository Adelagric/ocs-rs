//! Exact, matrix-free Optimum Contribution Selection (OCS) at genomic scale.
//!
//! The contribution is [`support_first`]: an active-set / column-generation solver
//! that exploits the tiny support of the OCS optimum, solves each fixed support in
//! closed form, and never forms the dense `n×n` relationship matrix (matrix-free
//! `G·c` from the genotype matrix `Z`). It reaches the same optimum as a conic
//! interior-point solver orders of magnitude faster; the sexed constraints
//! (`Σ_males = Σ_females = ½`) and per-candidate caps `c ≤ u` are supported.
//!
//! The project began as a go/no-go on the [Clarabel](https://clarabel.org) conic
//! solver (verdict **GO**, see `REPORT.md`); Clarabel is retained as an independent
//! cross-check oracle — [`socp`] assembles the cone, [`solve`] runs it.
//!
//! Pipeline: [`datagen`] → [`grm`] → {[`support_first`] | [`socp`] + [`solve`]} → [`report`].

pub mod datagen;
pub mod error;
pub mod grm;
pub mod report;
pub mod socp;
pub mod solve;
pub mod support_first;

pub use error::OcsError;

/// Feasibility of a returned contribution vector, checked against the original
/// problem data (never solver internals).
#[derive(Clone, Copy, Debug)]
pub struct Feasibility {
    /// `|Σcᵢ − 1|`.
    pub sum_err: f64,
    /// Most negative entry (0 if all non-negative).
    pub min_c: f64,
    /// Kinship quadratic form `cᵀGc`.
    pub quad: f64,
    /// Kinship bound `k`.
    pub k: f64,
}

impl Feasibility {
    /// True when every constraint holds within `tol`.
    pub fn ok(&self, tol: f64) -> bool {
        self.sum_err <= tol && self.min_c >= -tol && self.quad <= self.k + tol
    }
}

/// Assess feasibility from the contribution vector and the kinship form.
///
/// `quad` is `cᵀGc` computed by the caller via whichever route is in memory
/// ([`grm::quad_form`] for Route A, [`grm::quad_form_z`] for Route B).
pub fn feasibility(c: &[f64], quad: f64, k: f64) -> Feasibility {
    let sum: f64 = c.iter().sum();
    let min_c = c.iter().copied().fold(f64::INFINITY, f64::min);
    Feasibility {
        sum_err: (sum - 1.0).abs(),
        min_c: if min_c.is_finite() { min_c } else { 0.0 },
        quad,
        k,
    }
}
