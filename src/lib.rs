//! Optimum Contribution Selection (OCS) as a second-order cone program, solved
//! with the Clarabel interior-point method.
//!
//! This crate is a **go/no-go spike**. It exists to measure whether Clarabel
//! solves genomic-scale OCS reliably and quickly, given that the conic block it
//! must factor is effectively dense. See `README.md` for the formulation and
//! `REPORT.md` for the verdict with numbers.
//!
//! Pipeline: [`datagen`] → [`grm`] (Route A factor) → [`socp`] (cone assembly)
//! → [`solve`] (Clarabel) → [`report`] (feasibility + sweeps + verdict).

pub mod datagen;
pub mod error;
pub mod grm;
pub mod report;
pub mod socp;
pub mod solve;

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
