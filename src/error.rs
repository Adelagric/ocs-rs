//! Crate error type. Hand-written (no `thiserror` dependency) because the error
//! surface of a spike is tiny and stable.

use std::fmt;

/// Errors that can arise while building or solving an OCS instance.
#[derive(Debug)]
pub enum OcsError {
    /// Cholesky hit a non-positive pivot at `index`; the ridge was too small.
    /// `ridge` is the value that failed.
    NonPositivePivot { index: usize, ridge: f64 },
    /// Ridge escalation reached its cap without yielding a positive-definite
    /// matrix. `ridge` is the largest value tried.
    RidgeExhausted { ridge: f64 },
    /// Clarabel rejected the problem data (a constructor-level error, i.e. a
    /// bug in our assembly rather than a numerical outcome).
    SolverInit(String),
}

impl fmt::Display for OcsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OcsError::NonPositivePivot { index, ridge } => write!(
                f,
                "Cholesky non-positive pivot at index {index} with ridge {ridge:e}"
            ),
            OcsError::RidgeExhausted { ridge } => {
                write!(f, "ridge escalation exhausted (largest tried: {ridge:e})")
            }
            OcsError::SolverInit(msg) => write!(f, "Clarabel rejected problem data: {msg}"),
        }
    }
}

impl std::error::Error for OcsError {}
