//! Python bindings for the support-first genomic OCS solver.
//!
//! The native layer deliberately speaks **flat `f64` data**, not numpy: the thin
//! wrapper in `python/ocs_rs/__init__.py` accepts arrays and flattens them. That
//! keeps this crate's dependency surface to pyo3 + faer + the core solver, with no
//! version coupling to a numpy-interop crate.
//!
//! Each function mirrors one core entry point and returns a dict with the
//! contributions, the active support, the gain, the kinship reached, and status.

use faer::Mat;
use ocs_rs::support_first::{self, SfStatus, SupportFirstOutcome};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Rebuild the n×m genotype matrix from a row-major flat slice.
fn mat_from_flat(z: &[f64], n: usize, m: usize) -> PyResult<Mat<f64>> {
    let want = n
        .checked_mul(m)
        .ok_or_else(|| PyValueError::new_err("n*m overflows"))?;
    if z.len() != want {
        return Err(PyValueError::new_err(format!(
            "genotype matrix has {} entries, expected n*m = {want}",
            z.len()
        )));
    }
    Ok(Mat::from_fn(n, m, |i, j| z[i * m + j]))
}

fn check_len(name: &str, got: usize, want: usize) -> PyResult<()> {
    if got == want {
        Ok(())
    } else {
        Err(PyValueError::new_err(format!(
            "{name} has length {got}, expected {want}"
        )))
    }
}

fn to_dict(py: Python<'_>, o: SupportFirstOutcome) -> PyResult<Bound<'_, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("c", o.c)?;
    d.set_item("support", o.support)?;
    d.set_item("gain", o.gain)?;
    d.set_item("quad", o.quad)?;
    d.set_item("iterations", o.iterations)?;
    d.set_item("products", o.products)?;
    d.set_item(
        "status",
        match o.status {
            SfStatus::Solved => "Solved",
            SfStatus::MaxIter => "MaxIter",
        },
    )?;
    Ok(d)
}

/// Simplex OCS: maximise bᵀc s.t. Σc = 1, c ≥ 0, cᵀGc ≤ k.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve<'py>(
    py: Python<'py>,
    z: Vec<f64>,
    n: usize,
    m: usize,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    k: f64,
    max_iter: u32,
    tol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    check_len("b", b.len(), n)?;
    let zm = mat_from_flat(&z, n, m)?;
    let out = py.detach(|| support_first::solve(&zm, s, ridge, &b, k, max_iter, tol));
    to_dict(py, out)
}

/// Sexed OCS: the same with Σ_males = Σ_females = ½.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, male, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_sexed<'py>(
    py: Python<'py>,
    z: Vec<f64>,
    n: usize,
    m: usize,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    male: Vec<bool>,
    k: f64,
    max_iter: u32,
    tol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    check_len("b", b.len(), n)?;
    check_len("male", male.len(), n)?;
    let zm = mat_from_flat(&z, n, m)?;
    let out = py.detach(|| support_first::solve_sexed(&zm, s, ridge, &b, &male, k, max_iter, tol));
    to_dict(py, out)
}

/// Simplex OCS with per-candidate caps 0 ≤ c ≤ u.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, caps, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_capped<'py>(
    py: Python<'py>,
    z: Vec<f64>,
    n: usize,
    m: usize,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    caps: Vec<f64>,
    k: f64,
    max_iter: u32,
    tol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    check_len("b", b.len(), n)?;
    check_len("caps", caps.len(), n)?;
    let zm = mat_from_flat(&z, n, m)?;
    let out = py.detach(|| support_first::solve_capped(&zm, s, ridge, &b, &caps, k, max_iter, tol));
    to_dict(py, out)
}

/// Sexed OCS with per-candidate caps.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, male, caps, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_sexed_capped<'py>(
    py: Python<'py>,
    z: Vec<f64>,
    n: usize,
    m: usize,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    male: Vec<bool>,
    caps: Vec<f64>,
    k: f64,
    max_iter: u32,
    tol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    check_len("b", b.len(), n)?;
    check_len("male", male.len(), n)?;
    check_len("caps", caps.len(), n)?;
    let zm = mat_from_flat(&z, n, m)?;
    let out = py.detach(|| {
        support_first::solve_sexed_capped(&zm, s, ridge, &b, &male, &caps, k, max_iter, tol)
    });
    to_dict(py, out)
}

#[pymodule]
fn _ocs_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(solve_sexed, m)?)?;
    m.add_function(wrap_pyfunction!(solve_capped, m)?)?;
    m.add_function(wrap_pyfunction!(solve_sexed_capped, m)?)?;
    Ok(())
}
