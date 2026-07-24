//! Python bindings for the support-first genomic OCS solver.
//!
//! The native layer deliberately speaks **flat `f64` data**, not numpy: the thin
//! wrapper in `python/ocs_rs/__init__.py` accepts arrays and flattens them. That
//! keeps this crate's dependency surface to pyo3 + faer + the core solver, with no
//! version coupling to a numpy-interop crate.
//!
//! The genotype matrix — the only large argument — arrives as raw `bytes` from
//! `Z.tobytes()`, not as a `Vec<f64>`. Extracting a `Vec` iterates the array element
//! by element through the interpreter, which was measured as the dominant cost of a
//! call once the solver itself became milliseconds (≈700 ms at n=16000). `tobytes()`
//! is a single C-level `memcpy`, and pyo3 hands the result to Rust as a borrowed
//! `&[u8]` with no further copy. (The obvious alternative, pyo3's `PyBuffer`, is gated
//! on Python ≥ 3.11 under the stable ABI, and this wheel is `abi3-py39` to reach older
//! interpreters — so the buffer protocol is not available here.) The small length-`n`
//! vectors (`b`, `caps`, `male`) stay `Vec`, where an O(n) copy is negligible against
//! the O(n·m) matrix.
//!
//! Each function mirrors one core entry point and returns a dict with the
//! contributions, the active support, the gain, the kinship reached, and status.

use faer::Mat;
use ocs_rs::support_first::{self, SfStatus, SupportFirstOutcome};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};

/// Build the n×m column-major genotype matrix from row-major native-endian `f64`
/// bytes (`Z.tobytes()` on a C-contiguous `float64` array).
///
/// numpy is row-major and faer column-major, so a gather is unavoidable; the point is
/// that it is the *only* copy on the Rust side — the bytes are borrowed, and each
/// `f64` is read once, straight into faer's layout, with no per-element excursion into
/// Python. `copy_from_slice` from a length-8 window into an 8-byte array cannot fail,
/// so there is no `unwrap` hiding a panic.
fn mat_from_bytes(z: &[u8], n: usize, m: usize) -> PyResult<Mat<f64>> {
    let want = n
        .checked_mul(m)
        .and_then(|e| e.checked_mul(8))
        .ok_or_else(|| PyValueError::new_err("n*m*8 overflows"))?;
    if z.len() != want {
        return Err(PyValueError::new_err(format!(
            "genotype buffer is {} bytes, expected n*m*8 = {want}",
            z.len()
        )));
    }
    Ok(Mat::from_fn(n, m, |i, j| {
        let o = (i * m + j) * 8;
        let mut bytes = [0u8; 8];
        bytes.copy_from_slice(&z[o..o + 8]);
        f64::from_ne_bytes(bytes)
    }))
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
    z: &[u8],
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
    let zm = mat_from_bytes(z, n, m)?;
    let out = py.detach(|| support_first::solve(&zm, s, ridge, &b, k, max_iter, tol));
    to_dict(py, out)
}

/// Sexed OCS: the same with Σ_males = Σ_females = ½.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, male, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_sexed<'py>(
    py: Python<'py>,
    z: &[u8],
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
    let zm = mat_from_bytes(z, n, m)?;
    let out = py.detach(|| support_first::solve_sexed(&zm, s, ridge, &b, &male, k, max_iter, tol));
    to_dict(py, out)
}

/// Simplex OCS with per-candidate caps 0 ≤ c ≤ u.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, caps, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_capped<'py>(
    py: Python<'py>,
    z: &[u8],
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
    let zm = mat_from_bytes(z, n, m)?;
    let out = py.detach(|| support_first::solve_capped(&zm, s, ridge, &b, &caps, k, max_iter, tol));
    to_dict(py, out)
}

/// Sexed OCS with per-candidate caps.
#[pyfunction]
#[pyo3(signature = (z, n, m, s, ridge, b, male, caps, k, max_iter=10_000, tol=1e-9))]
#[allow(clippy::too_many_arguments)]
fn solve_sexed_capped<'py>(
    py: Python<'py>,
    z: &[u8],
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
    let zm = mat_from_bytes(z, n, m)?;
    let out = py.detach(|| {
        support_first::solve_sexed_capped(&zm, s, ridge, &b, &male, &caps, k, max_iter, tol)
    });
    to_dict(py, out)
}

/// Read a PLINK 1 trio (`<prefix>.bed`/`.bim`/`.fam`) into centred genotypes.
///
/// The genotypes come back as row-major native-endian `f64` `bytes`, the mirror of
/// the input path: the Python layer wraps them with `np.frombuffer`, so a real panel
/// (potentially tens of millions of entries) is never materialised as a Python list.
#[pyfunction]
fn read_plink<'py>(py: Python<'py>, prefix: &str) -> PyResult<Bound<'py, PyDict>> {
    let panel = py
        .detach(|| ocs_rs::plink::read_panel(std::path::Path::new(prefix)))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (n, m) = (panel.n, panel.m);
    let mut raw = Vec::with_capacity(n * m * 8);
    for i in 0..n {
        for j in 0..m {
            raw.extend_from_slice(&panel.z[(i, j)].to_ne_bytes());
        }
    }
    let d = PyDict::new(py);
    d.set_item("z", PyBytes::new(py, &raw))?;
    d.set_item("n", n)?;
    d.set_item("m", m)?;
    d.set_item("p", panel.p)?;
    d.set_item("s", panel.s)?;
    d.set_item("ids", panel.ids)?;
    Ok(d)
}

#[pymodule]
fn _ocs_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(read_plink, m)?)?;
    m.add_function(wrap_pyfunction!(solve_sexed, m)?)?;
    m.add_function(wrap_pyfunction!(solve_capped, m)?)?;
    m.add_function(wrap_pyfunction!(solve_sexed_capped, m)?)?;
    Ok(())
}
