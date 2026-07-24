//! R bindings for the support-first genomic OCS solver.
//!
//! These functions are the thin native layer; the user-facing API is the R code in
//! `R/ocs.R`. Two conventions are handled here rather than in R, where getting them
//! wrong is silent:
//!
//! - **Column-major.** R stores a matrix as a column-major vector, so element
//!   `(i, j)` is `z[j*n + i]`. Reading it in that order means neither side pays for
//!   a transpose.
//! - **1-based supports.** Indices leave this layer already incremented, so an R
//!   user never sees a 0-based index.

use extendr_api::prelude::*;
// `Result` is not in extendr's prelude, so the glob above leaves std's in scope; this
// explicit import is what makes `Result<List>` mean "or an R error".
use extendr_api::Result;
use faer::Mat;
use ocs_rs::support_first::{self, SfStatus, SupportFirstOutcome};

/// Rebuild the n×m genotype matrix from R's column-major vector.
fn mat_from_col_major(z: &[f64], n: usize, m: usize) -> Result<Mat<f64>> {
    let want = n
        .checked_mul(m)
        .ok_or_else(|| Error::Other("n * m overflows".to_string()))?;
    if z.len() != want {
        return Err(Error::Other(format!(
            "genotype matrix has {} entries, expected n * m = {want}",
            z.len()
        )));
    }
    Ok(Mat::from_fn(n, m, |i, j| z[j * n + i]))
}

fn check_len(name: &str, got: usize, want: usize) -> Result<()> {
    if got == want {
        Ok(())
    } else {
        Err(Error::Other(format!(
            "{name} has length {got}, expected {want}"
        )))
    }
}

fn as_bools(male: &[i32]) -> Vec<bool> {
    male.iter().map(|&x| x != 0).collect()
}

fn to_list(o: SupportFirstOutcome) -> List {
    let support: Vec<i32> = o.support.iter().map(|&i| i as i32 + 1).collect();
    list!(
        c = o.c,
        support = support,
        gain = o.gain,
        quad = o.quad,
        iterations = o.iterations as i32,
        products = o.products as i32,
        status = match o.status {
            SfStatus::Solved => "Solved",
            SfStatus::MaxIter => "MaxIter",
        }
    )
}

/// Simplex OCS: maximise bᵀc s.t. sum(c) = 1, c >= 0, cᵀGc <= k.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn sf_solve(
    z: Vec<f64>,
    n: i32,
    m: i32,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    k: f64,
    max_iter: i32,
    tol: f64,
) -> Result<List> {
    let (n, m) = (n as usize, m as usize);
    check_len("b", b.len(), n)?;
    let zm = mat_from_col_major(&z, n, m)?;
    Ok(to_list(support_first::solve(
        &zm,
        s,
        ridge,
        &b,
        k,
        max_iter as u32,
        tol,
    )))
}

/// Sexed OCS: the same with sum over males = sum over females = 1/2.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn sf_solve_sexed(
    z: Vec<f64>,
    n: i32,
    m: i32,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    male: Vec<i32>,
    k: f64,
    max_iter: i32,
    tol: f64,
) -> Result<List> {
    let (n, m) = (n as usize, m as usize);
    check_len("b", b.len(), n)?;
    check_len("male", male.len(), n)?;
    let zm = mat_from_col_major(&z, n, m)?;
    Ok(to_list(support_first::solve_sexed(
        &zm,
        s,
        ridge,
        &b,
        &as_bools(&male),
        k,
        max_iter as u32,
        tol,
    )))
}

/// Simplex OCS with per-candidate caps 0 <= c <= u.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn sf_solve_capped(
    z: Vec<f64>,
    n: i32,
    m: i32,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    caps: Vec<f64>,
    k: f64,
    max_iter: i32,
    tol: f64,
) -> Result<List> {
    let (n, m) = (n as usize, m as usize);
    check_len("b", b.len(), n)?;
    check_len("caps", caps.len(), n)?;
    let zm = mat_from_col_major(&z, n, m)?;
    Ok(to_list(support_first::solve_capped(
        &zm,
        s,
        ridge,
        &b,
        &caps,
        k,
        max_iter as u32,
        tol,
    )))
}

/// Sexed OCS with per-candidate caps.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn sf_solve_sexed_capped(
    z: Vec<f64>,
    n: i32,
    m: i32,
    s: f64,
    ridge: f64,
    b: Vec<f64>,
    male: Vec<i32>,
    caps: Vec<f64>,
    k: f64,
    max_iter: i32,
    tol: f64,
) -> Result<List> {
    let (n, m) = (n as usize, m as usize);
    check_len("b", b.len(), n)?;
    check_len("male", male.len(), n)?;
    check_len("caps", caps.len(), n)?;
    let zm = mat_from_col_major(&z, n, m)?;
    Ok(to_list(support_first::solve_sexed_capped(
        &zm,
        s,
        ridge,
        &b,
        &as_bools(&male),
        &caps,
        k,
        max_iter as u32,
        tol,
    )))
}

/// Read a PLINK 1 trio into centred genotypes, returned column-major for R.
#[extendr]
fn sf_read_plink(prefix: &str) -> Result<List> {
    let panel = ocs_rs::plink::read_panel(std::path::Path::new(prefix))
        .map_err(|e| Error::Other(e.to_string()))?;
    let (n, m) = (panel.n, panel.m);
    let mut z = Vec::with_capacity(n * m);
    for j in 0..m {
        for i in 0..n {
            z.push(panel.z[(i, j)]);
        }
    }
    Ok(list!(
        z = z,
        n = n as i32,
        m = m as i32,
        p = panel.p,
        s = panel.s,
        ids = panel.ids
    ))
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod ocsrs;
    fn sf_solve;
    fn sf_solve_sexed;
    fn sf_solve_capped;
    fn sf_solve_sexed_capped;
    fn sf_read_plink;
}
