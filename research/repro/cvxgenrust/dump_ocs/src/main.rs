//! Dump the Table-1 Route-A instance (L, b, k) for a third-party generated Clarabel
//! solver, and time the crate's own Route-A conic path on the same data.
//!
//!   cargo run --release -- <n> [m=20000] [seed=20240617] [k_frac=0.017]
//!
//! The instance is the one of the manuscript's Table 1 (`compare --n N --m 20000
//! --k-frac 0.017`, seed 20240617); `L` is written column-major, CVXPY's flattened
//! order, which the generated setters expect.
use std::io::Write;
use std::time::Instant;

use ocs_rs::datagen;
use ocs_rs::grm::Grm;
use ocs_rs::socp::{self, Factor};
use ocs_rs::solve::{self, SolveConfig};
use ocs_rs::support_first;

const RIDGE: f64 = 1e-5;

/// Same definition as `src/main.rs::mean_diag` (ridge excluded).
fn mean_diag(z: faer::MatRef<'_, f64>, s: f64) -> f64 {
    let (n, m) = (z.nrows(), z.ncols());
    let acc: f64 = (0..n)
        .map(|i| (0..m).map(|j| z[(i, j)] * z[(i, j)]).sum::<f64>() / s)
        .sum();
    acc / n as f64
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(1000);
    let m: usize = args.get(2).map(|s| s.parse()).transpose()?.unwrap_or(20000);
    let seed: u64 = args.get(3).map(|s| s.parse()).transpose()?.unwrap_or(20240617);
    let k_frac: f64 = args.get(4).map(|s| s.parse()).transpose()?.unwrap_or(0.017);

    let ds = datagen::generate(n, m, seed);
    let k = k_frac * mean_diag(ds.z.as_ref(), ds.s);
    let grm = Grm::build(ds.z.as_ref(), ds.s, RIDGE);
    let l = grm.cholesky_lower()?;

    // --- dump: header (n, k), b, then L column-major (CVXPY flattened order) ---
    let path = format!("ocs_n{n}.bin");
    let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
    f.write_all(&(n as u64).to_le_bytes())?;
    f.write_all(&k.to_le_bytes())?;
    for &v in &ds.b {
        f.write_all(&v.to_le_bytes())?;
    }
    for j in 0..n {
        for i in 0..n {
            f.write_all(&l[(i, j)].to_le_bytes())?;
        }
    }
    drop(f);
    eprintln!("[dump] wrote {path}  n={n} m={m} seed={seed} k={k:.6}");

    // --- the crate's own Route A on the identical data (cone assembly + IPM) ---
    let t = Instant::now();
    let prob = socp::build(Factor::Cholesky(&l), &ds.b, k, ds.s, None);
    let t_assemble = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let out = solve::solve(&prob, SolveConfig::default())?;
    let t_solve = t.elapsed().as_secs_f64();
    let supp = out.c.iter().filter(|&&c| c > 1e-6).count();
    println!(
        "ocs-rs RouteA  n={n}: assemble {t_assemble:.3}s  solve(wall) {t_solve:.3}s  \
         clarabel.solve_time {:.3}s  iters {}  gain {:.9}  |S|={supp}  status {:?}",
        out.solve_time, out.iterations, out.gain, out.status
    );

    // --- support-first, same instance, for the record ---
    let t = Instant::now();
    let sf = support_first::solve(ds.z.as_ref(), ds.s, RIDGE, &ds.b, k, 4000, 1e-7);
    let t_sf = t.elapsed().as_secs_f64();
    println!(
        "support-first  n={n}: {t_sf:.3}s  gain {:.9}  |S|={}  status {:?}  Δgain vs RouteA {:.2e}",
        sf.gain,
        sf.support.len(),
        sf.status,
        (sf.gain - out.gain).abs()
    );
    Ok(())
}
