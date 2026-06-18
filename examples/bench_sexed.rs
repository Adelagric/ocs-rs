//! Release-profile timing of the sexed support-first solver on a genomic-scale
//! synthetic instance, plus an export so the Python prototype can solve the
//! *identical* problem for a cross-language exactness check.
//!
//!   cargo run --release --example bench_sexed
//!   python3 /tmp/bench_sf.py 2000      # solves the same instance, prints Δgain
//!
//! The Rust solver enforces `cᵀ(ZZᵀ/s + εI)c ≤ k`; we export that exact ridged
//! matrix as `K` and `ub = k`, so the Python prototype (which constrains
//! `cᵀKc ≤ ub`) faces the byte-for-byte same problem. Stdout stays clean for the
//! data export note; the timing/feasibility report goes to stderr.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use ocs_rs::datagen;
use ocs_rs::grm::Grm;
use ocs_rs::support_first::solve_sexed;

fn main() -> std::io::Result<()> {
    let (n, m, seed) = (2000usize, 10000usize, 20240618u64);
    let ridge = 1e-5;
    let d = datagen::generate(n, m, seed);
    let g = Grm::build(&d.z, d.s, ridge);
    let mean_diag: f64 = (0..n).map(|i| g.g[(i, i)]).sum::<f64>() / n as f64;
    // Alternating sexes (a balanced split; the value of `b` decides the optimum).
    let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();

    // Sweep the kinship cap from slack to tight so the support grows and the
    // column-generation path is actually exercised; report median of five timed
    // release runs at each. Export the binding frac = 0.15 instance for Python.
    eprintln!("solve_sexed sweep  n={n}  m={m}  ridge={ridge}  (release, median of 5 runs)");
    let mut export: Option<(f64, _)> = None;
    for &frac in &[0.6_f64, 0.3, 0.15, 0.08] {
        let k = frac * mean_diag;
        let mut last = solve_sexed(&d.z, d.s, ridge, &d.b, &male, k, 5000, 1e-7); // warm
        let mut times = Vec::with_capacity(5);
        for _ in 0..5 {
            let t = Instant::now();
            last = solve_sexed(&d.z, d.s, ridge, &d.b, &male, k, 5000, 1e-7);
            times.push(t.elapsed().as_secs_f64());
        }
        times.sort_by(f64::total_cmp);
        let med = times[times.len() / 2] * 1e3;
        let sum_m: f64 = (0..n).filter(|&i| male[i]).map(|i| last.c[i]).sum();
        eprintln!(
            "  frac={frac:>4}  k={k:.5}  gain={:.5}  quad={:.5}  support={:>3}  iters={:>3}  \
             Σmale={sum_m:.4}  median={med:.3} ms  {:?}",
            last.gain,
            last.quad,
            last.support.len(),
            last.iterations,
            last.status
        );
        if (frac - 0.15).abs() < 1e-9 {
            export = Some((k, last));
        }
    }

    let Some((k, last)) = export else {
        return Ok(());
    };
    // Export the identical (binding) instance for the Python prototype.
    let mut fk = BufWriter::new(File::create(format!("/tmp/bench_K_{n}.csv"))?);
    for i in 0..n {
        for j in 0..n {
            if j > 0 {
                fk.write_all(b",")?;
            }
            write!(fk, "{}", g.g[(i, j)])?;
        }
        fk.write_all(b"\n")?;
    }
    fk.flush()?;
    let mut fb = BufWriter::new(File::create(format!("/tmp/bench_bc_{n}.csv"))?);
    writeln!(fb, "\"bv\",\"oc\",\"sex\"")?;
    for (i, &is_male) in male.iter().enumerate() {
        let sex = if is_male { "male" } else { "female" };
        writeln!(fb, "{},{},\"{sex}\"", d.b[i], last.c[i])?;
    }
    fb.flush()?;
    File::create(format!("/tmp/bench_ub_{n}.txt"))?.write_all(format!("{k}").as_bytes())?;

    println!(
        "exported /tmp/bench_{{K,bc,ub}}_{n}.*  ->  cross-check: python3 /tmp/bench_sf.py {n}"
    );
    Ok(())
}
