//! A PLINK panel to an OCS optimum without ever forming a matrix.
//!
//! `.bed` stores genotypes 2-bit, 4 per byte, column-major — the layout the solver's
//! packed store uses. `plink::read_packed` therefore remaps the file's bytes and
//! solves from them, while `plink::read_panel` expands the same file to a dense `f64`
//! `Z`, 32× larger. This prints both costs and checks the optimum is the same one.
//!
//!   cargo run --release --example plink_packed                 # generated pig-shape panel
//!   cargo run --release --example plink_packed -- <prefix>     # your own trio
//!   cargo run --release --example plink_packed -- <prefix> packed   # one route, for `time -l`
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use ocs_rs::plink;
use ocs_rs::support_first::{solve_sexed, solve_sexed_packed};

const RIDGE: f64 = 1e-5;

/// Write a synthetic trio, so the example runs with no data to hand. Genotypes are
/// HWE draws at per-marker frequencies, with a realistic sprinkling of missing calls.
fn generate_trio(dir: &Path, n: usize, m: usize, miss_rate: f64) -> std::io::Result<PathBuf> {
    let d = ocs_rs::datagen::generate(n, m, 20240617);
    std::fs::create_dir_all(dir)?;
    let prefix = dir.join("panel");
    let mut fam = std::fs::File::create(prefix.with_extension("fam"))?;
    for i in 0..n {
        writeln!(
            fam,
            "FAM{i} IND{i} 0 0 {} -9",
            if i % 2 == 0 { 1 } else { 2 }
        )?;
    }
    let mut bim = std::fs::File::create(prefix.with_extension("bim"))?;
    for j in 0..m {
        writeln!(bim, "1 rs{j} 0 {} A G", (j + 1) * 100)?;
    }
    let mut bed = std::fs::File::create(prefix.with_extension("bed"))?;
    bed.write_all(&[0x6c, 0x1b, 0x01])?;
    let row_bytes = n.div_ceil(4);
    // A cheap deterministic hash decides which calls go missing, so the file is
    // reproducible without threading an RNG through.
    let missing = |i: usize, j: usize| {
        let h = (i.wrapping_mul(2654435761) ^ j.wrapping_mul(40503)) % 10_000;
        (h as f64) < miss_rate * 10_000.0
    };
    for j in 0..m {
        let mut row = vec![0u8; row_bytes];
        for i in 0..n {
            let dosage = (d.z[(i, j)] + 2.0 * d.p[j]).round().clamp(0.0, 2.0) as u8;
            let code = if missing(i, j) {
                0b01
            } else {
                match dosage {
                    2 => 0b00,
                    1 => 0b10,
                    _ => 0b11,
                }
            };
            row[i >> 2] |= code << ((i & 3) * 2);
        }
        bed.write_all(&row)?;
    }
    Ok(prefix)
}

fn mb(bytes: usize) -> f64 {
    bytes as f64 / 1e6
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let route = args.get(1).cloned().unwrap_or_else(|| "both".to_string());
    let tmp;
    let prefix = match args.first() {
        Some(p) => PathBuf::from(p),
        None => {
            let (n, m) = (5000usize, 50_000usize); // a pig-scale panel
            tmp = std::env::temp_dir().join("ocs_rs_plink_packed_demo");
            eprintln!(
                "generating a {n}×{m} panel in {} (1% missing)...",
                tmp.display()
            );
            generate_trio(&tmp, n, m, 0.01)?
        }
    };
    let bed_bytes = std::fs::metadata(prefix.with_extension("bed"))?.len() as usize;
    println!("{}.bed: {:.1} MB on disk", prefix.display(), mb(bed_bytes));

    let mut dense_out = None;
    let mut packed_out = None;
    let mut n = 0usize;

    if route != "packed" {
        let t = Instant::now();
        let panel = plink::read_panel(&prefix)?;
        let t_read = t.elapsed().as_secs_f64();
        n = panel.n;
        let m = panel.m;
        let held = n * m * 8;
        println!(
            "  read_panel  : dense Z {:>9.1} MB ({:.0}× the file)  read {t_read:.2}s",
            mb(held),
            held as f64 / bed_bytes as f64
        );
        let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let b: Vec<f64> = (0..n)
            .map(|i| ((i * 2654435761) % 1000) as f64 / 1000.0)
            .collect();
        let mean_diag: f64 = (0..n)
            .map(|i| (0..m).map(|l| panel.z[(i, l)].powi(2)).sum::<f64>() / panel.s)
            .sum::<f64>()
            / n as f64;
        let k = 0.05 * mean_diag;
        let t = Instant::now();
        let out = solve_sexed(panel.z.as_ref(), panel.s, RIDGE, &b, &male, k, 50_000, 1e-9);
        println!(
            "                solve {:.2}s  |S|={}  gain {:.9}  {:?}",
            t.elapsed().as_secs_f64(),
            out.support.len(),
            out.gain,
            out.status
        );
        dense_out = Some((out, k));
    }

    if route != "dense" {
        let t = Instant::now();
        let panel = plink::read_packed(&prefix)?;
        let t_read = t.elapsed().as_secs_f64();
        n = panel.geno.n();
        let held = panel.geno.packed_bytes();
        println!(
            "  read_packed : packed  {:>9.1} MB ({:.2}× the file)  read {t_read:.2}s",
            mb(held),
            held as f64 / bed_bytes as f64
        );
        let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let b: Vec<f64> = (0..n)
            .map(|i| ((i * 2654435761) % 1000) as f64 / 1000.0)
            .collect();
        // The cap the dense route used, or the same rule computed matrix-free.
        let k = match &dense_out {
            Some((_, k)) => *k,
            None => {
                let diag: f64 = (0..n)
                    .map(|i| panel.geno.gram_entry(i, i, 0.0))
                    .sum::<f64>()
                    / n as f64;
                0.05 * diag
            }
        };
        let t = Instant::now();
        let out = solve_sexed_packed(&panel.geno, RIDGE, &b, &male, k, 50_000, 1e-9);
        println!(
            "                solve {:.2}s  |S|={}  gain {:.9}  {:?}",
            t.elapsed().as_secs_f64(),
            out.support.len(),
            out.gain,
            out.status
        );
        packed_out = Some(out);
    }

    if let (Some((dense, _)), Some(packed)) = (&dense_out, &packed_out) {
        println!(
            "  same optimum: support match {}  Δgain {:.2e}",
            dense.support == packed.support,
            (dense.gain - packed.gain).abs()
        );
    }
    println!(
        "  (a dense n×n G for this panel would be {:.1} MB)",
        mb(n * n * 8)
    );
    Ok(())
}
