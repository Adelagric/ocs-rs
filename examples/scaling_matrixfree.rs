//! Matrix-free vs dense-G scaling — the matrix-free justification.
//!
//! At m < n the dense n×n G costs O(n²) memory and O(n²·m) to build; the
//! matrix-free solver pays neither. This sweeps n at a fixed small m and shows
//! the dense G's memory and build time growing quadratically — to the point of
//! crossing the benchmark allocation policy — while the sexed matrix-free solve stays cheap
//! and the stored Z grows only linearly in n. This is the regime (large
//! populations, modest marker panel) where matrix-free is the enabler, not a
//! micro-optimisation.
//!
//!   cargo run --release --example scaling_matrixfree

use std::time::Instant;

use ocs_rs::datagen;
use ocs_rs::grm::Grm;
use ocs_rs::support_first::solve_sexed;

fn gib(elems: usize) -> f64 {
    elems as f64 * 8.0 / (1024.0 * 1024.0 * 1024.0)
}

fn main() {
    let m = 1000usize;
    let ridge = 1e-5;
    let build_ceiling = 30_000usize; // allocation policy: dense G through n=30000 (6.7 GiB)

    // Figure 1 reads these rows rather than carrying a hand-copied duplicate of them,
    // so regenerating the figure cannot silently lag a change to the solver.
    let output_dir = std::env::var("REPRO_OUTPUT_DIR").unwrap_or_else(|_| "artifacts".into());
    let csv_path = std::path::Path::new(&output_dir).join("scaling_matrixfree.csv");
    let mut csv = String::from("n,m,z_gib,g_gib,g_build_s,mfree_solve_s,support,status\n");

    println!("matrix-free vs dense-G   m={m}  (sexed OCS, release profile)");
    println!(
        "{:>7} {:>9} {:>9} {:>12} {:>14} {:>8} {:>8}",
        "n", "Z GiB", "G GiB", "G build s", "mfree solve s", "support", "status"
    );
    for &n in &[1000usize, 2000, 5000, 10000, 20000, 30000, 40000] {
        let d = datagen::generate(n, m, 20240618 + n as u64);
        // Mean GRM diagonal without forming G: G_ii = (z_i·z_i)/s + ridge.
        let mean_diag: f64 = (0..n)
            .map(|i| (0..m).map(|l| d.z[(i, l)].powi(2)).sum::<f64>() / d.s + ridge)
            .sum::<f64>()
            / n as f64;
        let k = 0.1 * mean_diag; // binding: forces the support to grow off [½,½]
        let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();

        // Matrix-free sexed solve: warm-up, then one timed release run.
        let _ = solve_sexed(d.z.as_ref(), d.s, ridge, &d.b, &male, k, 5000, 1e-7);
        let t = Instant::now();
        let last = solve_sexed(d.z.as_ref(), d.s, ridge, &d.b, &male, k, 5000, 1e-7);
        let best = t.elapsed().as_secs_f64();
        if format!("{:?}", last.status) != "Solved" {
            eprintln!("[scaling] n={n}: solve failed: {:?}", last.status);
            std::process::exit(1);
        }

        // Dense G: build + time only within this example's allocation policy.
        let g_build = if n <= build_ceiling {
            let t = Instant::now();
            let g = Grm::build(d.z.as_ref(), d.s, ridge);
            let secs = t.elapsed().as_secs_f64();
            std::hint::black_box(&g); // keep the build from being optimised away
            format!("{secs:.3}")
        } else {
            "not_measured".to_string()
        };

        println!(
            "{n:>7} {:>9.3} {:>9.3} {g_build:>12} {best:>14.4} {:>8} {:>8}",
            gib(n * m),
            gib(n * n),
            last.support.len(),
            format!("{:?}", last.status)
        );
        csv.push_str(&format!(
            "{n},{m},{:.4},{:.4},{g_build},{best:.4},{},{:?}\n",
            gib(n * m),
            gib(n * n),
            last.support.len(),
            last.status
        ));
    }

    if let Err(e) =
        std::fs::create_dir_all(&output_dir).and_then(|()| std::fs::write(&csv_path, &csv))
    {
        eprintln!("[scaling] could not write {}: {e}", csv_path.display());
        std::process::exit(1);
    } else {
        eprintln!("[scaling] wrote {}", csv_path.display());
    }
    println!(
        "G/Z GiB are theoretical float64 array footprints, excluding workspace. \
         Dense builds above n={build_ceiling} are not measured by allocation policy. \
         Build cost grows O(n²·m); this loose-cap solve is not the structured operating-point sweep."
    );
}
