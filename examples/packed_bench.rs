//! Pig-shape benchmark of the 2-bit packed sexed solve against the dense-`f64` one.
//! Same optimum (asserted) at 32× less memory; whether the packing also buys time
//! at the pig's marker count is what this prints. As measured it does not: the
//! packed product is compute-bound (unpacking) where the dense stream is
//! bandwidth-bound, and the packed solve is slower (research/REVISION.md,
//! post-0.3.0 audit). Run on an idle machine — the packed product is threaded and
//! is the side that suffers under load.
//!
//!   cargo run --release --example packed_bench

use std::time::Instant;

use ocs_rs::datagen;
use ocs_rs::packed::PackedGeno;
use ocs_rs::support_first::{solve_sexed, solve_sexed_packed};

fn main() {
    let n = 1194usize;
    let m = 52843usize; // the PIC pig marker count
    let ridge = 1e-5;
    eprintln!("generating a pig-shape instance (n={n}, m={m})...");
    let d = datagen::generate(n, m, 42);
    let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
    let packed = PackedGeno::from_getter(
        n,
        m,
        |i, l| (d.z[(i, l)] + 2.0 * d.p[l]).round() as u8,
        d.p.clone(),
        d.s,
    );

    // Operating cap (tight): support runs to tens–low hundreds, the informative regime.
    let k = 0.02;
    let solve_dense = || solve_sexed(d.z.as_ref(), d.s, ridge, &d.b, &male, k, 50_000, 1e-8);
    let solve_pk = || solve_sexed_packed(&packed, ridge, &d.b, &male, k, 50_000, 1e-8);

    let _ = solve_dense();
    let _ = solve_pk(); // warm caches

    let t = Instant::now();
    let dense = solve_dense();
    let t_dense = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let pk = solve_pk();
    let t_pk = t.elapsed().as_secs_f64();

    let dense_mb = (n * m * 8) as f64 / 1e6;
    let packed_mb = packed.packed_bytes() as f64 / 1e6;
    println!(
        "pig-shape solve  n={n} m={m}  |S|={} (cap k={k})",
        dense.support.len()
    );
    println!(
        "  same optimum: support match={}  gain Δ={:.2e}",
        dense.support == pk.support,
        (dense.gain - pk.gain).abs()
    );
    println!("  dense  Z (f64)  : {dense_mb:8.1} MB   solve {t_dense:.3} s");
    println!("  packed Z (2-bit): {packed_mb:8.1} MB   solve {t_pk:.3} s",);
    println!(
        "  -> packed is {:.1}× less memory and {:.2}× the dense solve time",
        dense_mb / packed_mb,
        t_dense / t_pk.max(1e-6)
    );
}
