//! Where does the time in a matrix-free `G·c` product actually go?
//!
//! The product is `ε·c + Z(Zᵀc)/s`. Both halves are O(n·m) and both stream `Z`, so
//! the naive model says they cost the same and skipping one should halve the product.
//! This measures the three pieces separately instead of assuming that, across two
//! aspect ratios — a wide panel (n < m, like the mouse data) and a tall one (n > m).
//!
//! Run: cargo run --release --example bench_matvec

use faer::Mat;
use ocs_rs::datagen;
use std::hint::black_box;
use std::time::Instant;

/// `Zᵀc` as a dense GEMV — reads every entry of `Z`.
fn zt_c_dense(z: &Mat<f64>, c: &[f64]) -> Mat<f64> {
    let cm = Mat::from_fn(z.nrows(), 1, |i, _| c[i]);
    z.transpose() * cm.as_ref()
}

/// `Zᵀc` gathered from the nonzero rows of `c` only.
fn zt_c_sparse(z: &Mat<f64>, c: &[f64], nz: &[usize]) -> Mat<f64> {
    let m = z.ncols();
    let mut w = Mat::<f64>::zeros(m, 1);
    for l in 0..m {
        let mut acc = 0.0;
        for &i in nz {
            acc += c[i] * z[(i, l)];
        }
        w[(l, 0)] = acc;
    }
    w
}

fn time<T>(label: &str, reps: u32, mut f: impl FnMut() -> T) -> f64 {
    // One warm-up pass so the first run's page faults are not charged to the timing.
    black_box(f());
    let t0 = Instant::now();
    for _ in 0..reps {
        black_box(f());
    }
    let per = t0.elapsed().as_secs_f64() / reps as f64;
    println!("    {label:34} {:8.2} ms", per * 1000.0);
    per
}

fn main() {
    let reps = 10;
    for (n, m, shape) in [
        (1814, 10346, "wide  (n < m, mouse-like)"),
        (16000, 2000, "tall  (n > m)"),
    ] {
        let data = datagen::generate(n, m, 7);
        let z = &data.z;
        // A support-sized iterate: 20 nonzeros, as the active set produces.
        let mut c = vec![0.0; n];
        let nz: Vec<usize> = (0..20).map(|j| j * (n / 20)).collect();
        for &i in &nz {
            c[i] = 1.0 / nz.len() as f64;
        }
        let bytes = (n * m * 8) as f64 / (1024.0 * 1024.0);
        println!(
            "\n{shape}: n={n} m={m}, Z = {bytes:.0} MiB, |support| = {}",
            nz.len()
        );

        let t_dense = time("Zᵀc  dense GEMV", reps, || zt_c_dense(z, &c));
        let t_sparse = time("Zᵀc  sparse gather", reps, || zt_c_sparse(z, &c, &nz));
        let w = zt_c_dense(z, &c);
        let t_zw = time("Z·w  (irreducible: prices all n)", reps, || z * w.as_ref());

        // Sanity: the two Zᵀc routes must agree, or the comparison is meaningless.
        let ws = zt_c_sparse(z, &c, &nz);
        let diff = (0..m)
            .map(|l| (w[(l, 0)] - ws[(l, 0)]).abs())
            .fold(0.0, f64::max);
        println!("    max |dense - sparse| = {diff:.2e}");
        println!(
            "    product: dense route {:.2} ms -> sparse route {:.2} ms  ({:.2}x)",
            (t_dense + t_zw) * 1000.0,
            (t_sparse + t_zw) * 1000.0,
            (t_dense + t_zw) / (t_sparse + t_zw)
        );
        println!(
            "    share of the product spent in Zᵀc: {:.0}%",
            100.0 * t_dense / (t_dense + t_zw)
        );
    }
}
