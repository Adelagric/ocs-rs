//! Throwaway-but-kept API spike. Exercises every external API the real crate
//! depends on, in one place, against a problem whose answer is known by hand.
//!
//! Hand solution (n=2, G=I, b=[1,2], k=0.6): maximise b·c s.t. sum c = 1,
//! c >= 0, ||c||_2^2 <= 0.6. Substituting c1 = 1 - c2 gives c2 in the binding
//! case c2 = (1 + sqrt(0.2)) / 2 = 0.723606..., c1 = 0.276393..., gain 1.7236.
//!
//! Run: cargo run --release --example spike

use clarabel::algebra::CscMatrix;
use clarabel::solver::{
    DefaultSettingsBuilder, DefaultSolver, IPSolver, NonnegativeConeT, SecondOrderConeT,
    SolverStatus, ZeroConeT,
};
use faer::{Mat, Side};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Binomial, Distribution, Normal, Uniform};

fn main() {
    // --- faer: build a small SPD matrix, matmul, and Cholesky ---------------
    let z = Mat::<f64>::from_fn(3, 2, |i, j| (i + 2 * j) as f64 + 1.0);
    let g = &z * z.transpose(); // 3x3 PSD (rank 2) -> ridge before llt
    let mut g_ridge = g.clone();
    for i in 0..3 {
        g_ridge[(i, i)] += 1e-6;
    }
    let llt = g_ridge.llt(Side::Lower).expect("ridge keeps it PD");
    let l = llt.L();
    // Reconstruct lower triangle of L*L^T and compare to g_ridge.
    let mut max_err = 0.0_f64;
    for i in 0..3 {
        for j in 0..=i {
            let mut acc = 0.0;
            for t in 0..=j {
                acc += l[(i, t)] * l[(j, t)];
            }
            max_err = max_err.max((acc - g_ridge[(i, j)]).abs());
        }
    }
    println!("faer: L L^T reconstruction max abs err = {max_err:.3e}");
    assert!(max_err < 1e-9, "Cholesky reconstruction failed");

    // --- rand 0.10 / rand_distr 0.6: seeded sampling ------------------------
    let mut rng = StdRng::seed_from_u64(42);
    let unif = Uniform::new(0.05, 0.5).unwrap();
    let p: f64 = unif.sample(&mut rng);
    let binom = Binomial::new(2, p).unwrap();
    let dosage: u64 = binom.sample(&mut rng);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let gebv: f64 = normal.sample(&mut rng);
    println!("rand: p={p:.4} dosage={dosage} gebv={gebv:.4}");
    assert!((0.05..0.5).contains(&p));
    assert!(dosage <= 2);

    // --- clarabel: the 2-variable OCS SOCP with a known optimum -------------
    // Layout matches the real socp module: rows = [eq(1); nonneg(n); soc(d+1)].
    let n = 2usize;
    let d = 2usize; // G = I  =>  F = L = I, SOC dimension d+1 = 3
    let k = 0.6_f64;
    let r = k.sqrt();
    let bvec = [1.0_f64, 2.0_f64];

    // Triplets (row, col, val). new_from_triplets sorts/consolidates for us;
    // the real crate builds CSC columns directly, but for the spike this is
    // the clearest way to pin the cone semantics.
    let mut ri = Vec::new();
    let mut ci = Vec::new();
    let mut vi = Vec::new();
    // eq row 0: all ones
    for j in 0..n {
        ri.push(0);
        ci.push(j);
        vi.push(1.0);
    }
    // nonneg rows 1..=n: -I
    for j in 0..n {
        ri.push(1 + j);
        ci.push(j);
        vi.push(-1.0);
    }
    // soc rows: r-row is (1+n); then -F^T in rows (2+n .. 2+n+d). F = I here.
    for j in 0..n {
        ri.push(2 + n + j);
        ci.push(j);
        vi.push(-1.0);
    }
    let a = CscMatrix::new_from_triplets(1 + n + (d + 1), n, ri, ci, vi);

    let mut bfull = vec![1.0]; // eq rhs
    bfull.extend(std::iter::repeat_n(0.0, n)); // nonneg rhs
    bfull.push(r); // soc r
    bfull.extend(std::iter::repeat_n(0.0, d)); // soc tail

    let q: Vec<f64> = bvec.iter().map(|x| -x).collect();
    let p_mat = CscMatrix::<f64>::zeros((n, n));

    let cones = [ZeroConeT(1), NonnegativeConeT(n), SecondOrderConeT(d + 1)];
    let settings = DefaultSettingsBuilder::default()
        .verbose(false)
        .build()
        .unwrap();

    let mut solver = DefaultSolver::new(&p_mat, &q, &a, &bfull, &cones, settings).unwrap();
    solver.solve();

    let x = &solver.solution.x;
    println!(
        "clarabel: status={:?} iters={} time={:.4}ms x=[{:.6}, {:.6}] gain={:.6}",
        solver.solution.status,
        solver.solution.iterations,
        solver.solution.solve_time * 1e3,
        x[0],
        x[1],
        -solver.solution.obj_val,
    );

    assert_eq!(solver.solution.status, SolverStatus::Solved);
    let c2_expected = (1.0 + 0.2_f64.sqrt()) / 2.0; // 0.7236068
    assert!(
        (x[1] - c2_expected).abs() < 1e-5,
        "expected c2 ~= {c2_expected:.6}, got {:.6}",
        x[1]
    );
    // Feasibility from the ORIGINAL data (not solver internals): the real check.
    let sum: f64 = x.iter().sum();
    let quad: f64 = x.iter().map(|v| v * v).sum(); // c^T I c
    assert!((sum - 1.0).abs() < 1e-6, "sum c = {sum}");
    assert!(x.iter().all(|&v| v > -1e-7), "nonneg violated");
    assert!(quad <= k + 1e-6, "kinship violated: {quad} > {k}");

    println!("\nSPIKE PASSED: faer + rand + clarabel compose and the SOCP optimum is correct.");
}
