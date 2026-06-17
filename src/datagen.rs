//! Seeded synthetic genotype and breeding-value generation.
//!
//! Reproducibility contract: a given `(n, m, seed)` triple always yields the
//! exact same data, because every random draw comes from a single
//! [`StdRng`](rand::rngs::StdRng) consumed in a fixed order:
//!
//! 1. for each marker column `j = 0..m`: draw the allele frequency `p[j]`, then
//!    the `n` dosages of that column (individual-major within the column);
//! 2. then the `n` breeding values `b[i]`.
//!
//! Genotypes are stored already centred as `Z = M - 2p`, because nothing in the
//! pipeline needs the raw dosages `M` once the relationship matrix is built.

use faer::Mat;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Binomial, Distribution, Normal, Uniform};

/// A synthetic breeding scenario.
pub struct Dataset {
    /// Number of candidate individuals.
    pub n: usize,
    /// Number of biallelic markers.
    pub m: usize,
    /// Column-centred genotype matrix `Z = M - 2p`, shape `n x m`.
    pub z: Mat<f64>,
    /// Allele frequencies, one per marker. Each in `[0.05, 0.5)`.
    pub p: Vec<f64>,
    /// Genomic estimated breeding values (genetic gain coefficients), length `n`.
    pub b: Vec<f64>,
    /// VanRaden scaling `s = 2 * sum_j p_j (1 - p_j)`.
    pub s: f64,
    /// Seed used, retained for provenance in dumped artifacts.
    pub seed: u64,
}

/// Generate a reproducible dataset of `n` individuals and `m` markers.
///
/// Allele frequencies are `Uniform(0.05, 0.5)`, dosages are `Binomial(2, p_j)`,
/// breeding values are `Normal(0, 1)`.
pub fn generate(n: usize, m: usize, seed: u64) -> Dataset {
    assert!(n > 0 && m > 0, "n and m must be positive");
    let mut rng = StdRng::seed_from_u64(seed);

    // The supports are constant and valid, so the distribution constructors
    // cannot fail here; the explicit panics document that invariant rather than
    // hiding it behind `?` in a function with no other error path.
    let freq = Uniform::new(0.05_f64, 0.5_f64).expect("0.05 < 0.5 is a valid uniform support");
    let gain = Normal::new(0.0_f64, 1.0_f64).expect("unit normal is valid");

    let mut p = vec![0.0_f64; m];
    let mut z = Mat::<f64>::zeros(n, m);
    let mut s = 0.0_f64;

    for j in 0..m {
        let pj = freq.sample(&mut rng);
        p[j] = pj;
        s += pj * (1.0 - pj);
        let dosage = Binomial::new(2, pj).expect("0 <= p <= 1 by construction");
        let centre = 2.0 * pj;
        for i in 0..n {
            let d = dosage.sample(&mut rng) as f64;
            z[(i, j)] = d - centre;
        }
    }
    s *= 2.0;

    let b: Vec<f64> = (0..n).map(|_| gain.sample(&mut rng)).collect();

    Dataset {
        n,
        m,
        z,
        p,
        b,
        s,
        seed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproducible_for_same_seed() {
        let a = generate(20, 50, 7);
        let b = generate(20, 50, 7);
        assert_eq!(a.b, b.b);
        assert_eq!(a.p, b.p);
        assert!((a.s - b.s).abs() < 1e-12);
        for i in 0..20 {
            for j in 0..50 {
                assert_eq!(a.z[(i, j)], b.z[(i, j)]);
            }
        }
    }

    #[test]
    fn different_seed_differs() {
        let a = generate(20, 50, 1);
        let b = generate(20, 50, 2);
        assert!(a.b != b.b);
    }

    #[test]
    fn columns_are_centred() {
        // Z = M - 2p, and E[M_ij] = 2 p_j, so each column of Z averages ~0.
        let d = generate(2000, 4, 123);
        for j in 0..d.m {
            let mean: f64 = (0..d.n).map(|i| d.z[(i, j)]).sum::<f64>() / d.n as f64;
            assert!(mean.abs() < 0.1, "column {j} mean {mean} not near 0");
        }
    }

    #[test]
    fn dosages_in_range() {
        // Z + 2p must reconstruct integer dosages in {0,1,2}.
        let d = generate(100, 10, 99);
        for j in 0..d.m {
            for i in 0..d.n {
                let dosage = d.z[(i, j)] + 2.0 * d.p[j];
                let rounded = dosage.round();
                assert!((dosage - rounded).abs() < 1e-9);
                assert!((0.0..=2.0).contains(&rounded));
            }
        }
    }

    #[test]
    fn vanraden_scale_positive() {
        let d = generate(10, 100, 5);
        assert!(d.s > 0.0);
    }
}
