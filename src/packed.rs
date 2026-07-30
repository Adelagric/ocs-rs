//! Genotypes packed at 2 bits each, with the VanRaden kinship products formed on the
//! fly — 32× less memory than the `f64` centred `Z`, reproducing the dense-`f64` path
//! to machine precision. This is what makes the matrix-free memory claim hold at
//! population scale (`n ≈ 5·10⁵`, `m = 5·10⁴`: `Z` packed 6.25 GB vs dense `G` 2 TB).
//!
//! Raw genotypes `M ∈ {0,1,2}` are stored 4-to-a-byte, column-major (matching the
//! solver's column-major layout). Centring, `Z = M − 2p`, is never materialised: it
//! factors out of both products the solver needs as cheap rank-1 corrections, using
//! the per-candidate scalars `r_i = p·M_i` and `‖p‖²` precomputed once.
//!
//!   Zᵀc      = Mᵀc − 2p·(1ᵀc)
//!   Z(Zᵀc)   = M·t − 2p·(pᵀt)                       (t = Zᵀc)
//!   zᵢ·zⱼ    = Mᵢ·Mⱼ − 2rᵢ − 2rⱼ + 4‖p‖²             (a Gram entry)
//!
//! so the hot loops touch only the packed raw genotypes. This module provides the two
//! primitives the active-set solver uses — the full product `G·c` and the restricted
//! Gram entry `G_S[i,j]` — verified against the `f64` reference in the tests. Wiring it
//! into `support_first` is a separate step; kept out here so the core cannot break.

/// Column-major 2-bit genotype store plus the scalars that let centring be applied on
/// the fly. `G = ZZᵀ/s + ridge·I` with `Z = M − 2p`, `M` the raw dosages.
pub struct PackedGeno {
    /// 2-bit genotypes, 4 per byte, column-major: candidate `i` of marker `j` lives in
    /// byte `j·bytes_per_col + i/4`, bits `(i%4)·2`.
    data: Vec<u8>,
    n: usize,
    m: usize,
    bytes_per_col: usize,
    /// Allele frequencies `p[j]` (the centring vector is `2p`).
    p: Vec<f64>,
    /// VanRaden scale `s = 2 Σ p(1−p)`.
    s: f64,
    /// Per-candidate `r_i = p · M_i = (M p)_i`, precomputed once for the Gram centring.
    r: Vec<f64>,
    /// `‖p‖²`, the constant in the Gram centring.
    p_sq: f64,
}

impl PackedGeno {
    /// Pack from any raw-genotype accessor `get(i, j) → {0,1,2}`, precomputing the
    /// centring scalars. `p` and `s` are the panel's allele frequencies and VanRaden
    /// scale (as the dense path already computes them).
    pub fn from_getter(
        n: usize,
        m: usize,
        get: impl Fn(usize, usize) -> u8,
        p: Vec<f64>,
        s: f64,
    ) -> Self {
        let bytes_per_col = n.div_ceil(4);
        let mut data = vec![0u8; bytes_per_col * m];
        let mut r = vec![0.0f64; n];
        for (j, &pj) in p.iter().enumerate() {
            let base = j * bytes_per_col;
            for (i, ri) in r.iter_mut().enumerate() {
                let g = get(i, j) & 0b11;
                data[base + (i >> 2)] |= g << ((i & 3) * 2);
                *ri += pj * g as f64;
            }
        }
        let p_sq = p.iter().map(|x| x * x).sum();
        Self {
            data,
            n,
            m,
            bytes_per_col,
            p,
            s,
            r,
            p_sq,
        }
    }

    pub fn n(&self) -> usize {
        self.n
    }
    pub fn m(&self) -> usize {
        self.m
    }
    /// Bytes held by the packed genotypes (the object whose size the memory claim is about).
    pub fn packed_bytes(&self) -> usize {
        self.data.len()
    }

    /// Raw dosage `M[i,j] ∈ {0,1,2}` as `f64`.
    #[inline]
    fn get(&self, i: usize, j: usize) -> f64 {
        let byte = self.data[j * self.bytes_per_col + (i >> 2)];
        ((byte >> ((i & 3) * 2)) & 0b11) as f64
    }

    /// `G·c = ridge·c + Z(Zᵀc)/s`, matrix-free from the packed genotypes. `O(n·m)`.
    ///
    /// Both marker passes are split across `std::thread`s over disjoint column ranges
    /// (no external runtime, no `unsafe`): phase 1 fills disjoint slices of `t`, phase 2
    /// accumulates a per-chunk partial `out` that is summed after the scope joins. On a
    /// serial machine (`available_parallelism` = 1) or a small panel it runs single-threaded.
    pub fn g_matvec(&self, c: &[f64], ridge: f64) -> Vec<f64> {
        debug_assert_eq!(c.len(), self.n);
        let sum_c: f64 = c.iter().sum();
        let hw = std::thread::available_parallelism()
            .map(|x| x.get())
            .unwrap_or(1);
        let nthreads = hw.min((self.m / 512).max(1));
        let chunk = self.m.div_ceil(nthreads);

        // Phase 1: t[j] = (Mᵀc)[j] − 2p[j]·(1ᵀc), over disjoint column chunks.
        let mut t = vec![0.0f64; self.m];
        std::thread::scope(|scope| {
            for (k, t_chunk) in t.chunks_mut(chunk).enumerate() {
                let (this, c) = (&*self, &*c);
                let j0 = k * chunk;
                scope.spawn(move || {
                    for (dj, tj) in t_chunk.iter_mut().enumerate() {
                        let j = j0 + dj;
                        let base = j * this.bytes_per_col;
                        let mut acc = 0.0;
                        for (i, &ci) in c.iter().enumerate() {
                            let g = (this.data[base + (i >> 2)] >> ((i & 3) * 2)) & 0b11;
                            acc += g as f64 * ci;
                        }
                        *tj = acc - 2.0 * this.p[j] * sum_c;
                    }
                });
            }
        });

        // Phase 2: out[i] = Σ_j M[i,j]·t[j] = (Z·t)[i] + 2p·(pᵀt) correction; each thread
        // builds a partial over its column chunk, summed after the scope.
        let p_dot_t: f64 = self.p.iter().zip(&t).map(|(&pj, &tj)| pj * tj).sum();
        let nchunks = self.m.div_ceil(chunk);
        let mut partials = vec![vec![0.0f64; self.n]; nchunks];
        std::thread::scope(|scope| {
            for (k, part) in partials.iter_mut().enumerate() {
                let (this, t) = (&*self, &t);
                let j0 = k * chunk;
                scope.spawn(move || {
                    let j1 = (j0 + chunk).min(this.m);
                    for (dj, &tj) in t[j0..j1].iter().enumerate() {
                        let base = (j0 + dj) * this.bytes_per_col;
                        for (i, o) in part.iter_mut().enumerate() {
                            let g = (this.data[base + (i >> 2)] >> ((i & 3) * 2)) & 0b11;
                            *o += g as f64 * tj;
                        }
                    }
                });
            }
        });
        let mut out = vec![0.0f64; self.n];
        for part in &partials {
            for (o, &p) in out.iter_mut().zip(part) {
                *o += p;
            }
        }
        for (o, &ci) in out.iter_mut().zip(c) {
            *o = ridge * ci + (*o - 2.0 * p_dot_t) / self.s;
        }
        out
    }

    /// Restricted Gram entry `G[i,j] = (zᵢ·zⱼ)/s + ridge·[i=j]`, matrix-free.
    pub fn gram_entry(&self, i: usize, j: usize, ridge: f64) -> f64 {
        let raw: f64 = (0..self.m).map(|l| self.get(i, l) * self.get(j, l)).sum();
        let z_ij = raw - 2.0 * self.r[i] - 2.0 * self.r[j] + 4.0 * self.p_sq;
        z_ij / self.s + if i == j { ridge } else { 0.0 }
    }
}

/// The two kinship operations the active-set solver needs from a genotype
/// representation: the full product `G·c` (pricing) and a restricted Gram entry
/// `G[i,j]` (the closed-form solve). Implementing it for both a dense `f64` `Z` and
/// the 2-bit [`PackedGeno`] lets one solver run on either, so the packed path reuses
/// the exact algorithm — only the representation changes.
pub trait Kinship {
    /// Number of candidates `n`.
    fn dim(&self) -> usize;
    /// `G·c = ridge·c + Z(Zᵀc)/s`.
    fn matvec(&self, c: &[f64], ridge: f64) -> Vec<f64>;
    /// `G[i,j] = (zᵢ·zⱼ)/s + ridge·[i=j]`.
    fn gram(&self, i: usize, j: usize, ridge: f64) -> f64;
}

impl Kinship for PackedGeno {
    fn dim(&self) -> usize {
        self.n
    }
    fn matvec(&self, c: &[f64], ridge: f64) -> Vec<f64> {
        self.g_matvec(c, ridge)
    }
    fn gram(&self, i: usize, j: usize, ridge: f64) -> f64 {
        self.gram_entry(i, j, ridge)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rand_distr::{Distribution, Uniform};

    /// Build raw M ∈ {0,1,2}, p, s and the centred f64 Z, plus a PackedGeno on the same M.
    fn fixture(n: usize, m: usize, seed: u64) -> (Vec<Vec<u8>>, Vec<f64>, f64, PackedGeno) {
        let mut rng = StdRng::seed_from_u64(seed);
        let gd = Uniform::new(0u8, 3u8).expect("0<3");
        let raw: Vec<Vec<u8>> = (0..n)
            .map(|_| (0..m).map(|_| gd.sample(&mut rng)).collect())
            .collect();
        let p: Vec<f64> = (0..m)
            .map(|j| (0..n).map(|i| raw[i][j] as f64).sum::<f64>() / (2.0 * n as f64))
            .collect();
        let s = 2.0 * p.iter().map(|pj| pj * (1.0 - pj)).sum::<f64>();
        let packed = PackedGeno::from_getter(n, m, |i, j| raw[i][j], p.clone(), s);
        (raw, p, s, packed)
    }

    /// Dense f64 reference: Gc = ridge·c + Z(Zᵀc)/s with Z = M − 2p.
    fn gc_ref(raw: &[Vec<u8>], p: &[f64], s: f64, c: &[f64], ridge: f64) -> Vec<f64> {
        let n = raw.len();
        let m = p.len();
        let z = |i: usize, j: usize| raw[i][j] as f64 - 2.0 * p[j];
        let t: Vec<f64> = (0..m)
            .map(|j| (0..n).map(|i| z(i, j) * c[i]).sum())
            .collect();
        (0..n)
            .map(|i| ridge * c[i] + (0..m).map(|j| z(i, j) * t[j]).sum::<f64>() / s)
            .collect()
    }

    #[test]
    fn matvec_matches_f64() {
        let (raw, p, s, packed) = fixture(200, 500, 7);
        let mut rng = StdRng::seed_from_u64(99);
        let cd = Uniform::new(-1.0f64, 1.0).expect("valid");
        let c: Vec<f64> = (0..200).map(|_| cd.sample(&mut rng)).collect();
        let ridge = 1e-5;
        let got = packed.g_matvec(&c, ridge);
        let want = gc_ref(&raw, &p, s, &c, ridge);
        let max_rel = (0..200)
            .map(|i| (got[i] - want[i]).abs())
            .fold(0.0, f64::max)
            / want.iter().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(
            max_rel < 1e-12,
            "packed matvec diverges: max_rel={max_rel:e}"
        );
    }

    #[test]
    fn gram_entry_matches_f64() {
        let (raw, p, s, packed) = fixture(120, 400, 3);
        let ridge = 1e-5;
        let z = |i: usize, j: usize| raw[i][j] as f64 - 2.0 * p[j];
        let mut max_rel = 0.0f64;
        for i in [0usize, 5, 60, 119] {
            for j in [1usize, 5, 60, 119] {
                let want: f64 = (0..400).map(|l| z(i, l) * z(j, l)).sum::<f64>() / s
                    + if i == j { ridge } else { 0.0 };
                let got = packed.gram_entry(i, j, ridge);
                max_rel = max_rel.max((got - want).abs() / want.abs().max(1e-12));
            }
        }
        assert!(max_rel < 1e-12, "packed Gram diverges: max_rel={max_rel:e}");
    }

    #[test]
    fn memory_is_two_bits_per_genotype() {
        let (_raw, _p, _s, packed) = fixture(1000, 800, 1);
        // 4 genotypes per byte, column-major: m · ceil(n/4).
        assert_eq!(packed.packed_bytes(), 800 * 250);
        // 32× smaller than the f64 store (8 bytes each).
        assert_eq!(1000 * 800 * 8 / packed.packed_bytes(), 32);
    }
}
