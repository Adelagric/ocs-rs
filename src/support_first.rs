//! Support-first exact OCS solver (Rust port of the `research/` prototype).
//!
//! The OCS optimum activates only a tiny support `S`. For a fixed support the
//! problem `max b_Sᵀc_S s.t. 1ᵀc_S=1, c_Sᵀ G_SS c_S = k` is "maximise a linear
//! form over an ellipsoid", whose multiplier `μ` solves a scalar quadratic — so
//! each restricted solve is two Cholesky back-substitutions on the small `G_SS`,
//! no iteration. The whole cost is identifying `S`, done by column generation:
//! add the best reduced-cost candidate, or — if the support cannot yet satisfy
//! the kinship bound — the least related one. Every full `G·c` is formed
//! matrix-free as `ridge·c + Z(Zᵀc)/s`; `G` is never materialised.
//!
//! `G` here means the *ridged* `G+εI` that Clarabel's cone actually enforces, so
//! the two solvers are compared on the same constraint.

use faer::linalg::solvers::Solve;
use faer::{Mat, Side};

/// Terminal status of the active-set loop.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SfStatus {
    /// KKT-optimal: feasible and no candidate has a positive reduced cost.
    Solved,
    /// Iteration cap hit (should not happen for well-posed instances).
    MaxIter,
}

/// What the solve produced.
///
/// `c`/`gain` are a usable optimum **only when `status == Solved`**. On
/// `MaxIter` (degenerate inputs the active set could not resolve) the iterate may
/// be infeasible — check `quad <= k` before trusting it.
#[derive(Clone, Debug)]
pub struct SupportFirstOutcome {
    /// Optimal contributions (length `n`); valid only when `status == Solved`.
    pub c: Vec<f64>,
    /// Active support (indices with `cᵢ > 0`), sorted.
    pub support: Vec<usize>,
    /// Outer iterations.
    pub iterations: u32,
    /// Number of full `G·c` products (the dominant `O(n·m)` cost).
    pub products: u32,
    /// Genetic gain `bᵀc`.
    pub gain: f64,
    /// Kinship `cᵀ(G+εI)c` (compared against `k`).
    pub quad: f64,
    /// Terminal status.
    pub status: SfStatus,
}

/// `G·c = ridge·c + Z(Zᵀc)/s`, never forming `G`. Cost `O(n·m)`.
fn g_matvec(z: &Mat<f64>, s: f64, ridge: f64, c: &[f64]) -> Vec<f64> {
    let n = z.nrows();
    let cm = Mat::from_fn(n, 1, |i, _| c[i]);
    let t = z.transpose() * cm.as_ref(); // m×1  = Zᵀc
    let u = z * t.as_ref(); // n×1  = Z(Zᵀc)
    let inv_s = 1.0 / s;
    (0..n).map(|i| u[(i, 0)] * inv_s + ridge * c[i]).collect()
}

/// `G_SS = Z_S Z_Sᵀ / s + ridge·I` (small, `|S|×|S|`).
fn build_gss(z: &Mat<f64>, s: f64, ridge: f64, support: &[usize]) -> Mat<f64> {
    let ns = support.len();
    let m = z.ncols();
    let inv_s = 1.0 / s;
    Mat::from_fn(ns, ns, |i, j| {
        let (si, sj) = (support[i], support[j]);
        let mut acc = 0.0;
        for l in 0..m {
            acc += z[(si, l)] * z[(sj, l)];
        }
        acc * inv_s + if i == j { ridge } else { 0.0 }
    })
}

/// Real roots of `a x² + b x + c = 0` (handles the near-linear case). The
/// degeneracy test is relative to the coefficient scale so it does not depend on
/// the magnitude of `G`.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    let scale = (a.abs() + b.abs() + c.abs()).max(1.0);
    if a.abs() < 1e-14 * scale {
        if b.abs() < 1e-14 * scale {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    vec![(-b + sq) / (2.0 * a), (-b - sq) / (2.0 * a)]
}

/// Closed-form restricted solve on a support. Returns `(c_S, μ, λ)`, or `None`
/// if the kinship ellipsoid does not meet the affine hull on `S` (infeasible).
fn closed_form(g_ss: &Mat<f64>, b_s: &[f64], k: f64) -> Option<(Vec<f64>, f64, f64)> {
    let ns = b_s.len();
    if ns == 1 {
        // c = [1]; feasible iff the single self-relationship is within k.
        return if g_ss[(0, 0)] <= k {
            Some((vec![1.0], b_s[0], 0.0))
        } else {
            None
        };
    }
    let llt = g_ss.llt(Side::Lower).ok()?;
    // Solve G_SS [w v] = [1 b_S].
    let rhs = Mat::from_fn(ns, 2, |i, j| if j == 0 { 1.0 } else { b_s[i] });
    let sol = llt.solve(rhs.as_ref());
    let w: Vec<f64> = (0..ns).map(|i| sol[(i, 0)]).collect();
    let v: Vec<f64> = (0..ns).map(|i| sol[(i, 1)]).collect();
    let alpha: f64 = w.iter().sum(); // 1ᵀw = 1ᵀG⁻¹1
    let beta: f64 = v.iter().sum(); // 1ᵀv = 1ᵀG⁻¹b
    let delta: f64 = b_s.iter().zip(&v).map(|(b, vi)| b * vi).sum(); // bᵀG⁻¹b

    // (kα−1)(α μ² − 2β μ) + kβ² − δ = 0
    let qa = (k * alpha - 1.0) * alpha;
    let qb = -2.0 * beta * (k * alpha - 1.0);
    let qc = k * beta * beta - delta;

    let mut best: Option<(Vec<f64>, f64, f64, f64)> = None; // (c_S, μ, λ, gain)
    for mu in solve_quadratic(qa, qb, qc) {
        let two_lam = beta - mu * alpha;
        if two_lam <= 0.0 {
            continue; // need λ > 0 (binding constraint, correct multiplier sign)
        }
        let cs: Vec<f64> = (0..ns).map(|i| (v[i] - mu * w[i]) / two_lam).collect();
        let gain: f64 = b_s.iter().zip(&cs).map(|(b, c)| b * c).sum();
        if best.as_ref().is_none_or(|b| gain > b.3) {
            best = Some((cs, mu, two_lam / 2.0, gain));
        }
    }
    best.map(|(cs, mu, lam, _)| (cs, mu, lam))
}

fn argmax(v: &[f64]) -> usize {
    v.iter()
        .enumerate()
        .fold(0, |best, (i, &x)| if x > v[best] { i } else { best })
}

/// Solve OCS by support-first column generation.
///
/// `z` centred genotypes (`n×m`), `s` VanRaden scale, `ridge` ε, `b` GEBV, `k`
/// kinship bound. The kinship enforced is `cᵀ(G+εI)c ≤ k` (same as Clarabel).
pub fn solve(
    z: &Mat<f64>,
    s: f64,
    ridge: f64,
    b: &[f64],
    k: f64,
    max_iter: u32,
    tol: f64,
) -> SupportFirstOutcome {
    let n = b.len();
    let mut support = vec![argmax(b)];
    let mut c_cur = vec![0.0; n]; // last feasible iterate (basis for diversification)
    c_cur[support[0]] = 1.0;
    let mut products = 0u32;
    // Anti-cycling: indices dropped for negativity are taboo for re-addition
    // until the next genuine progress step (a reduced-cost add), which clears
    // the set. Without this, a drop→infeasible→re-add stall could loop to
    // `max_iter`. In the VanRaden regime no drop recurs, so this is inert.
    let mut dropped = vec![false; n];

    for it in 0..max_iter {
        let b_s: Vec<f64> = support.iter().map(|&i| b[i]).collect();
        let g_ss = build_gss(z, s, ridge, &support);

        match closed_form(&g_ss, &b_s, k) {
            None => {
                // Support cannot satisfy k yet: add the least related candidate.
                let gc = g_matvec(z, s, ridge, &c_cur);
                products += 1;
                let mut best_j = None;
                let mut best_val = f64::INFINITY;
                for (j, &gj) in gc.iter().enumerate() {
                    if !support.contains(&j) && !dropped[j] && gj < best_val {
                        best_val = gj;
                        best_j = Some(j);
                    }
                }
                match best_j {
                    Some(j) => support.push(j),
                    None => break, // no fresh candidate can diversify: bail to MaxIter
                }
            }
            Some((cs, mu, lam)) => {
                if cs.iter().any(|&x| x < -tol) {
                    // Drop negative contributions, re-solve on the smaller
                    // support; the dropped indices become taboo (anti-cycling).
                    let mut keep = Vec::with_capacity(support.len());
                    for idx in 0..support.len() {
                        if cs[idx] > tol {
                            keep.push(support[idx]);
                        } else {
                            dropped[support[idx]] = true;
                        }
                    }
                    support = keep;
                    if support.is_empty() {
                        support.push(argmax(b));
                    }
                    continue;
                }
                let mut c = vec![0.0; n];
                for idx in 0..support.len() {
                    c[support[idx]] = cs[idx];
                }
                c_cur = c;

                // Reduced costs rⱼ = bⱼ − μ − 2λ (G c)ⱼ ; add the best violator.
                let gc = g_matvec(z, s, ridge, &c_cur);
                products += 1;
                let mut best_j = None;
                let mut best_r = tol;
                for (j, &gj) in gc.iter().enumerate() {
                    if support.contains(&j) {
                        continue;
                    }
                    let r = b[j] - mu - 2.0 * lam * gj;
                    if r > best_r {
                        best_r = r;
                        best_j = Some(j);
                    }
                }
                match best_j {
                    None => {
                        let mut sorted = support.clone();
                        sorted.sort_unstable();
                        let gain = b.iter().zip(&c_cur).map(|(b, c)| b * c).sum();
                        let quad = c_cur.iter().zip(&gc).map(|(c, g)| c * g).sum();
                        return SupportFirstOutcome {
                            c: c_cur,
                            support: sorted,
                            iterations: it + 1,
                            products,
                            gain,
                            quad,
                            status: SfStatus::Solved,
                        };
                    }
                    Some(j) => {
                        // Genuine progress (gain strictly improves): clear taboos.
                        dropped.iter_mut().for_each(|d| *d = false);
                        support.push(j);
                    }
                }
            }
        }
    }

    // Cap hit: return the best feasible iterate found.
    let gc = g_matvec(z, s, ridge, &c_cur);
    let mut sorted: Vec<usize> = (0..n).filter(|&i| c_cur[i] > tol).collect();
    sorted.sort_unstable();
    SupportFirstOutcome {
        gain: b.iter().zip(&c_cur).map(|(b, c)| b * c).sum(),
        quad: c_cur.iter().zip(&gc).map(|(c, g)| c * g).sum(),
        c: c_cur,
        support: sorted,
        iterations: max_iter,
        products: products + 1,
        status: SfStatus::MaxIter,
    }
}

/// Index of the highest-`b` candidate whose male flag equals `want`, if any.
fn argmax_masked(b: &[f64], male: &[bool], want: bool) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, &x) in b.iter().enumerate() {
        if male[i] == want && best.is_none_or(|j| x > b[j]) {
            best = Some(i);
        }
    }
    best
}

/// Closed-form restricted solve for the **sexed** OCS, where contributions split
/// half to males and half to females: `max b_Sᵀc_S s.t. A_S c_S = [½,½],
/// c_Sᵀ G_SS c_S = k`, with `A_S` the 2×|S| sex-incidence (row 0 males, row 1
/// females). Returns `(c_S, μ_male, μ_female, λ)`, or `None` if the support lacks
/// a sex or the kinship ellipsoid does not meet the affine hull on `S`.
///
/// For `|S| > 2` the two equalities leave a non-trivial null space. Rather than
/// form it (the prototype used an SVD), the equalities are eliminated through
/// `P = A_S G⁻¹ A_Sᵀ` (2×2): with one Cholesky solve `G_SS [W | u] = [A_Sᵀ | b_S]`
/// the optimum is `c_S = g/(2λ) + h` where `g = u − W P⁻¹q`, `h = W P⁻¹d`, and the
/// active ellipsoid reduces to the scalar quadratic `4(C−k)λ² + 4Bλ + A = 0`.
/// For the fully-determined `|S| = 2` (one male, one female) the contributions
/// are forced to `[½,½]`; `λ = 0` is accepted there (the binding case is
/// measure-zero, mirroring the single-support branch of [`closed_form`]).
fn closed_form_sexed(
    g_ss: &Mat<f64>,
    b_s: &[f64],
    male_s: &[bool],
    k: f64,
) -> Option<(Vec<f64>, f64, f64, f64)> {
    let ns = b_s.len();
    let n_male = male_s.iter().filter(|&&m| m).count();
    if n_male == 0 || n_male == ns {
        return None; // both sexes must be present to meet Σ = ½ on each row
    }

    if ns == 2 {
        // A_S is 2×2 of full rank ⇒ c_S is forced to ½ each.
        let cs = vec![0.5_f64, 0.5_f64];
        let qv = 0.25 * (g_ss[(0, 0)] + g_ss[(1, 1)] + 2.0 * g_ss[(0, 1)]);
        if qv > k + 1e-9 {
            return None;
        }
        let (mu_m, mu_f) = if male_s[0] {
            (b_s[0], b_s[1])
        } else {
            (b_s[1], b_s[0])
        };
        return Some((cs, mu_m, mu_f, 0.0));
    }

    let llt = g_ss.llt(Side::Lower).ok()?;
    // Solve G_SS [W | u] = [A_Sᵀ | b_S]: W = G⁻¹A_Sᵀ (cols 0,1), u = G⁻¹b_S (col 2).
    let rhs = Mat::from_fn(ns, 3, |i, j| match j {
        0 => f64::from(male_s[i]),
        1 => f64::from(!male_s[i]),
        _ => b_s[i],
    });
    let sol = llt.solve(rhs.as_ref());

    // P = A_S W (2×2) and q = A_S u (2): sum the W/u rows by sex.
    let (mut p00, mut p01, mut p10, mut p11, mut q0, mut q1) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..ns {
        let (w0, w1, ui) = (sol[(i, 0)], sol[(i, 1)], sol[(i, 2)]);
        if male_s[i] {
            p00 += w0;
            p01 += w1;
            q0 += ui;
        } else {
            p10 += w0;
            p11 += w1;
            q1 += ui;
        }
    }
    let det = p00 * p11 - p01 * p10;
    let scale = (p00.abs() + p01.abs() + p10.abs() + p11.abs()).max(1.0);
    if det.abs() < 1e-14 * scale {
        return None;
    }
    // e = P⁻¹q, f = P⁻¹d (d = [½,½]), via the explicit 2×2 inverse.
    let apply_pinv = |x0: f64, x1: f64| ((p11 * x0 - p01 * x1) / det, (p00 * x1 - p10 * x0) / det);
    let (e0, e1) = apply_pinv(q0, q1);
    let (f0, f1) = apply_pinv(0.5, 0.5);

    // g = u − W e, h = W f; their images G g = b_S − A_Sᵀe and G h = A_Sᵀf give
    // A = gᵀGg, B = gᵀGh, C = hᵀGh without a second solve.
    let mut g = vec![0.0; ns];
    let mut h = vec![0.0; ns];
    let (mut a_coef, mut b_coef, mut c_coef) = (0.0, 0.0, 0.0);
    for i in 0..ns {
        let (w0, w1, ui) = (sol[(i, 0)], sol[(i, 1)], sol[(i, 2)]);
        let (ei, fi) = if male_s[i] { (e0, f0) } else { (e1, f1) };
        g[i] = ui - (w0 * e0 + w1 * e1);
        h[i] = w0 * f0 + w1 * f1;
        let gg = b_s[i] - ei; // (G g)_i
        let gh = fi; // (G h)_i
        a_coef += g[i] * gg;
        b_coef += g[i] * gh;
        c_coef += h[i] * gh;
    }

    let mut best: Option<(Vec<f64>, f64, f64, f64, f64)> = None; // (c,μm,μf,λ,gain)
    for lam in solve_quadratic(4.0 * (c_coef - k), 4.0 * b_coef, a_coef) {
        if lam <= 0.0 {
            continue;
        }
        let cs: Vec<f64> = (0..ns).map(|i| g[i] / (2.0 * lam) + h[i]).collect();
        let gain: f64 = b_s.iter().zip(&cs).map(|(b, c)| b * c).sum();
        if best.as_ref().is_none_or(|bb| gain > bb.4) {
            best = Some((cs, e0 - 2.0 * lam * f0, e1 - 2.0 * lam * f1, lam, gain));
        }
    }
    best.map(|(cs, mm, mf, lam, _)| (cs, mm, mf, lam))
}

/// Solve the **sexed** OCS by support-first column generation.
///
/// Identical machinery to [`solve`], but the simplex `1ᵀc = 1` is replaced by the
/// two sex equalities `Σ_{males} c = Σ_{females} c = ½` (the *true* OCS, mates
/// drawn one of each sex). `male[i]` is the sex flag. Seeds the support with the
/// best male and best female; reduced costs price candidates against the
/// sex-specific multiplier `μ_{sex(j)}`.
#[allow(clippy::too_many_arguments)]
pub fn solve_sexed(
    z: &Mat<f64>,
    s: f64,
    ridge: f64,
    b: &[f64],
    male: &[bool],
    k: f64,
    max_iter: u32,
    tol: f64,
) -> SupportFirstOutcome {
    let n = b.len();
    let (best_m, best_f) = match (argmax_masked(b, male, true), argmax_masked(b, male, false)) {
        (Some(m), Some(f)) => (m, f),
        _ => {
            // A whole sex is missing: the sexed problem is infeasible.
            return SupportFirstOutcome {
                c: vec![0.0; n],
                support: Vec::new(),
                iterations: 0,
                products: 0,
                gain: 0.0,
                quad: 0.0,
                status: SfStatus::MaxIter,
            };
        }
    };
    let mut support = vec![best_m, best_f];
    let mut c_cur = vec![0.0; n];
    c_cur[best_m] = 0.5;
    c_cur[best_f] = 0.5;
    let mut products = 0u32;
    let mut dropped = vec![false; n];

    for it in 0..max_iter {
        let b_s: Vec<f64> = support.iter().map(|&i| b[i]).collect();
        let male_s: Vec<bool> = support.iter().map(|&i| male[i]).collect();
        let g_ss = build_gss(z, s, ridge, &support);

        match closed_form_sexed(&g_ss, &b_s, &male_s, k) {
            None => {
                let gc = g_matvec(z, s, ridge, &c_cur);
                products += 1;
                let mut best_j = None;
                let mut best_val = f64::INFINITY;
                for (j, &gj) in gc.iter().enumerate() {
                    if !support.contains(&j) && !dropped[j] && gj < best_val {
                        best_val = gj;
                        best_j = Some(j);
                    }
                }
                match best_j {
                    Some(j) => support.push(j),
                    None => break,
                }
            }
            Some((cs, mu_m, mu_f, lam)) => {
                if cs.iter().any(|&x| x < -tol) {
                    let mut keep = Vec::with_capacity(support.len());
                    for idx in 0..support.len() {
                        if cs[idx] > tol {
                            keep.push(support[idx]);
                        } else {
                            dropped[support[idx]] = true;
                        }
                    }
                    support = keep;
                    // A sex must stay represented (Σ = ½ on each row): re-seed it.
                    if !support.iter().any(|&i| male[i]) {
                        support.push(best_m);
                    }
                    if !support.iter().any(|&i| !male[i]) {
                        support.push(best_f);
                    }
                    continue;
                }
                let mut c = vec![0.0; n];
                for idx in 0..support.len() {
                    c[support[idx]] = cs[idx];
                }
                c_cur = c;

                let gc = g_matvec(z, s, ridge, &c_cur);
                products += 1;
                let mut best_j = None;
                let mut best_r = tol;
                for (j, &gj) in gc.iter().enumerate() {
                    if support.contains(&j) {
                        continue;
                    }
                    let mu = if male[j] { mu_m } else { mu_f };
                    let r = b[j] - mu - 2.0 * lam * gj;
                    if r > best_r {
                        best_r = r;
                        best_j = Some(j);
                    }
                }
                match best_j {
                    None => {
                        let mut sorted = support.clone();
                        sorted.sort_unstable();
                        let gain = b.iter().zip(&c_cur).map(|(b, c)| b * c).sum();
                        let quad = c_cur.iter().zip(&gc).map(|(c, g)| c * g).sum();
                        return SupportFirstOutcome {
                            c: c_cur,
                            support: sorted,
                            iterations: it + 1,
                            products,
                            gain,
                            quad,
                            status: SfStatus::Solved,
                        };
                    }
                    Some(j) => {
                        dropped.iter_mut().for_each(|d| *d = false);
                        support.push(j);
                    }
                }
            }
        }
    }

    let gc = g_matvec(z, s, ridge, &c_cur);
    let mut sorted: Vec<usize> = (0..n).filter(|&i| c_cur[i] > tol).collect();
    sorted.sort_unstable();
    SupportFirstOutcome {
        gain: b.iter().zip(&c_cur).map(|(b, c)| b * c).sum(),
        quad: c_cur.iter().zip(&gc).map(|(c, g)| c * g).sum(),
        c: c_cur,
        support: sorted,
        iterations: max_iter,
        products: products + 1,
        status: SfStatus::MaxIter,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::datagen;
    use crate::grm;
    use crate::socp::{self, Factor};
    use crate::solve as conic;

    #[test]
    fn closed_form_two_variable() {
        // G_SS = I, b = [1,2], k = 0.6 → c = [0.276393, 0.723607] (hand-checked).
        let g = Mat::<f64>::identity(2, 2);
        let (cs, _mu, lam) = closed_form(&g, &[1.0, 2.0], 0.6).unwrap();
        assert!(lam > 0.0);
        assert!((cs[0] - 0.276393).abs() < 1e-5);
        assert!((cs[1] - 0.723607).abs() < 1e-5);
        assert!((cs[0] + cs[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn infeasible_single_support_is_none() {
        // One candidate with self-relationship above k cannot satisfy it.
        let g = Mat::<f64>::from_fn(1, 1, |_, _| 1.0);
        assert!(closed_form(&g, &[1.0], 0.6).is_none());
    }

    #[test]
    fn matches_clarabel_small() {
        // The crux test: support-first must reach the same optimum as the conic
        // IPM on the same data.
        let d = datagen::generate(60, 2000, 20240617);
        let ridge = 1e-5;
        let grm = grm::Grm::build(&d.z, d.s, ridge);
        let l = grm.cholesky_lower().unwrap();
        let mean_diag: f64 = (0..d.n).map(|i| grm.g[(i, i)]).sum::<f64>() / d.n as f64;
        let k = 0.6 * mean_diag;

        let prob = socp::build(Factor::Cholesky(&l), &d.b, k, d.s, None);
        let clar = conic::solve(&prob, conic::SolveConfig::default()).unwrap();

        let sf = solve(&d.z, d.s, ridge, &d.b, k, 500, 1e-7);

        assert_eq!(sf.status, SfStatus::Solved);
        assert!(
            (sf.gain - clar.gain).abs() < 1e-6,
            "support-first gain {} vs clarabel {}",
            sf.gain,
            clar.gain
        );
        assert!(sf.quad <= k + 1e-6, "kinship {} > k {}", sf.quad, k);
        assert!(sf.c.iter().all(|&c| c >= -1e-7));
        assert!((sf.c.iter().sum::<f64>() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn solved_status_implies_feasible_across_k() {
        // Safety invariant (guards the active-set anti-cycling fix): across a
        // range of k including very tight values, the solver must terminate and
        // a `Solved` result must be feasible and on the simplex. `MaxIter` is
        // allowed (degenerate), but must never be reported as `Solved`.
        let d = datagen::generate(50, 1000, 11);
        let ridge = 1e-5;
        let grm = grm::Grm::build(&d.z, d.s, ridge);
        let mean_diag: f64 = (0..d.n).map(|i| grm.g[(i, i)]).sum::<f64>() / d.n as f64;
        for frac in [0.03, 0.1, 0.3, 0.6, 0.9] {
            let k = frac * mean_diag;
            let sf = solve(&d.z, d.s, ridge, &d.b, k, 2000, 1e-7);
            if sf.status == SfStatus::Solved {
                assert!(
                    sf.quad <= k + 1e-6,
                    "Solved but infeasible at frac={frac}: quad {} > k {}",
                    sf.quad,
                    k
                );
                assert!((sf.c.iter().sum::<f64>() - 1.0).abs() < 1e-6);
                assert!(sf.c.iter().all(|&c| c >= -1e-7));
            }
        }
    }

    #[test]
    fn sexed_two_support_forced_half() {
        // One male + one female ⇒ contributions forced to [½, ½]; λ = 0.
        let g = Mat::<f64>::identity(2, 2);
        let (cs, mu_m, mu_f, lam) =
            closed_form_sexed(&g, &[1.0, 2.0], &[true, false], 10.0).unwrap();
        assert!((cs[0] - 0.5).abs() < 1e-12 && (cs[1] - 0.5).abs() < 1e-12);
        assert_eq!(lam, 0.0);
        assert!((mu_m - 1.0).abs() < 1e-12 && (mu_f - 2.0).abs() < 1e-12);
    }

    #[test]
    fn sexed_single_sex_is_none() {
        // No female ⇒ Σ_females = ½ is unsatisfiable.
        let g = Mat::<f64>::identity(2, 2);
        assert!(closed_form_sexed(&g, &[1.0, 2.0], &[true, true], 10.0).is_none());
    }

    #[test]
    fn sexed_three_support_kkt() {
        // G_SS = I, sexes [M,F,M], k binding. The closed form must satisfy the two
        // sex sums, the active ellipsoid, and stationarity bᵢ − μ_sex − 2λcᵢ = 0.
        let g = Mat::<f64>::identity(3, 3);
        let b = [3.0, 1.0, 2.0];
        let male = [true, false, true];
        let k = 0.4;
        let (cs, mu_m, mu_f, lam) = closed_form_sexed(&g, &b, &male, k).unwrap();
        assert!(((cs[0] + cs[2]) - 0.5).abs() < 1e-9, "male sum");
        assert!((cs[1] - 0.5).abs() < 1e-9, "female sum");
        let quad: f64 = cs.iter().map(|c| c * c).sum();
        assert!(
            (quad - k).abs() < 1e-9,
            "ellipsoid not active: {quad} vs {k}"
        );
        assert!(lam > 0.0);
        for i in 0..3 {
            let mu = if male[i] { mu_m } else { mu_f };
            let res = b[i] - mu - 2.0 * lam * cs[i];
            assert!(res.abs() < 1e-7, "stationarity i={i}: {res}");
        }
    }

    #[test]
    fn sexed_solve_kkt_certificate() {
        // Full solve on synthetic data: a Solved result must be feasible, on the
        // simplex, split ½/½ by sex, and non-negative.
        let d = datagen::generate(60, 2000, 7);
        let ridge = 1e-5;
        let grm = grm::Grm::build(&d.z, d.s, ridge);
        let mean_diag: f64 = (0..d.n).map(|i| grm.g[(i, i)]).sum::<f64>() / d.n as f64;
        let k = 0.5 * mean_diag;
        let male: Vec<bool> = (0..d.n).map(|i| i % 2 == 0).collect();

        let sf = solve_sexed(&d.z, d.s, ridge, &d.b, &male, k, 1000, 1e-7);
        assert_eq!(sf.status, SfStatus::Solved);
        assert!(sf.quad <= k + 1e-6, "kinship {} > k {}", sf.quad, k);
        assert!(sf.c.iter().all(|&c| c >= -1e-7));
        let sum_m: f64 = (0..d.n).filter(|&i| male[i]).map(|i| sf.c[i]).sum();
        let sum_f: f64 = (0..d.n).filter(|&i| !male[i]).map(|i| sf.c[i]).sum();
        assert!((sum_m - 0.5).abs() < 1e-6, "male sum {sum_m}");
        assert!((sum_f - 0.5).abs() < 1e-6, "female sum {sum_f}");
    }

    #[test]
    fn sexed_solved_feasible_across_k() {
        // Anti-cycling / termination invariant for the sexed loop: across tight to
        // loose k, a `Solved` result is always feasible and correctly sex-split.
        let d = datagen::generate(50, 1000, 11);
        let ridge = 1e-5;
        let grm = grm::Grm::build(&d.z, d.s, ridge);
        let mean_diag: f64 = (0..d.n).map(|i| grm.g[(i, i)]).sum::<f64>() / d.n as f64;
        let male: Vec<bool> = (0..d.n).map(|i| i % 2 == 0).collect();
        for frac in [0.05, 0.1, 0.3, 0.6, 0.9] {
            let k = frac * mean_diag;
            let sf = solve_sexed(&d.z, d.s, ridge, &d.b, &male, k, 2000, 1e-7);
            if sf.status == SfStatus::Solved {
                assert!(sf.quad <= k + 1e-6, "Solved infeasible at frac={frac}");
                let sum_m: f64 = (0..d.n).filter(|&i| male[i]).map(|i| sf.c[i]).sum();
                let sum_f: f64 = (0..d.n).filter(|&i| !male[i]).map(|i| sf.c[i]).sum();
                assert!((sum_m - 0.5).abs() < 1e-6 && (sum_f - 0.5).abs() < 1e-6);
                assert!(sf.c.iter().all(|&c| c >= -1e-7));
            }
        }
    }
}
