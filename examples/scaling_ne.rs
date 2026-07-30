//! Honest operating-point scaling.
//!
//! Annotate the current Figure-1A points with the breeder quantities they actually
//! sit at — the effective number of contributing parents Ne and the implied rate of
//! inbreeding ΔF — and re-measure the same n-sweep at a FIXED breeder-defensible
//! operating point (a target Ne) rather than at a fixed absolute cap.
//!
//!   Ne_c = 1 / Σcᵢ²   (effective number of parents; = K for K equal contributions
//!                       summing to one — the sexed budget has Σcᵢ = 1)
//!   ΔF   ≈ 1 / (2 Ne_c)
//!   f    = cᵀG c       (realised offspring coancestry, matrix-free ‖Zᵀc‖²/s)
//!
//! Caveat, kept in the header so it cannot be lost: the synthetic generator is HWE
//! with independent loci — no family structure, no LD, no founders. Growing n at
//! fixed m therefore *dilutes* relatedness, the most favourable possible case for a
//! small support. And fixing Ne bounds |S| almost by construction. So this file
//! measures the COST at a defensible operating point and quantifies where Figure 1A
//! actually sits; the "support bounded as n grows" claim needs a STRUCTURED
//! population (fixed base Ne) and is tested separately.
//!
//!   cargo run --release --example scaling_ne

use std::time::Instant;

use ocs_rs::datagen;
use ocs_rs::grm::quad_form_z;
use ocs_rs::support_first::{SfStatus, SupportFirstOutcome, solve_sexed};

fn sum_sq(c: &[f64]) -> f64 {
    c.iter().map(|x| x * x).sum()
}

/// Ne_c = 1 / Σcᵢ² (Σcᵢ = 1 for the sexed budget).
fn ne_c(c: &[f64]) -> f64 {
    1.0 / sum_sq(c)
}

fn support_by_sex(support: &[usize], male: &[bool]) -> (usize, usize) {
    let males = support.iter().filter(|&&i| male[i]).count();
    (males, support.len() - males)
}

fn main() {
    let m = 1000usize;
    let ridge = 1e-5;
    let ns = [1000usize, 2000, 5000, 10000, 20000, 40000];
    let targets_ne = [100.0_f64, 50.0, 25.0]; // ΔF = 0.5%, 1%, 2%

    let mut csv = String::from(
        "mode,n,m,k,realised_f,support,support_male,support_female,ne_c,delta_f_pct,solve_s,status\n",
    );

    println!("honest operating-point scaling   m={m}  (sexed OCS, release)");
    println!(
        "{:>14} {:>7} {:>10} {:>11} {:>6} {:>6} {:>8} {:>8} {:>9}",
        "mode", "n", "k", "realised_f", "|S|", "|S_m|", "Ne_c", "dF %", "solve s"
    );

    for &n in &ns {
        let d = datagen::generate(n, m, 20240618 + n as u64);
        let male: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let n_m = male.iter().filter(|&&x| x).count();
        let n_f = n - n_m;

        // Mean GRM diagonal without forming G: G_ii = (z_i·z_i)/s + ridge.
        let mean_diag: f64 = (0..n)
            .map(|i| (0..m).map(|l| d.z[(i, l)].powi(2)).sum::<f64>() / d.s + ridge)
            .sum::<f64>()
            / n as f64;

        // Feasible tight bracket for bisection: coancestry of the uniform
        // within-sex allocation (always feasible, high Ne).
        let mut c_unif = vec![0.0f64; n];
        for (i, cu) in c_unif.iter_mut().enumerate() {
            *cu = if male[i] {
                0.5 / n_m as f64
            } else {
                0.5 / n_f as f64
            };
        }
        let k_tight = quad_form_z(d.z.as_ref(), d.s, &c_unif);
        let k_loose = mean_diag; // loose: collapses toward the 2-support truncation

        // Warm-up (caches) at a mid cap.
        let _ = solve_sexed(
            d.z.as_ref(),
            d.s,
            ridge,
            &d.b,
            &male,
            0.5 * (k_tight + k_loose),
            50_000,
            1e-8,
        );

        let emit = |mode: &str, k: f64, out: &SupportFirstOutcome, secs: f64, csv: &mut String| {
            let f = quad_form_z(d.z.as_ref(), d.s, &out.c);
            let ne = ne_c(&out.c);
            let df = 50.0 / ne; // ΔF% = 100 / (2 Ne)
            let (sm, sf) = support_by_sex(&out.support, &male);
            println!(
                "{mode:>14} {n:>7} {k:>10.5} {f:>11.6} {:>6} {sm:>6} {ne:>8.1} {df:>8.3} {secs:>9.4}",
                out.support.len()
            );
            csv.push_str(&format!(
                "{mode},{n},{m},{k:.6},{f:.6},{},{sm},{sf},{ne:.2},{df:.3},{secs:.4},{:?}\n",
                out.support.len(),
                out.status
            ));
        };

        // Mode A — the current Figure-1A cap: fixed ABSOLUTE k = 0.1·mean_diag.
        {
            let k = 0.1 * mean_diag;
            let t = Instant::now();
            let out = solve_sexed(d.z.as_ref(), d.s, ridge, &d.b, &male, k, 50_000, 1e-8);
            let secs = t.elapsed().as_secs_f64();
            emit("fixed_abs_k", k, &out, secs, &mut csv);
        }

        // Mode B — fixed operating point: bisect k to hit each target Ne.
        for &t_ne in &targets_ne {
            let (mut lo, mut hi) = (k_tight, k_loose);
            let mut chosen_k = 0.5 * (lo + hi);
            for _ in 0..40 {
                let k = 0.5 * (lo + hi);
                let out = solve_sexed(d.z.as_ref(), d.s, ridge, &d.b, &male, k, 50_000, 1e-8);
                chosen_k = k;
                if out.status != SfStatus::Solved {
                    lo = k; // degenerate ⇒ treat as too tight, loosen
                    continue;
                }
                let ne = ne_c(&out.c);
                // Ne_c decreases as k loosens: too many effective parents ⇒ raise k.
                if ne > t_ne {
                    lo = k;
                } else {
                    hi = k;
                }
                if (ne - t_ne).abs() / t_ne < 0.01 {
                    break;
                }
            }
            let t = Instant::now();
            let out = solve_sexed(
                d.z.as_ref(),
                d.s,
                ridge,
                &d.b,
                &male,
                chosen_k,
                50_000,
                1e-8,
            );
            let secs = t.elapsed().as_secs_f64();
            emit(
                &format!("target_ne_{}", t_ne as u32),
                chosen_k,
                &out,
                secs,
                &mut csv,
            );
        }
    }

    let path = "artifacts/scaling_ne.csv";
    if let Err(e) = std::fs::create_dir_all("artifacts").and_then(|()| std::fs::write(path, &csv)) {
        eprintln!("[scaling_ne] could not write {path}: {e}");
    } else {
        eprintln!("[scaling_ne] wrote {path}");
    }
}
