//! OCS × Clarabel spike driver.
//!
//! stdout carries machine-readable data (the `scale` CSV row); all human
//! progress goes to stderr. Subcommands:
//!
//! - `correctness` — n=50, solve, dump artifacts, assert feasibility.
//! - `frontier` — sweep k, assert gain is monotone in k, write CSV.
//! - `scale` — one (n, m, route) point; print one CSV row to stdout.
//! - `all` — orchestrate the whole sweep (spawns `scale` children under
//!   `/usr/bin/time -l` for peak RSS) and write REPORT.md.

use ocs_rs::report::{self, FrontierPoint, ScalingRecord, Verdict, VerdictInputs, markers_for};
use ocs_rs::socp::{self, Factor, Route};
use ocs_rs::solve::{self, SolveConfig};
use ocs_rs::{datagen, feasibility, grm, support_first};
use std::env;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::time::Instant;

const FEAS_TOL: f64 = 1e-6;
const FRONTIER_TOL: f64 = 1e-6;
const INITIAL_RIDGE: f64 = 1e-5;
const MAX_RIDGE_TRIES: u32 = 12;
const RIDGE_FLAG: f64 = 1e-3; // above this, conditioning is "investigate"
const MEM_GB: f64 = 36.0;
const ARTIFACTS: &str = "artifacts";

fn main() {
    std::process::exit(real_main());
}

fn real_main() -> i32 {
    let args: Vec<String> = env::args().collect();
    let cmd = args.get(1).map(String::as_str).unwrap_or("all");

    let result = match cmd {
        "correctness" => cmd_correctness(
            flag_usize(&args, "--n", 50),
            flag_usize(&args, "--m", 2000),
            flag_u64(&args, "--seed", 20240617),
            flag_f64(&args, "--k-frac", 0.6),
        ),
        "frontier" => cmd_frontier(
            flag_usize(&args, "--n", 250),
            flag_usize(&args, "--m", 4000),
            flag_u64(&args, "--seed", 20240617),
            flag_usize(&args, "--points", 14),
        ),
        "scale" => {
            let n = flag_usize(&args, "--n", 1000);
            cmd_scale(
                n,
                flag_usize(&args, "--m", markers_for(n)),
                flag_u64(&args, "--seed", 20240617),
                flag_f64(&args, "--k-frac", 0.6),
                flag_route(&args, "--route", Route::Cholesky),
                flag_u32(&args, "--max-iter", 200),
                flag_f64(&args, "--time-limit", 600.0),
            )
        }
        "all" => cmd_all(
            flag_u64(&args, "--seed", 20240617),
            flag_usize(&args, "--max-n", 10000),
            flag_f64(&args, "--k-frac", 0.6),
            flag_u32(&args, "--max-iter", 200),
            flag_f64(&args, "--time-limit", 600.0),
        ),
        // Re-derive REPORT.md from existing artifacts/ without re-solving.
        "report" => cmd_report(),
        // Head-to-head: Clarabel vs support-first on identical data.
        "compare" => {
            let n = flag_usize(&args, "--n", 2000);
            cmd_compare(
                n,
                flag_usize(&args, "--m", markers_for(n)),
                flag_u64(&args, "--seed", 20240617),
                flag_f64(&args, "--k-frac", 0.6),
            )
        }
        "-h" | "--help" | "help" => {
            eprintln!(
                "usage: ocs_rs [correctness|frontier|scale|all] [--n N] [--m M] [--seed S] \
                 [--k-frac F] [--route a|b] [--max-iter I] [--time-limit SEC] [--max-n N]"
            );
            Ok(())
        }
        other => {
            eprintln!("unknown command: {other}");
            return 2;
        }
    };

    match result {
        Ok(()) => 0,
        Err(code) => code,
    }
}

// ----------------------------------------------------------------------------
// Subcommands
// ----------------------------------------------------------------------------

/// n=50 correctness: solve, dump for external cross-check, assert feasibility.
fn cmd_correctness(n: usize, m: usize, seed: u64, k_frac: f64) -> Result<(), i32> {
    eprintln!("[correctness] n={n} m={m} seed={seed}");
    let ds = datagen::generate(n, m, seed);
    let (grm, l, tries) = grm::build_and_factor(&ds.z, ds.s, INITIAL_RIDGE, MAX_RIDGE_TRIES)
        .map_err(|e| {
            eprintln!("[correctness] factorization failed: {e}");
            1
        })?;
    let k = k_frac * mean_diag(&ds.z, ds.s);
    let prob = socp::build(Factor::Cholesky(&l), &ds.b, k, ds.s, None);
    let out = solve::solve(&prob, SolveConfig::default()).map_err(|e| {
        eprintln!("[correctness] solve error: {e}");
        1
    })?;

    let quad = grm::quad_form(&grm, &out.c);
    let feas = feasibility(&out.c, quad, k);
    let ok = feas.ok(FEAS_TOL) && out.is_solved();

    eprintln!(
        "[correctness] status={:?} iters={} ridge={:e} tries={}",
        out.status, out.iterations, grm.ridge, tries
    );
    eprintln!(
        "[correctness] sum_c-1={:.2e} min_c={:.2e} cᵀGc={:.6} k={:.6} margin={:.3e} gain={:.6}",
        feas.sum_err,
        feas.min_c,
        feas.quad,
        feas.k,
        feas.k - feas.quad,
        out.gain
    );

    // Dump artifacts for independent cross-check.
    let dir = Path::new(ARTIFACTS).join("correctness");
    fs::create_dir_all(&dir).map_err(io_fail)?;
    report::write_true_grm_csv(&dir.join("grm_true.csv"), &grm).map_err(io_fail)?;
    report::write_vector_csv(&dir.join("b.csv"), "gebv", &ds.b).map_err(io_fail)?;
    report::write_vector_csv(&dir.join("c.csv"), "contribution", &out.c).map_err(io_fail)?;
    let meta = format!(
        "n,{n}\nm,{m}\nseed,{seed}\nk,{k:.10e}\nridge,{:e}\nvanraden_s,{:.10e}\n\
         status,{:?}\niterations,{}\ngain,{:.10e}\nquad_cGc,{:.10e}\nsum_err,{:.3e}\nmin_c,{:.3e}\n\
         note,grm_true.csv is G WITHOUT ridge; Clarabel solved with G+ridge*I\n",
        grm.ridge, ds.s, out.status, out.iterations, out.gain, feas.quad, feas.sum_err, feas.min_c,
    );
    fs::write(dir.join("meta.csv"), meta).map_err(io_fail)?;
    eprintln!("[correctness] artifacts -> {}", dir.display());

    if ok {
        println!("correctness PASS (feasible, {:?})", out.status);
        Ok(())
    } else {
        println!(
            "correctness FAIL (status={:?}, feasible={})",
            out.status,
            feas.ok(FEAS_TOL)
        );
        Err(1)
    }
}

/// Frontier: sweep k tight→loose, assert gain monotone non-decreasing.
fn cmd_frontier(n: usize, m: usize, seed: u64, points: usize) -> Result<(), i32> {
    eprintln!("[frontier] n={n} m={m} seed={seed} points={points}");
    let ds = datagen::generate(n, m, seed);
    let (grm, l, _tries) = grm::build_and_factor(&ds.z, ds.s, INITIAL_RIDGE, MAX_RIDGE_TRIES)
        .map_err(|e| {
            eprintln!("[frontier] factorization failed: {e}");
            1
        })?;

    // k range: from just above the uniform-contribution kinship (always
    // feasible) up past the gain-greedy vertex kinship (constraint inactive).
    let (k_lo, k_hi) = kinship_range(&grm);
    eprintln!("[frontier] k in [{k_lo:.4}, {k_hi:.4}]");

    let mut fpts = Vec::with_capacity(points);
    for i in 0..points {
        let t = i as f64 / (points - 1).max(1) as f64;
        let k = k_lo + t * (k_hi - k_lo);
        let prob = socp::build(Factor::Cholesky(&l), &ds.b, k, ds.s, None);
        let out = solve::solve(&prob, SolveConfig::default()).map_err(|e| {
            eprintln!("[frontier] solve error: {e}");
            1
        })?;
        let quad = grm::quad_form(&grm, &out.c);
        let feas = feasibility(&out.c, quad, k);
        fpts.push(FrontierPoint {
            k,
            quad,
            gain: out.gain,
            feasible: feas.ok(FEAS_TOL) && out.is_solved(),
        });
        eprintln!(
            "[frontier] k={k:.4} -> cᵀGc={quad:.4} gain={:.4} {:?}",
            out.gain, out.status
        );
    }

    fs::create_dir_all(ARTIFACTS).map_err(io_fail)?;
    let mut csv = String::from("k,diversity_cGc,gain,feasible\n");
    for p in &fpts {
        csv.push_str(&format!(
            "{:.6e},{:.6e},{:.6e},{}\n",
            p.k, p.quad, p.gain, p.feasible
        ));
    }
    fs::write(Path::new(ARTIFACTS).join("frontier.csv"), csv).map_err(io_fail)?;

    let monotone = report::frontier_is_monotone(&fpts, FRONTIER_TOL);
    let violation = report::frontier_max_violation(&fpts);
    eprintln!("[frontier] monotone={monotone} max_violation={violation:.2e}");
    if monotone {
        println!("frontier PASS (gain monotone in k, max violation {violation:.2e})");
        Ok(())
    } else {
        println!("frontier FAIL (max violation {violation:.2e})");
        Err(1)
    }
}

/// One scaling point: emit exactly one CSV row on stdout.
fn cmd_scale(
    n: usize,
    m: usize,
    seed: u64,
    k_frac: f64,
    route: Route,
    max_iter: u32,
    time_limit: f64,
) -> Result<(), i32> {
    let cfg = SolveConfig {
        max_iter,
        time_limit,
        verbose: false,
    };
    let rec = solve_one(n, m, seed, k_frac, route, cfg);
    // stdout: the data row, and only the data row.
    println!("{}", rec.to_csv_row());
    if rec.solved && rec.feasible {
        Ok(())
    } else {
        Err(1)
    }
}

/// Full orchestration: scaling sweep (per-point peak RSS via child processes),
/// plus correctness + frontier, then the GO/NO-GO REPORT.md.
fn cmd_all(
    seed: u64,
    max_n: usize,
    k_frac: f64,
    max_iter: u32,
    time_limit: f64,
) -> Result<(), i32> {
    fs::create_dir_all(ARTIFACTS).map_err(io_fail)?;

    // 1) correctness + frontier (in-process, cheap).
    eprintln!("=== n=50 correctness ===");
    let _ = cmd_correctness(50, 2000, seed, k_frac); // non-fatal: still emit a report
    let correctness_margin = read_correctness_margin().unwrap_or(f64::NAN);

    eprintln!("=== frontier ===");
    let _ = cmd_frontier(250, 4000, seed, 14);
    let (frontier_monotone, frontier_violation) =
        read_frontier_result().unwrap_or((false, f64::NAN));

    // 2) scaling sweep plan.
    let exe = env::current_exe().map_err(io_fail)?;
    let mut plan: Vec<(usize, usize, Route)> = Vec::new();
    for &n in report::SWEEP_NS.iter().filter(|&&n| n <= max_n) {
        plan.push((n, markers_for(n), Route::Cholesky));
    }
    for &n in report::BOTH_ROUTES_NS.iter().filter(|&&n| n <= max_n) {
        plan.push((n, markers_for(n), Route::Raw));
    }

    let csv_path = Path::new(ARTIFACTS).join("scaling.csv");
    let mut csv_body = format!("{},peak_rss_mb\n", ScalingRecord::CSV_HEADER);
    fs::write(&csv_path, &csv_body).map_err(io_fail)?;

    let mut records: Vec<(ScalingRecord, Option<f64>)> = Vec::new();
    for (n, m, route) in plan {
        eprintln!(
            "=== scale n={n} m={m} route={} (child + /usr/bin/time -l) ===",
            route.tag()
        );
        let (rec, rss) = run_scale_child(&exe, n, m, seed, k_frac, route, max_iter, time_limit);
        eprintln!(
            "    -> status={} solved={} feasible={} iters={} solve={:.2}s rss={}",
            rec.status,
            rec.solved,
            rec.feasible,
            rec.iterations,
            rec.t_solve_s,
            rss.map(|v| format!("{v:.0} MB"))
                .unwrap_or_else(|| "n/a".into()),
        );
        // Persist incrementally so a later OOM cannot lose earlier rows.
        csv_body.push_str(&rec.to_csv_row());
        csv_body.push(',');
        csv_body.push_str(&rss.map(|v| format!("{v:.1}")).unwrap_or_default());
        csv_body.push('\n');
        fs::write(&csv_path, &csv_body).map_err(io_fail)?;
        records.push((rec, rss));
    }

    // 3) verdict.
    let inputs = VerdictInputs {
        records: &records,
        frontier_monotone,
        frontier_violation,
        correctness_margin,
        tol: FEAS_TOL,
        ridge_flag: RIDGE_FLAG,
        mem_gb: MEM_GB,
    };
    let (verdict, md) = report::generate_verdict(&inputs);
    fs::write("REPORT.md", &md).map_err(io_fail)?;
    eprintln!("=== wrote REPORT.md ===");
    println!(
        "VERDICT: {}",
        if verdict == Verdict::Go {
            "GO"
        } else {
            "NO-GO / INVESTIGATE"
        }
    );
    Ok(())
}

/// Re-derive REPORT.md from existing artifacts/ without re-solving anything.
fn cmd_report() -> Result<(), i32> {
    let records = read_scaling_records().ok_or_else(|| {
        eprintln!("[report] cannot read artifacts/scaling.csv — run `all` or `scale` first");
        1
    })?;
    let correctness_margin = read_correctness_margin().unwrap_or(f64::NAN);
    let (frontier_monotone, frontier_violation) =
        read_frontier_result().unwrap_or((false, f64::NAN));
    let inputs = VerdictInputs {
        records: &records,
        frontier_monotone,
        frontier_violation,
        correctness_margin,
        tol: FEAS_TOL,
        ridge_flag: RIDGE_FLAG,
        mem_gb: MEM_GB,
    };
    let (verdict, md) = report::generate_verdict(&inputs);
    fs::write("REPORT.md", &md).map_err(io_fail)?;
    eprintln!("[report] regenerated REPORT.md from artifacts/scaling.csv");
    println!(
        "VERDICT: {}",
        if verdict == Verdict::Go {
            "GO"
        } else {
            "NO-GO / INVESTIGATE"
        }
    );
    Ok(())
}

/// Head-to-head: the conic IPM (Clarabel, Route A) vs the support-first solver
/// on identical data. Reports per-stage timing, gain agreement, support sizes.
fn cmd_compare(n: usize, m: usize, seed: u64, k_frac: f64) -> Result<(), i32> {
    eprintln!("[compare] n={n} m={m} seed={seed}");
    let (ds, t_datagen) = timed(|| datagen::generate(n, m, seed));
    let k = k_frac * mean_diag(&ds.z, ds.s);

    // --- Clarabel, Route A (the path the spike validated) ---
    let (grm, t_grm) = timed(|| grm::Grm::build(&ds.z, ds.s, INITIAL_RIDGE));
    let (l_res, t_chol) = timed(|| grm.cholesky_lower());
    let l = l_res.map_err(|e| {
        eprintln!("[compare] cholesky: {e}");
        1
    })?;
    let (prob, t_assemble) = timed(|| socp::build(Factor::Cholesky(&l), &ds.b, k, ds.s, None));
    let (clar_res, t_clar) = timed(|| solve::solve(&prob, SolveConfig::default()));
    let clar = clar_res.map_err(|e| {
        eprintln!("[compare] clarabel: {e}");
        1
    })?;
    let clar_supp = clar.c.iter().filter(|&&c| c > 1e-6).count();
    let t_clar_pipeline = t_grm + t_chol + t_assemble + t_clar;

    // --- support-first (no GRM, no Cholesky; matrix-free over Z) ---
    let (sf, t_sf) =
        timed(|| support_first::solve(&ds.z, ds.s, INITIAL_RIDGE, &ds.b, k, 4000, 1e-7));

    let gain_match = (sf.gain - clar.gain).abs();
    let speedup_solve = if t_sf > 0.0 {
        t_clar / t_sf
    } else {
        f64::INFINITY
    };
    let speedup_pipe = if t_sf > 0.0 {
        t_clar_pipeline / t_sf
    } else {
        f64::INFINITY
    };

    // stdout: the comparison.
    println!("n={n}  m={m}  k={k:.5}  (datagen {t_datagen:.3}s, shared)\n");
    println!("  stage                 Clarabel(RouteA)    support-first");
    println!("  GRM build             {t_grm:>10.3}s          —");
    println!("  Cholesky              {t_chol:>10.3}s          —");
    println!("  assemble              {t_assemble:>10.3}s          —");
    println!(
        "  solve                 {t_clar:>10.3}s   {t_sf:>10.3}s   ({} G·c products)",
        sf.products
    );
    println!("  ----------------------------------------------------------");
    println!("  total (excl datagen)  {t_clar_pipeline:>10.3}s   {t_sf:>10.3}s");
    println!(
        "  gain                  {:>11.6}    {:>11.6}   (Δ={gain_match:.2e})",
        clar.gain, sf.gain
    );
    println!(
        "  |support|             {clar_supp:>11}    {:>11}",
        sf.support.len()
    );
    println!(
        "  feasible cᵀGc≤k       {:>11}    {:>11}",
        clar.is_solved(),
        sf.quad <= k + 1e-6
    );
    println!(
        "  status                {:>11?}    {:?}",
        clar.status, sf.status
    );
    println!("\n  speedup: solve {speedup_solve:.0}×   |   full pipeline {speedup_pipe:.1}×");

    if gain_match < 1e-6 && sf.status == support_first::SfStatus::Solved {
        Ok(())
    } else {
        eprintln!("[compare] WARNING: gain mismatch {gain_match:.2e} or non-Solved status");
        Err(1)
    }
}

// ----------------------------------------------------------------------------
// Core pipeline (timed)
// ----------------------------------------------------------------------------

fn solve_one(
    n: usize,
    m: usize,
    seed: u64,
    k_frac: f64,
    route: Route,
    cfg: SolveConfig,
) -> ScalingRecord {
    let (ds, t_datagen) = timed(|| datagen::generate(n, m, seed));
    let k = k_frac * mean_diag(&ds.z, ds.s);

    let dense_work_gb = match route {
        Route::Cholesky => 8.0 * (n as f64 * m as f64 + 2.0 * n as f64 * n as f64) / 1e9,
        Route::Raw => 8.0 * (2.0 * n as f64 * m as f64) / 1e9,
    };

    match route {
        Route::Cholesky => {
            let (grm, t_grm) = timed(|| grm::Grm::build(&ds.z, ds.s, INITIAL_RIDGE));
            // Factor with ridge escalation; escalation cost belongs to
            // "obtaining a usable factor".
            let t0 = Instant::now();
            let mut grm = grm;
            let mut tries = 0u32;
            let mut ridge = INITIAL_RIDGE;
            let l = loop {
                match grm.cholesky_lower() {
                    Ok(l) => break Some(l),
                    Err(_) if tries < MAX_RIDGE_TRIES => {
                        tries += 1;
                        ridge *= 10.0;
                        grm = grm::Grm::build(&ds.z, ds.s, ridge);
                    }
                    Err(_) => break None,
                }
            };
            let t_factor = t0.elapsed().as_secs_f64();

            let Some(l) = l else {
                return failed_record(n, m, route, k, t_datagen, dense_work_gb, "ridge_exhausted");
            };

            let (prob, t_assemble) =
                timed(|| socp::build(Factor::Cholesky(&l), &ds.b, k, ds.s, None));
            assemble_and_solve(
                ProblemCtx {
                    n,
                    m,
                    route,
                    k,
                    ridge,
                    tries,
                    t_datagen,
                    t_grm,
                    t_factor,
                    t_assemble,
                    dense_work_gb,
                },
                prob,
                cfg,
                |c| grm::quad_form(&grm, c),
            )
        }
        Route::Raw => {
            // No GRM, no Cholesky — that is the route's whole point.
            let (prob, t_assemble) =
                timed(|| socp::build(Factor::Raw(&ds.z), &ds.b, k, ds.s, None));
            assemble_and_solve(
                ProblemCtx {
                    n,
                    m,
                    route,
                    k,
                    ridge: 0.0,
                    tries: 0,
                    t_datagen,
                    t_grm: 0.0,
                    t_factor: 0.0,
                    t_assemble,
                    dense_work_gb,
                },
                prob,
                cfg,
                |c| grm::quad_form_z(&ds.z, ds.s, c),
            )
        }
    }
}

struct ProblemCtx {
    n: usize,
    m: usize,
    route: Route,
    k: f64,
    ridge: f64,
    tries: u32,
    t_datagen: f64,
    t_grm: f64,
    t_factor: f64,
    t_assemble: f64,
    dense_work_gb: f64,
}

fn assemble_and_solve(
    ctx: ProblemCtx,
    prob: socp::ConicProblem,
    cfg: SolveConfig,
    quad_of: impl Fn(&[f64]) -> f64,
) -> ScalingRecord {
    let a_rows = prob.a.m;
    let a_nnz = prob.a_nnz;
    let t0 = Instant::now();
    let solved = solve::solve(&prob, cfg);
    let t_solve = t0.elapsed().as_secs_f64();

    match solved {
        Ok(out) => {
            let quad = quad_of(&out.c);
            let feas = feasibility(&out.c, quad, ctx.k);
            ScalingRecord {
                n: ctx.n,
                m: ctx.m,
                route: ctx.route,
                status: format!("{:?}", out.status),
                solved: out.is_solved(),
                feasible: feas.ok(FEAS_TOL),
                iterations: out.iterations,
                ridge: ctx.ridge,
                ridge_tries: ctx.tries,
                t_datagen_s: ctx.t_datagen,
                t_grm_s: ctx.t_grm,
                t_factor_s: ctx.t_factor,
                t_assemble_s: ctx.t_assemble,
                t_solve_s: t_solve,
                clarabel_solve_s: out.solve_time,
                gain: out.gain,
                quad,
                k: ctx.k,
                sum_err: feas.sum_err,
                min_c: feas.min_c,
                a_rows,
                a_nnz,
                dense_work_gb: ctx.dense_work_gb,
            }
        }
        Err(e) => {
            eprintln!("[scale] solver error: {e}");
            let mut r = failed_record(
                ctx.n,
                ctx.m,
                ctx.route,
                ctx.k,
                ctx.t_datagen,
                ctx.dense_work_gb,
                "solver_error",
            );
            r.t_grm_s = ctx.t_grm;
            r.t_factor_s = ctx.t_factor;
            r.t_assemble_s = ctx.t_assemble;
            r.t_solve_s = t_solve;
            r.a_rows = a_rows;
            r.a_nnz = a_nnz;
            r
        }
    }
}

fn failed_record(
    n: usize,
    m: usize,
    route: Route,
    k: f64,
    t_datagen: f64,
    dense_work_gb: f64,
    why: &str,
) -> ScalingRecord {
    ScalingRecord {
        n,
        m,
        route,
        status: why.to_string(),
        solved: false,
        feasible: false,
        iterations: 0,
        ridge: 0.0,
        ridge_tries: 0,
        t_datagen_s: t_datagen,
        t_grm_s: 0.0,
        t_factor_s: 0.0,
        t_assemble_s: 0.0,
        t_solve_s: 0.0,
        clarabel_solve_s: 0.0,
        gain: 0.0,
        quad: 0.0,
        k,
        sum_err: 0.0,
        min_c: 0.0,
        a_rows: 0,
        a_nnz: 0,
        dense_work_gb,
    }
}

// ----------------------------------------------------------------------------
// Child-process scaling run (peak RSS via /usr/bin/time -l)
// ----------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn run_scale_child(
    exe: &Path,
    n: usize,
    m: usize,
    seed: u64,
    k_frac: f64,
    route: Route,
    max_iter: u32,
    time_limit: f64,
) -> (ScalingRecord, Option<f64>) {
    let route_arg = match route {
        Route::Cholesky => "a",
        Route::Raw => "b",
    };
    let output = Command::new("/usr/bin/time")
        .arg("-l")
        .arg(exe)
        .arg("scale")
        .arg("--n")
        .arg(n.to_string())
        .arg("--m")
        .arg(m.to_string())
        .arg("--seed")
        .arg(seed.to_string())
        .arg("--k-frac")
        .arg(k_frac.to_string())
        .arg("--route")
        .arg(route_arg)
        .arg("--max-iter")
        .arg(max_iter.to_string())
        .arg("--time-limit")
        .arg(time_limit.to_string())
        .output();

    let output = match output {
        Ok(o) => o,
        Err(e) => {
            eprintln!("    !! failed to spawn child: {e}");
            return (
                failed_record(n, m, route, f64::NAN, 0.0, 0.0, "spawn_failed"),
                None,
            );
        }
    };

    let stderr = String::from_utf8_lossy(&output.stderr);
    let rss = parse_peak_rss_mb(&stderr);

    // The child prints exactly one CSV row to stdout on success.
    let stdout = String::from_utf8_lossy(&output.stdout);
    let row = stdout.lines().rev().find(|l| !l.trim().is_empty());
    let rec = row
        .and_then(ScalingRecord::from_csv_row)
        .unwrap_or_else(|| {
            // No parseable row: the child was likely OOM-killed.
            failed_record(n, m, route, f64::NAN, 0.0, 0.0, "killed_or_oom")
        });
    (rec, rss)
}

/// Parse "maximum resident set size" (bytes, macOS) from `/usr/bin/time -l`.
fn parse_peak_rss_mb(time_stderr: &str) -> Option<f64> {
    for line in time_stderr.lines() {
        if line.contains("maximum resident set size") {
            let first = line.split_whitespace().next()?;
            let bytes: f64 = first.parse().ok()?;
            return Some(bytes / 1_048_576.0);
        }
    }
    None
}

// ----------------------------------------------------------------------------
// Numeric helpers
// ----------------------------------------------------------------------------

/// Mean genomic self-relationship `(1/n) Σ_i (Σ_j Z[i,j]²)/s`, i.e. the mean of
/// the GRM diagonal, computed without forming G.
fn mean_diag(z: &faer::Mat<f64>, s: f64) -> f64 {
    let n = z.nrows();
    let m = z.ncols();
    let mut acc = 0.0;
    for i in 0..n {
        let mut row = 0.0;
        for j in 0..m {
            let v = z[(i, j)];
            row += v * v;
        }
        acc += row / s;
    }
    acc / n as f64
}

/// Frontier k endpoints from the (ridged) GRM: uniform-contribution kinship up
/// to a little past the gain-greedy vertex kinship.
fn kinship_range(grm: &grm::Grm) -> (f64, f64) {
    let n = grm.n;
    // uniform c = 1/n: cᵀGc = (1/n²) Σᵢⱼ G_ij (true G).
    let mut sum_all = 0.0;
    let mut max_diag = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            sum_all += grm.g[(i, j)];
        }
        let dii = grm.g[(i, i)] - grm.ridge;
        max_diag = max_diag.max(dii);
    }
    sum_all -= grm.ridge * n as f64; // remove ridge from the diagonal
    let k_uniform = sum_all / (n as f64 * n as f64);
    let k_lo = k_uniform * 1.05;
    let k_hi = max_diag * 1.10;
    (k_lo, k_hi.max(k_lo * 1.5))
}

fn timed<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let t0 = Instant::now();
    let v = f();
    (v, t0.elapsed().as_secs_f64())
}

// ----------------------------------------------------------------------------
// Artifact readers (for the `all` orchestrator)
// ----------------------------------------------------------------------------

fn read_scaling_records() -> Option<Vec<(ScalingRecord, Option<f64>)>> {
    let txt = fs::read_to_string(Path::new(ARTIFACTS).join("scaling.csv")).ok()?;
    let mut out = Vec::new();
    for line in txt.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        // peak_rss_mb is the final column appended by the orchestrator; split it
        // off and parse the rest with the record's own (round-trip-tested) parser.
        let (row, rss_str) = line.rsplit_once(',')?;
        let rec = ScalingRecord::from_csv_row(row)?;
        let rss = rss_str.trim().parse::<f64>().ok();
        out.push((rec, rss));
    }
    (!out.is_empty()).then_some(out)
}

fn read_correctness_margin() -> Option<f64> {
    let path = Path::new(ARTIFACTS).join("correctness").join("meta.csv");
    let txt = fs::read_to_string(path).ok()?;
    let mut k = None;
    let mut quad = None;
    for line in txt.lines() {
        let (key, val) = line.split_once(',')?;
        match key {
            "k" => k = val.parse::<f64>().ok(),
            "quad_cGc" => quad = val.parse::<f64>().ok(),
            _ => {}
        }
    }
    Some(k? - quad?)
}

fn read_frontier_result() -> Option<(bool, f64)> {
    let path = Path::new(ARTIFACTS).join("frontier.csv");
    let txt = fs::read_to_string(path).ok()?;
    let mut pts = Vec::new();
    for line in txt.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() < 4 {
            continue;
        }
        pts.push(FrontierPoint {
            k: f[0].parse().ok()?,
            quad: f[1].parse().ok()?,
            gain: f[2].parse().ok()?,
            feasible: f[3].trim().parse().unwrap_or(false),
        });
    }
    Some((
        report::frontier_is_monotone(&pts, FRONTIER_TOL),
        report::frontier_max_violation(&pts),
    ))
}

// ----------------------------------------------------------------------------
// Tiny arg parsing (edge of the program; core logic stays testable)
// ----------------------------------------------------------------------------

fn flag<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
}

fn flag_usize(args: &[String], name: &str, default: usize) -> usize {
    flag(args, name)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}
fn flag_u32(args: &[String], name: &str, default: u32) -> u32 {
    flag(args, name)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}
fn flag_u64(args: &[String], name: &str, default: u64) -> u64 {
    flag(args, name)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}
fn flag_f64(args: &[String], name: &str, default: f64) -> f64 {
    flag(args, name)
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}
fn flag_route(args: &[String], name: &str, default: Route) -> Route {
    match flag(args, name) {
        Some("a") | Some("A") => Route::Cholesky,
        Some("b") | Some("B") => Route::Raw,
        _ => default,
    }
}

fn io_fail(e: std::io::Error) -> i32 {
    eprintln!("io error: {e}");
    1
}
