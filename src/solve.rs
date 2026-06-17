//! Thin wrapper over Clarabel: turn a [`ConicProblem`] into a structured
//! outcome. No business logic here beyond mapping settings and reading the
//! solution fields the report cares about.

use crate::error::OcsError;
use crate::socp::ConicProblem;
use clarabel::solver::{DefaultSettingsBuilder, DefaultSolver, IPSolver, SolverStatus};

/// Solver knobs exposed by the spike.
#[derive(Clone, Copy, Debug)]
pub struct SolveConfig {
    /// Iteration cap (Clarabel default is 200).
    pub max_iter: u32,
    /// Wall-clock cap in seconds; `f64::INFINITY` disables it.
    pub time_limit: f64,
    /// Stream Clarabel's per-iteration table to stderr.
    pub verbose: bool,
}

impl Default for SolveConfig {
    fn default() -> Self {
        SolveConfig {
            max_iter: 200,
            time_limit: f64::INFINITY,
            verbose: false,
        }
    }
}

/// What the solve produced, distilled for reporting.
#[derive(Clone, Debug)]
pub struct SolveOutcome {
    /// Optimal contributions `c` (length `n`).
    pub c: Vec<f64>,
    /// Terminal solver status.
    pub status: SolverStatus,
    /// Iterations taken.
    pub iterations: u32,
    /// Clarabel's own measured solve time (seconds).
    pub solve_time: f64,
    /// Genetic gain `bᵀc = -obj_val`.
    pub gain: f64,
    /// Primal residual at termination.
    pub r_prim: f64,
    /// Dual residual at termination.
    pub r_dual: f64,
}

impl SolveOutcome {
    /// Whether Clarabel reached a usable optimum (full or reduced accuracy).
    pub fn is_solved(&self) -> bool {
        matches!(
            self.status,
            SolverStatus::Solved | SolverStatus::AlmostSolved
        )
    }
}

/// Run Clarabel on an assembled problem.
pub fn solve(prob: &ConicProblem, config: SolveConfig) -> Result<SolveOutcome, OcsError> {
    let settings = DefaultSettingsBuilder::<f64>::default()
        .verbose(config.verbose)
        .max_iter(config.max_iter)
        .time_limit(config.time_limit)
        .build()
        .map_err(|e| OcsError::SolverInit(format!("settings: {e:?}")))?;

    let mut solver = DefaultSolver::new(&prob.p, &prob.q, &prob.a, &prob.b, &prob.cones, settings)
        .map_err(|e| OcsError::SolverInit(format!("{e:?}")))?;

    solver.solve();

    let sol = &solver.solution;
    Ok(SolveOutcome {
        c: sol.x.clone(),
        status: sol.status,
        iterations: sol.iterations,
        solve_time: sol.solve_time,
        gain: -sol.obj_val,
        r_prim: sol.r_prim,
        r_dual: sol.r_dual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::socp::{Factor, build};
    use faer::Mat;

    #[test]
    fn known_two_variable_optimum() {
        // Same hand-checked instance as examples/spike.rs: n=2, G=I, b=[1,2],
        // k=0.6 → c₂ = (1+√0.2)/2 = 0.7236068.
        let n = 2;
        let l = Mat::<f64>::identity(n, n);
        let b = vec![1.0, 2.0];
        let prob = build(Factor::Cholesky(&l), &b, 0.6, 1.0, None);
        let out = solve(&prob, SolveConfig::default()).unwrap();
        assert!(out.is_solved());
        let c2_expected = (1.0 + 0.2_f64.sqrt()) / 2.0;
        assert!((out.c[1] - c2_expected).abs() < 1e-5);
        assert!((out.c[0] + out.c[1] - 1.0).abs() < 1e-6);
        assert!((out.gain - (1.0 * out.c[0] + 2.0 * out.c[1])).abs() < 1e-6);
    }

    #[test]
    fn per_candidate_cap_binds() {
        // Same instance, but cap every contribution at 0.5. With Σc=1 and n=2
        // the only feasible point is c=[0.5,0.5] (‖c‖²=0.5 ≤ 0.6), so the cap —
        // not the cone — sets the optimum. Exercises the c ≤ u path end to end.
        let n = 2;
        let l = Mat::<f64>::identity(n, n);
        let b = vec![1.0, 2.0];
        let u = vec![0.5, 0.5];
        let prob = build(Factor::Cholesky(&l), &b, 0.6, 1.0, Some(&u));
        let out = solve(&prob, SolveConfig::default()).unwrap();
        assert!(out.is_solved());
        assert!((out.c[0] - 0.5).abs() < 1e-5, "c0={}", out.c[0]);
        assert!((out.c[1] - 0.5).abs() < 1e-5, "c1={}", out.c[1]);
        assert!(out.c.iter().all(|&ci| ci <= 0.5 + 1e-7), "cap violated");
    }
}
