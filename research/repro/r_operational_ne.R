# Operational-point re-benchmark: support-first vs optiSel on REAL panels, at a
# breeder-defensible operating point (target effective size Ne / rate of inbreeding
# ΔF), not at the loose cap that produced the headline. Answers review point 1 on
# real data.
#
#   Ne_c = 1 / sum(c^2)   (effective number of contributing parents)
#   ΔF   ≈ 1 / (2 Ne_c)
#
# Conventions identical to r_binding_real_panels.R: optiSel constrains c' sKin c <= ub
# with sKin = G/2 + 1e-5 I, so ocs_solve() solves the identical problem with
# k = 2*ub on G + 2e-5 I. The cap ub is tuned (bisection on the solution's Ne) to
# each target, then both solvers are timed at that operating point.
#
# Reported: R BLAS is OpenBLAS (multithreaded), so the optiSel baseline is NOT on a
# slow reference BLAS — the speed-up is conservative w.r.t. BLAS choice.
#
# Usage: Rscript research/repro/r_operational_ne.R [wheat|mouse|both]

suppressMessages({
  library(BGLR)
  library(optiSel)
  library(ocsrs)
})

which_panel <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "both"
}

targets_ne <- c(25, 50, 100) # ΔF = 2%, 1%, 0.5%
ne_of <- function(c) 1 / sum(c^2)

rows <- list()

bench_panel <- function(label, z, bv, sx, s, sKin) {
  n <- length(bv)
  ids <- rownames(sKin)
  phen <- data.frame(
    Indiv = ids, Born = 1L, Breed = "X", BV = bv, Sex = sx,
    isCandidate = TRUE, stringsAsFactors = FALSE
  )
  cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)

  # Brackets on ub: the 2-support truncation plan (Ne~2, loose) and the uniform
  # within-sex allocation (high Ne, tight but always feasible).
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(n)
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  ub_loose <- as.numeric(t(cg) %*% sKin %*% cg)
  cu <- numeric(n)
  cu[sx == "male"] <- 0.5 / sum(sx == "male")
  cu[sx == "female"] <- 0.5 / sum(sx == "female")
  ub_tight <- as.numeric(t(cu) %*% sKin %*% cu)

  cat(sprintf(
    "\n=== %s: n=%d markers=%d ===  (ub bracket [%.3e, %.3e])\n",
    label, n, ncol(z), ub_tight, ub_loose
  ))

  for (t_ne in targets_ne) {
    lo <- ub_tight
    hi <- ub_loose
    chosen <- 0.5 * (lo + hi)
    for (it in 1:40) {
      ub <- 0.5 * (lo + hi)
      chosen <- ub
      res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
      ne <- ne_of(res$c)
      if (ne > t_ne) lo <- ub else hi <- ub # Ne decreases as ub loosens
      if (abs(ne - t_ne) / t_ne < 0.01) break
    }
    ub <- chosen

    # Timed head-to-head at the operational cap.
    t_os <- system.time(
      fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
    )
    oc <- fit$parent$oc
    t_rs <- system.time(
      res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
    )

    ne <- ne_of(res$c)
    df <- 50 / ne
    kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
    kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)
    su <- t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4)
    gap <- sum(bv * res$c) - sum(bv * oc)

    cat(sprintf(
      "  Ne=%3.0f dF=%.2f%% | |S|=%3d (optiSel %3d) | optiSel %6.3fs  ocsrs %6.3fs  speedup %5.0fx | gaingap %+.2e  coan=%.5f\n",
      ne, df, length(res$support), sum(oc > 1e-6),
      t_os[["elapsed"]], t_rs[["elapsed"]], su, gap, kin_rs
    ))

    rows[[length(rows) + 1]] <<- data.frame(
      panel = label, target_ne = t_ne, achieved_ne = ne, delta_f_pct = df,
      ub = ub, support_ocsrs = length(res$support), support_optisel = sum(oc > 1e-6),
      optisel_s = t_os[["elapsed"]], ocsrs_s = t_rs[["elapsed"]], speedup = su,
      gain_gap = gap, coancestry = kin_rs, stringsAsFactors = FALSE
    )
  }
}

if (which_panel %in% c("wheat", "both")) {
  data(wheat)
  x <- wheat.X # 599 x 1279 DArT 0/1
  p <- colMeans(x)
  z <- sweep(x, 2, p)
  s <- sum(p * (1 - p))
  g <- (z %*% t(z)) / s
  sKin <- g / 2
  diag(sKin) <- diag(sKin) + 1e-5
  n <- nrow(x)
  rownames(sKin) <- colnames(sKin) <- paste0("L", seq_len(n))
  bv <- wheat.Y[, 1]
  set.seed(1)
  sx <- rep(c("male", "female"), length.out = n)
  bench_panel("wheat (CIMMYT)", z, bv, sx, s, sKin)
}

if (which_panel %in% c("mouse", "both")) {
  data(mice, package = "BGLR")
  x <- mice.X
  ids <- rownames(x)
  gender <- as.character(mice.pheno$GENDER)
  sx <- ifelse(gender %in% c("M", "Male", "male", "1"), "male", "female")
  bv <- as.numeric(mice.pheno$Obesity.BMI)
  ok <- !is.na(bv)
  x <- x[ok, ]
  ids <- ids[ok]
  sx <- sx[ok]
  bv <- bv[ok]
  p <- colMeans(x) / 2
  z <- sweep(x, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))
  g <- tcrossprod(z) / s
  sKin <- g / 2
  diag(sKin) <- diag(sKin) + 1e-5
  rownames(sKin) <- colnames(sKin) <- ids
  bench_panel("HS mouse (real sex)", z, bv, sx, s, sKin)
}

out <- do.call(rbind, rows)
dir.create("artifacts", showWarnings = FALSE)
write.csv(out, "artifacts/operational_ne_realpanels.csv", row.names = FALSE)
cat("\n[operational_ne] wrote artifacts/operational_ne_realpanels.csv\n")
