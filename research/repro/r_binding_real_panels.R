# The matrix-free solver against optiSel on REAL panels, both in one R session.
#
# Why this script exists: the timings in the manuscript's Table 2 come from a NumPy
# prototype that consumes the dense kinship matrix K and extracts K[S, S] by
# indexing. The Rust solver is matrix-free — it recomputes those entries from Z at
# O(|S| n m) per iteration — so its cost depends on the support in a way the dense
# prototype's does not. The exports in this directory only ever wrote K, so the
# matrix-free path had never been timed on real data. The R binding takes Z
# directly, which makes the measurement a single call.
#
# Instances and caps are reproduced exactly from wheat_export.R and mouse_export.R,
# including their deliberate choice of a bound relative to the gain-greedy plan.
# Conventions: optiSel constrains c' sKin c <= ub with sKin = G/2 + 1e-5 I, so the
# identical problem for ocs_solve() is k = 2*ub on G + 2e-5 I.
#
# Usage: Rscript research/repro/r_binding_real_panels.R [wheat|mouse|both]

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

report <- function(label, panel, z, bv, sx, s, ub, sKin) {
  n <- length(bv)
  ids <- rownames(sKin)
  phen <- data.frame(
    Indiv = ids, Born = 1L, Breed = label, BV = bv, Sex = sx,
    isCandidate = TRUE, stringsAsFactors = FALSE
  )
  cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)

  t_os <- system.time(
    fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )
  oc <- fit$parent$oc

  # Isolate marshalling from the solve: max_iter = 1 returns immediately, so the
  # difference is the R -> Rust copy of an n x m matrix.
  t_marshal <- system.time(
    ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male", max_iter = 1)
  )
  t_rs <- system.time(
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  )

  kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
  kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)

  cat(sprintf("\n=== %s: n=%d markers=%d, sexed, ub=%.8f ===\n", label, n, ncol(z), ub))
  cat(sprintf("optiSel   : gain=%.6f  coancestry=%.8f (%.3f%% of ub)  support=%3d  %.3fs\n",
              sum(bv * oc), kin_os, 100 * kin_os / ub, sum(oc > 1e-6), t_os[["elapsed"]]))
  cat(sprintf("ocsrs     : gain=%.6f  coancestry=%.8f (%.3f%% of ub)  support=%3d  %.3fs"
              , sum(bv * res$c), kin_rs, 100 * kin_rs / ub, length(res$support),
              t_rs[["elapsed"]]))
  cat(sprintf("  (of which %.3fs is the R->Rust copy)\n", t_marshal[["elapsed"]]))
  cat(sprintf("            iterations=%d  products=%d  status=%s\n",
              res$iterations, res$products, res$status))
  cat(sprintf("speed-up  : %.0fx   gain gap: %+.6f\n",
              t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4),
              sum(bv * res$c) - sum(bv * oc)))
  cat(sprintf("feasible  : sum=%.10f  min=%.2e  males=%.8f  within ub=%s\n",
              sum(res$c), min(res$c), sum(res$c[sx == "male"]),
              kin_rs <= ub * (1 + 1e-9)))
  invisible(NULL)
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
  ids <- paste0("L", seq_len(n))
  rownames(sKin) <- colnames(sKin) <- ids
  bv <- wheat.Y[, 1]
  # Wheat lines are autogamous; the two-group split exists only because optiSel
  # requires two sexes. It benchmarks both solvers on one well-defined problem.
  set.seed(1)
  sx <- rep(c("male", "female"), length.out = n)
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(n)
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  ub <- 0.12 * as.numeric(t(cg) %*% sKin %*% cg)
  report("wheat (CIMMYT)", panel = NULL, z = z, bv = bv, sx = sx, s = s, ub = ub, sKin = sKin)
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
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(length(ids))
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  ub <- 0.15 * as.numeric(t(cg) %*% sKin %*% cg)
  report("HS mouse (real sex)", panel = NULL, z = z, bv = bv, sx = sx, s = s, ub = ub, sKin = sKin)
}
