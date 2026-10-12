# Structured scaling at a contribution-Ne target; empirical support measurement.
#
# The synthetic generator behind Figure 1A is HWE with independent loci, so growing
# n at fixed m dilutes relatedness (min coancestry ~ 1/n → 0) and a bounded support
# is partly an artifact of a relaxing constraint. Here the candidates descend from a
# FIXED pool of parents (coalescent founders with real LD), so as n grows the
# population retains family structure. Panel-wise centering nevertheless forces
# the uniform-allocation kinship floor toward zero; this is not a fixed-reference
# genomic coancestry experiment. At contribution Ne_c=50 (ΔF proxy=1%),
# does the optimal support |S| stay bounded as n grows, or does it grow with n?
#
# Ne_c = 1/sum(c^2); ΔF_proxy = 1/(2 Ne_c), not realized genomic inbreeding.
# The legacy k_floor_uniform column is c_uniform' G c_uniform, a feasible
# reference value rather than a certified minimum. Group coancestry uses G/2.
# (matrix-free ||Z'c||^2/s plus ridge). With panel-wise centering it approaches
# zero even when candidates share parents; it is not a fixed-reference measure
# of persistent population structure.
#
# optiSel is run only at small n (it must form the dense n x n G); the dense-G GiB is
# reported as a theoretical array footprint, excluding solver workspace.
#
# Usage: Rscript research/repro/scaling_structured.R

source("research/repro/repro_common.R")
unlink(file.path(repro_out, "scaling_structured.csv"))
suppressMessages({
  library(AlphaSimR)
  library(ocsrs)
})
set.seed(1)

# Fixed base: 100 founders with coalescent LD and a fixed historical Ne; a fixed
# parent pool from which every candidate set is bred.
founderPop <- repro_founders("structured_founders", nInd = 100, nChr = 10,
                             segSites = 1100, seed = 1)
SP <- SimParam$new(founderPop)
SP$nThreads <- repro_sim_threads
SP$setSexes("yes_rand")
SP$addTraitA(nQtlPerChr = 100)
SP$addSnpChip(nSnpPerChr = 100) # 1000 SNP markers, matching Figure 1's m
parents <- newPop(founderPop)

ne_of <- function(c) 1 / sum(c^2)
quad_mf <- function(Z, s, c) sum((as.numeric(crossprod(Z, c)))^2) / s + 2e-5 * sum(c^2) # includes the solver ridge
target_ne <- 50
optisel_upto <- 2000L
ns <- c(1000L, 2000L, 5000L, 10000L, 20000L, 40000L)

rows <- list()
cat("structured scaling  (fixed 100-parent pool, contribution proxy=1% / Ne_c=50)\n")
for (n in ns) {
  off <- randCross(parents, nCrosses = n, nProgeny = 1)
  M <- pullSnpGeno(off)
  sx <- ifelse(off@sex == "M", "male", "female")
  bv <- as.numeric(off@gv[, 1])
  p <- colMeans(M) / 2
  Z <- sweep(M, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))

  # Uniform within-sex feasible allocation: a tight bracket, not the exact minimum.
  cu <- numeric(n)
  cu[sx == "male"] <- 0.5 / sum(sx == "male")
  cu[sx == "female"] <- 0.5 / sum(sx == "female")
  k_floor <- quad_mf(Z, s, cu)
  # Loose bracket: the 2-support truncation plan.
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(n)
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  k_loose <- quad_mf(Z, s, cg)

  # Bisect the cap to the target Ne.
  lo <- k_floor
  hi <- k_loose
  chosen <- 0.5 * (lo + hi)
  for (it in 1:40) {
    k <- 0.5 * (lo + hi)
    chosen <- k
    res <- checked_solve(Z, bv, k = k, s = s, ridge = 2e-5, male = sx == "male")
    ne <- ne_of(res$c)
    if (ne > target_ne) lo <- k else hi <- k
    if (abs(ne - target_ne) / target_ne < 0.01) break
  }
  require_target(ne, target_ne)
  k <- chosen
  invisible(checked_solve(Z, bv, k = k, s = s, ridge = 2e-5, male = sx == "male")) # warm
  t_rs <- system.time(res <- checked_solve(Z, bv, k = k, s = s, ridge = 2e-5, male = sx == "male"))
  ne <- ne_of(res$c)
  coan <- quad_mf(Z, s, res$c)
  g_gib <- n * n * 8 / 2^30

  # optiSel cross-check at small n (same optimum on structured data).
  os_txt <- ""
  if (n <= optisel_upto) {
    suppressMessages(library(optiSel))
    sKin <- (Z %*% t(Z)) / s / 2
    diag(sKin) <- diag(sKin) + 1e-5
    ids <- paste0("I", seq_len(n))
    rownames(sKin) <- colnames(sKin) <- ids
    phen <- data.frame(
      Indiv = ids, Born = 1L, Breed = "X", BV = bv, Sex = sx,
      isCandidate = TRUE, stringsAsFactors = FALSE
    )
    cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
    ub <- k / 2 # optiSel bounds c' sKin c = c'(G/2) c ≤ ub, ours is c'Gc ≤ k
    t_os <- system.time(
      fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
    )
    oc <- fit$parent$oc
    gap <- sum(bv * res$c) - sum(bv * oc)
    os_txt <- sprintf(
      "  | optiSel %.2fs (|S|=%d, gaingap %+.2e, speedup %.0fx)",
      t_os[["elapsed"]], sum(oc > 1e-6), gap, t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4)
    )
  }

  cat(sprintf(
    "n=%6d | Ne_c=%4.1f proxy_dF=%.2f%% | |S|=%4d | k_floor=%.4e | coan=%.4e | solve %6.3fs | dense G %6.2f GiB%s\n",
    n, ne, 50 / ne, length(res$support), k_floor, coan, t_rs[["elapsed"]], g_gib, os_txt
  ))
  rows[[length(rows) + 1]] <- data.frame(
    n = n, m = ncol(Z), target_ne = target_ne, ne = ne, ridge = 2e-5, status = res$status,
    simulation_threads = repro_sim_threads, proxy_delta_f_pct = 50 / ne, support = length(res$support),
    k_floor_uniform = k_floor, coancestry = coan, optisel_s = if (n <= optisel_upto) t_os[["elapsed"]] else NA_real_,
    gain_gap = if (n <= optisel_upto) gap else NA_real_, solve_s = t_rs[["elapsed"]], g_gib = g_gib
  )
}

out <- do.call(rbind, rows)
write.csv(out, file.path(repro_out, "scaling_structured.csv"), row.names = FALSE)
cat("\n[scaling_structured] wrote ", file.path(repro_out, "scaling_structured.csv"), "\n", sep = "")
