# Multi-generation validation (review point 8, the domain-standard check): recurrent
# selection over generations, support-first OCS at a fixed operating point (ΔF=1%,
# Ne=50) against truncation selection. The expected, textbook result: OCS gives up a
# little gain to hold the rate of coancestry down, truncation runs faster gain into
# runaway coancestry. Criterion is the true genetic value (isolates the OCS mechanism
# from GEBV noise). support-first is called through the shipped R binding each
# generation, so this also exercises the deployed tool end to end.
#
# Usage: Rscript research/repro/multigen_validation.R
suppressMessages({
  library(AlphaSimR)
  library(ocsrs)
})

n_gen <- 12L
n_rep <- 3L
n_progeny <- 250L      # population size held constant per generation
target_ne <- 50        # OCS operating point (ΔF = 1%)
n_trunc_per_sex <- 25L # truncation: top-25 of each sex by gv

ne_of <- function(c) 1 / sum(c^2)

# One generation of OCS mating: bisect the cap to Ne, sample crosses ∝ contributions.
ocs_cross <- function(pop) {
  M <- pullSnpGeno(pop)
  p <- colMeans(M) / 2
  z <- sweep(M, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))
  bv <- as.numeric(gv(pop)[, 1])
  sx <- pop@sex
  male <- sx == "M"
  iM <- which(male)[which.max(bv[male])]
  iF <- which(!male)[which.max(bv[!male])]
  n <- length(bv)
  sKin <- tcrossprod(z) / s / 2
  diag(sKin) <- diag(sKin) + 1e-5
  cg <- numeric(n); cg[iM] <- 0.5; cg[iF] <- 0.5
  ub_hi <- as.numeric(t(cg) %*% sKin %*% cg)
  cu <- numeric(n)
  cu[male] <- 0.5 / sum(male); cu[!male] <- 0.5 / sum(!male)
  ub_lo <- as.numeric(t(cu) %*% sKin %*% cu)
  ub <- 0.5 * (ub_lo + ub_hi)
  res <- NULL
  for (it in 1:25) {
    ub <- 0.5 * (ub_lo + ub_hi)
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = male)
    ne <- ne_of(res$c)
    if (ne > target_ne) ub_lo <- ub else ub_hi <- ub
    if (abs(ne - target_ne) / target_ne < 0.02) break
  }
  c <- pmax(res$c, 0)
  # Sample n_progeny crosses: dam ∝ c on females, sire ∝ c on males.
  f_idx <- which(!male); m_idx <- which(male)
  fw <- c[f_idx]; mw <- c[m_idx]
  dams <- sample(f_idx, n_progeny, replace = TRUE, prob = fw)
  sires <- sample(m_idx, n_progeny, replace = TRUE, prob = mw)
  makeCross(pop, cbind(dams, sires))
}

# One generation of truncation mating: top-k of each sex by gv, crossed at random.
trunc_cross <- function(pop) {
  bv <- as.numeric(gv(pop)[, 1])
  male <- pop@sex == "M"
  f_top <- which(!male)[order(bv[!male], decreasing = TRUE)][seq_len(n_trunc_per_sex)]
  m_top <- which(male)[order(bv[male], decreasing = TRUE)][seq_len(n_trunc_per_sex)]
  dams <- sample(f_top, n_progeny, replace = TRUE)
  sires <- sample(m_top, n_progeny, replace = TRUE)
  makeCross(pop, cbind(dams, sires))
}

# Population summaries: mean genetic value (gain) and observed heterozygosity (the
# fraction of heterozygous genotype calls). Heterozygosity is base-free and falls as
# inbreeding accumulates and loci fix — unlike a VanRaden GRM re-centred each
# generation, whose mean coancestry is ~0 by construction.
summ <- function(pop) {
  M <- pullSnpGeno(pop)
  list(gain = mean(gv(pop)[, 1]),
       het = mean(M == 1))
}

rows <- list()
for (rep in seq_len(n_rep)) {
  set.seed(1000 + rep)
  founderPop <- runMacs(nInd = 120, nChr = 10, segSites = 1100, species = "GENERIC")
  SP <- SimParam$new(founderPop)
  SP$setSexes("yes_rand")
  SP$addTraitA(nQtlPerChr = 100)
  SP$addSnpChip(nSnpPerChr = 100)
  base <- newPop(founderPop)
  for (method in c("ocs", "truncation")) {
    pop <- base
    for (g in 0:n_gen) {
      st <- summ(pop)
      rows[[length(rows) + 1]] <- data.frame(
        rep = rep, method = method, gen = g, gain = st$gain, het = st$het
      )
      if (g < n_gen) {
        pop <- if (method == "ocs") ocs_cross(pop) else trunc_cross(pop)
      }
    }
    cat(sprintf("rep %d %-10s done (gen %d)\n", rep, method, n_gen))
    flush.console()
  }
}

out <- do.call(rbind, rows)
dir.create("artifacts", showWarnings = FALSE)
write.csv(out, "artifacts/multigen_validation.csv", row.names = FALSE)

# Trajectories averaged over reps.
agg <- aggregate(cbind(gain, het) ~ method + gen, out, mean)
cat("\n=== gain and heterozygosity by generation (mean over reps) ===\n")
for (m in c("ocs", "truncation")) {
  cat(sprintf("%-10s gain: %s\n", m,
              paste(sprintf("%.2f", agg$gain[agg$method == m]), collapse = " ")))
  cat(sprintf("%-10s het : %s\n", m,
              paste(sprintf("%.3f", agg$het[agg$method == m]), collapse = " ")))
}
cat("\n[multigen] wrote artifacts/multigen_validation.csv\n")
