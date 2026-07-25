# One clean, serial measurement of the shipped matrix-free solver vs optiSel across
# every reproducible panel, in a single R session so provenance is uniform. Nothing
# else must run concurrently: optiSel timings are sensitive to CPU/memory contention
# (a parallel job inflates them severalfold), so this is deliberately one job.
#
# For each panel it prints optiSel time, the ocsrs end-to-end time, the R->Rust copy
# alone (max_iter=1), support, gain and realised coancestry. Conventions match the
# export scripts: optiSel constrains c' sKin c <= ub on sKin = G/2 + 1e-5 I, so the
# identical problem for ocs_solve is k = 2*ub with ridge = 2e-5.
#
# Usage: Rscript research/repro/r_binding_bench_all.R

suppressMessages({
  library(BGLR)
  library(optiSel)
  library(ocsrs)
  library(data.table)
})

bench <- function(label, z, bv, sx, s, ub, sKin) {
  ids <- rownames(sKin)
  phen <- data.frame(
    Indiv = ids, Born = 1L, Breed = "x", BV = bv, Sex = sx,
    isCandidate = TRUE, stringsAsFactors = FALSE
  )
  cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
  t_os <- system.time(
    fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )
  oc <- fit$parent$oc
  t_copy <- system.time(
    ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male", max_iter = 1)
  )
  t_rs <- system.time(
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  )
  kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)
  cat(sprintf(
    "%-26s n=%-5d m=%-6d | optiSel %8.3fs | ocsrs %7.3fs (copy %6.3fs) | support %3d | gain %+.5f (gap %+.5f) | %.0fx\n",
    label, length(bv), ncol(z), t_os[["elapsed"]], t_rs[["elapsed"]], t_copy[["elapsed"]],
    length(res$support), sum(bv * res$c), sum(bv * res$c) - sum(bv * oc),
    t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4)
  ))
  flush.console()
}

# --- synthetic structured populations (tight cap = large support, the hard case) ---
for (n in c(1000L, 2000L, 5000L)) {
  m <- 500L
  set.seed(1)
  p <- runif(m, 0.05, 0.5)
  n_fam <- max(2L, n %/% 20L)
  base <- matrix(rbinom(n_fam * m, 2, rep(p, each = n_fam)), n_fam, m)
  fam <- sample(n_fam, n, replace = TRUE)
  mm <- base[fam, ]
  idx <- which(matrix(runif(as.numeric(n) * m), n, m) < 0.05)
  mm[idx] <- rbinom(length(idx), 2, p[((idx - 1L) %/% n) + 1L])
  z <- sweep(mm, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))
  g <- tcrossprod(z) / s
  sKin <- g / 2
  diag(sKin) <- diag(sKin) + 1e-5
  ids <- paste0("I", seq_len(n))
  rownames(sKin) <- colnames(sKin) <- ids
  bv <- rnorm(n)
  sx <- sample(rep(c("male", "female"), length.out = n))
  ub <- 1.04 * mean(sKin)
  bench("synthetic (structured)", z, bv, sx, s, ub, sKin)
}

# --- CIMMYT wheat ---
data(wheat)
x <- wheat.X
p <- colMeans(x)
z <- sweep(x, 2, p)
s <- sum(p * (1 - p))
g <- (z %*% t(z)) / s
sKin <- g / 2
diag(sKin) <- diag(sKin) + 1e-5
ids <- paste0("L", seq_len(nrow(x)))
rownames(sKin) <- colnames(sKin) <- ids
bv <- wheat.Y[, 1]
set.seed(1)
sx <- rep(c("male", "female"), length.out = nrow(x))
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(nrow(x))
cg[iM] <- 0.5
cg[iF] <- 0.5
bench("CIMMYT wheat", z, bv, sx, s, 0.12 * as.numeric(t(cg) %*% sKin %*% cg), sKin)

# --- HS mouse (real sex) ---
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
bench("HS mouse (real sex)", z, bv, sx, s, 0.15 * as.numeric(t(cg) %*% sKin %*% cg), sKin)

# --- PIC pig (52k SNP) ---
if (file.exists("/tmp/pig/FileS1/genotypes.txt")) {
  g_raw <- fread("/tmp/pig/FileS1/genotypes.txt", showProgress = FALSE)
  ids <- as.character(g_raw[[1]])
  mm <- as.matrix(g_raw[, -1])
  rm(g_raw)
  gc()
  p <- colMeans(mm) / 2
  z <- sweep(mm, 2, 2 * p)
  rm(mm)
  gc()
  s <- 2 * sum(p * (1 - p))
  g <- tcrossprod(z) / s
  sKin <- g / 2
  diag(sKin) <- diag(sKin) + 1e-5
  rm(g)
  gc()
  rownames(sKin) <- colnames(sKin) <- ids
  ebv <- fread("/tmp/pig/FileS1/ebvs.txt")
  bv <- ebv$ebv3[match(as.integer(ids), ebv$Id)]
  set.seed(1)
  sx <- sample(rep(c("male", "female"), length.out = length(ids)))
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(length(ids))
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  bench("PIC pig (52k SNP)", z, bv, sx, s, 0.15 * as.numeric(t(cg) %*% sKin %*% cg), sKin)
} else {
  cat("PIC pig: /tmp/pig/FileS1/genotypes.txt absent, skipped\n")
}
