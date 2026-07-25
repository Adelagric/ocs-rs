# Repeated timing: median (min–max) over REPS serial repetitions, machine idle.
# A single sub-second timing carries run-to-run jitter (GC, scheduler); for a paper
# table the median of several reps is the honest summary. Gains and supports are
# deterministic, so only the times are re-sampled. Each panel's data is built once,
# then optiSel and the matrix-free solve are each timed REPS times back to back.
#
# Usage: Rscript research/repro/r_binding_bench_reps.R [reps]

suppressMessages({
  library(BGLR)
  library(optiSel)
  library(ocsrs)
  library(data.table)
})

reps <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  5L
}

fmt <- function(v) {
  sprintf("median %.3fs  [min %.3f, max %.3f]", median(v), min(v), max(v))
}

bench <- function(label, z, bv, sx, s, ub, sKin, reps_os_max = 3L) {
  ids <- rownames(sKin)
  phen <- data.frame(
    Indiv = ids, Born = 1L, Breed = "x", BV = bv, Sex = sx,
    isCandidate = TRUE, stringsAsFactors = FALSE
  )
  cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)

  # optiSel: time once; repeat (up to reps_os_max) only when it is cheap enough that
  # repeating is affordable. Its slow cells are demonstrably stable across runs, so a
  # single sample there is honest; the run-to-run jitter lives in the sub-second ocsrs
  # timings, which are always repeated `reps` times.
  t_os <- system.time(
    opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )[["elapsed"]]
  if (t_os < 20) {
    for (r in seq_len(reps_os_max - 1L)) {
      t_os <- c(t_os, system.time(
        opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
      )[["elapsed"]])
    }
  }

  t_rs <- numeric(reps)
  gain <- NA
  supp <- NA
  for (r in seq_len(reps)) {
    t_rs[r] <- system.time(
      res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
    )[["elapsed"]]
    gain <- sum(bv * res$c)
    supp <- length(res$support)
  }

  spd <- median(t_os) / median(t_rs)
  cat(sprintf(
    "%-26s n=%-5d m=%-6d support=%d gain=%+.5f\n  optiSel (%dx) %s\n  ocsrs   (%dx) %s\n  speed-up (medians) %.0fx\n",
    label, length(bv), ncol(z), supp, gain, length(t_os), fmt(t_os), reps, fmt(t_rs), spd
  ))
  flush.console()
}

# --- synthetic ---
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
  rownames(sKin) <- colnames(sKin) <- paste0("I", seq_len(n))
  bv <- rnorm(n)
  sx <- sample(rep(c("male", "female"), length.out = n))
  bench("synthetic (structured)", z, bv, sx, s, 1.04 * mean(sKin), sKin)
}

# --- wheat ---
data(wheat)
x <- wheat.X
p <- colMeans(x)
z <- sweep(x, 2, p)
s <- sum(p * (1 - p))
sKin <- (z %*% t(z)) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- paste0("L", seq_len(nrow(x)))
bv <- wheat.Y[, 1]
set.seed(1)
sx <- rep(c("male", "female"), length.out = nrow(x))
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(nrow(x))
cg[iM] <- 0.5
cg[iF] <- 0.5
bench("CIMMYT wheat", z, bv, sx, s, 0.12 * as.numeric(t(cg) %*% sKin %*% cg), sKin)

# --- mouse ---
data(mice, package = "BGLR")
x <- mice.X
gender <- as.character(mice.pheno$GENDER)
sx <- ifelse(gender %in% c("M", "Male", "male", "1"), "male", "female")
bv <- as.numeric(mice.pheno$Obesity.BMI)
ok <- !is.na(bv)
x <- x[ok, ]
sx <- sx[ok]
bv <- bv[ok]
p <- colMeans(x) / 2
z <- sweep(x, 2, 2 * p)
s <- 2 * sum(p * (1 - p))
sKin <- tcrossprod(z) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- rownames(x)
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(length(bv))
cg[iM] <- 0.5
cg[iF] <- 0.5
bench("HS mouse (real sex)", z, bv, sx, s, 0.15 * as.numeric(t(cg) %*% sKin %*% cg), sKin)

# --- pig (built once, then timed reps) ---
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
  sKin <- tcrossprod(z) / s / 2
  diag(sKin) <- diag(sKin) + 1e-5
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
}
