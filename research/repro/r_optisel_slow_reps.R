# optiSel timed 3x on the two slow cells (synthetic n=5000, PIC pig), the only ones
# the reps run sampled just once. Same instances/caps as r_binding_bench_reps.R.
# ocsrs is already pinned (5x, sub-% spread); this fills the last gap in Table 2.
#
# Usage: Rscript research/repro/r_optisel_slow_reps.R

suppressMessages({
  library(BGLR)
  library(optiSel)
  library(data.table)
})

med <- function(v) sprintf("median %.3fs  [min %.3f, max %.3f]  runs %s",
                           median(v), min(v), max(v), paste(sprintf("%.2f", v), collapse = " "))

time_optisel <- function(label, sKin, bv, sx, ub, reps = 3L) {
  ids <- rownames(sKin)
  phen <- data.frame(
    Indiv = ids, Born = 1L, Breed = "x", BV = bv, Sex = sx,
    isCandidate = TRUE, stringsAsFactors = FALSE
  )
  cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
  t <- numeric(reps)
  for (r in seq_len(reps)) {
    t[r] <- system.time(
      opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
    )[["elapsed"]]
    cat(sprintf("  %s rep %d: %.3fs\n", label, r, t[r]))
    flush.console()
  }
  cat(sprintf("%s optiSel %s\n\n", label, med(t)))
}

# --- synthetic n=5000 ---
n <- 5000L
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
sKin <- tcrossprod(z) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- paste0("I", seq_len(n))
bv <- rnorm(n)
sx <- sample(rep(c("male", "female"), length.out = n))
time_optisel("synthetic n=5000", sKin, bv, sx, 1.04 * mean(sKin))

# --- PIC pig ---
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
rm(z)
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
time_optisel("PIC pig", sKin, bv, sx, 0.15 * as.numeric(t(cg) %*% sKin %*% cg))
