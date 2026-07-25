# Export the dense kinship K, bv/sex, and cap ub for each panel, at the *exact*
# instances and caps benchmarked in r_binding_bench_all.R — so the NumPy prototype
# (sf_at_ub.py, the active set given a dense G) can be timed on the same problems the
# shipped matrix-free solver was, for a matched two-column comparison. optiSel is not
# run here (its contribution vector is not needed: sf_at_ub.py reads only bv and sex),
# so the export is fast.
#
# Writes /tmp/bench_{K,bc,ub}_{n}.{csv,txt}. Usage: Rscript research/repro/r_export_K.R

suppressMessages({
  library(BGLR)
  library(data.table)
})

emit <- function(n, sKin, bv, sx, ub) {
  write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",",
              row.names = FALSE, col.names = FALSE)
  # column 2 (oc) is unused by sf_at_ub.py; write zeros.
  write.table(data.frame(bv = bv, oc = 0, sex = sx), sprintf("/tmp/bench_bc_%d.csv", n),
              sep = ",", row.names = FALSE)
  writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_%d.txt", n))
  cat(sprintf("exported n=%d m via ub=%.8f\n", n, ub))
  flush.console()
}

# --- synthetic structured populations (identical construction to r_binding_bench_all.R) ---
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
  bv <- rnorm(n)
  sx <- sample(rep(c("male", "female"), length.out = n))
  emit(n, sKin, bv, sx, 1.04 * mean(sKin))
}

# --- CIMMYT wheat ---
data(wheat)
x <- wheat.X
p <- colMeans(x)
z <- sweep(x, 2, p)
s <- sum(p * (1 - p))
sKin <- (z %*% t(z)) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
bv <- wheat.Y[, 1]
set.seed(1)
sx <- rep(c("male", "female"), length.out = nrow(x))
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(nrow(x))
cg[iM] <- 0.5
cg[iF] <- 0.5
emit(nrow(x), sKin, bv, sx, 0.12 * as.numeric(t(cg) %*% sKin %*% cg))

# --- HS mouse (real sex) ---
data(mice, package = "BGLR")
x <- mice.X
ids <- rownames(x)
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
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(length(bv))
cg[iM] <- 0.5
cg[iF] <- 0.5
emit(length(bv), sKin, bv, sx, 0.15 * as.numeric(t(cg) %*% sKin %*% cg))

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
  sKin <- tcrossprod(z) / s / 2
  diag(sKin) <- diag(sKin) + 1e-5
  rm(z)
  gc()
  ebv <- fread("/tmp/pig/FileS1/ebvs.txt")
  bv <- ebv$ebv3[match(as.integer(ids), ebv$Id)]
  set.seed(1)
  sx <- sample(rep(c("male", "female"), length.out = length(ids)))
  iM <- which(sx == "male")[which.max(bv[sx == "male"])]
  iF <- which(sx == "female")[which.max(bv[sx == "female"])]
  cg <- numeric(length(ids))
  cg[iM] <- 0.5
  cg[iF] <- 0.5
  emit(length(ids), sKin, bv, sx, 0.15 * as.numeric(t(cg) %*% sKin %*% cg))
} else {
  cat("PIC pig: genotypes absent, skipped\n")
}
