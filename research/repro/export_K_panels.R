# Export the dense kinship K (= sKin), bv and sex for wheat and mouse, so the
# given-G active-set prototype (sf_operational_ne.py) can be timed at operational Ne
# on the same instances. Same conventions as r_binding_real_panels.R.
# Usage: Rscript research/repro/export_K_panels.R
suppressMessages(library(BGLR))

export_one <- function(n, sKin, bv, sx) {
  write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",",
              row.names = FALSE, col.names = FALSE)
  write.table(data.frame(bv = bv, oc = 0, sex = sx), sprintf("/tmp/bench_bc_%d.csv", n),
              sep = ",", row.names = FALSE)
  cat(sprintf("exported n=%d (m in the GRM only)\n", n))
}

# wheat
data(wheat)
x <- wheat.X
p <- colMeans(x)
z <- sweep(x, 2, p)
s <- sum(p * (1 - p))
sKin <- (z %*% t(z)) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
set.seed(1)
sx <- rep(c("male", "female"), length.out = nrow(x))
export_one(nrow(x), sKin, wheat.Y[, 1], sx)

# mouse
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
export_one(nrow(x), sKin, bv, sx)
