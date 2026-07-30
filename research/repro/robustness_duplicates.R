# Numerical robustness under near-/exact-duplicate genotypes (review point 7b).
#
# Public panels contain clones, full sibs and duplicated samples. Two identical
# genotype rows make G have two identical rows — G_S is exactly singular, and a
# closed form that inverts it would explode. The claim to test: the ridge εI floors
# the spectrum, so κ(G_S) is bounded by ~(λmax+ε)/ε rather than infinite, and the
# active set with its degeneracy handling still terminates at a feasible KKT point.
#
# Usage: Rscript research/repro/robustness_duplicates.R

suppressMessages(library(ocsrs))
set.seed(7)

n0 <- 400L
m <- 800L
p <- runif(m, 0.05, 0.5)
# Structured base: 20 families sharing a base genotype + ~5% mutations.
n_fam <- 20L
base <- matrix(rbinom(n_fam * m, 2, rep(p, each = n_fam)), n_fam, m)
fam <- sample(n_fam, n0, replace = TRUE)
M <- base[fam, ]
idx <- which(matrix(runif(as.numeric(n0) * m), n0, m) < 0.05)
M[idx] <- rbinom(length(idx), 2, p[((idx - 1L) %/% n0) + 1L])
bv <- rnorm(n0)
sx <- rep(c("male", "female"), length.out = n0)

# Inject EXACT duplicates of five high-BV individuals (clones): identical genotype,
# identical BV, identical sex — the worst degeneracy for the closed form.
top <- order(bv, decreasing = TRUE)[1:5]
M <- rbind(M, M[top, ])
bv <- c(bv, bv[top])
sx <- c(sx, sx[top])
n <- nrow(M)
clone_pairs <- cbind(top, (n0 + 1L):(n0 + 5L)) # (original, clone) index pairs

p2 <- colMeans(M) / 2
Z <- sweep(M, 2, 2 * p2)
s <- 2 * sum(p2 * (1 - p2))
ne_of <- function(c) 1 / sum(c^2)

# G with and without ridge, to show the ridge is what bounds the conditioning.
G0 <- tcrossprod(Z) / s
G <- G0
diag(G) <- diag(G) + 1e-5

# Operating cap: bisect to Ne≈25.
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(n)
cg[iM] <- 0.5
cg[iF] <- 0.5
k_loose <- as.numeric(t(cg) %*% G %*% cg)
cu <- numeric(n)
cu[sx == "male"] <- 0.5 / sum(sx == "male")
cu[sx == "female"] <- 0.5 / sum(sx == "female")
k_tight <- as.numeric(t(cu) %*% G %*% cu)
lo <- k_tight
hi <- k_loose
for (it in 1:30) {
  k <- 0.5 * (lo + hi)
  res <- ocs_solve(Z, bv, k = k, s = s, ridge = 1e-5, male = sx == "male")
  if (ne_of(res$c) > 25) lo <- k else hi <- k
}
res <- ocs_solve(Z, bv, k = k, s = s, ridge = 1e-5, male = sx == "male")

S <- res$support
GS <- G[S, S]
evS <- eigen(GS, symmetric = TRUE, only.values = TRUE)$values
GS0 <- G0[S, S]
evS0 <- eigen(GS0, symmetric = TRUE, only.values = TRUE)$values
coan <- as.numeric(t(res$c) %*% G0 %*% res$c)

cat(sprintf("panel: n=%d (incl. 5 exact clones) m=%d\n", n, m))
cat(sprintf("solve: status=%s  |S|=%d  Ne=%.1f  feasible(coan<=k)=%s\n",
            res$status, length(S), ne_of(res$c), coan <= k * (1 + 1e-9)))
cat(sprintf("min eigenvalue of G_S WITHOUT ridge : %.3e  (≈0 ⇒ duplicates make it singular)\n", min(evS0)))
cat(sprintf("min eigenvalue of G_S WITH ridge    : %.3e  (floored at ε=1e-5)\n", min(evS)))
cat(sprintf("condition number κ(G_S) with ridge  : %.3e  (bounded ~λmax/ε, not infinite)\n",
            max(evS) / min(evS)))

# How were the clone pairs treated? Identical candidates should get a combined
# contribution that is stable (the split between them is immaterial and feasible).
cat("clone pairs (original,clone) -> contributions:\n")
for (r in seq_len(nrow(clone_pairs))) {
  a <- clone_pairs[r, 1]
  b <- clone_pairs[r, 2]
  cat(sprintf("  pair %d: c[orig]=%.6f  c[clone]=%.6f  sum=%.6f  (both in support: %s)\n",
              r, res$c[a], res$c[b], res$c[a] + res$c[b],
              (a %in% S) && (b %in% S)))
}
