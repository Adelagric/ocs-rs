suppressMessages({library(BGLR); library(optiSel)})
data(wheat)
X <- wheat.X                      # 599 x 1279, DArT 0/1
p <- colMeans(X)
Z <- sweep(X, 2, p)               # centred
denom <- sum(p * (1 - p))
G <- (Z %*% t(Z)) / denom         # marker-based relationship (diag ~1)
sKin <- G / 2                     # kinship
diag(sKin) <- diag(sKin) + 1e-5
n <- nrow(X)
ids <- paste0("L", seq_len(n))
rownames(sKin) <- colnames(sKin) <- ids
bv <- wheat.Y[, 1]                # phenotype env1 as breeding-value proxy
cat(sprintf("wheat: n=%d markers=%d | sKin diag~%.3f offdiag range[%.3f,%.3f] | bv range[%.2f,%.2f]\n",
            n, ncol(X), mean(diag(sKin)), min(sKin - diag(diag(sKin))), max(sKin - diag(diag(sKin))),
            min(bv), max(bv)))

# optiSel requires two sexes. The wheat lines have none (autogamous), so we
# impose an ARBITRARY 2-group partition only to satisfy that requirement. This
# benchmarks the two SOLVERS on an identical, well-defined problem over a REAL
# wheat relationship matrix and REAL phenotypes — it is not a statement about
# wheat biology.
set.seed(1)
sx <- rep(c("male", "female"), length.out = n)
phen <- data.frame(Indiv = ids, Born = 1L, Breed = "wheat", BV = bv,
                   Sex = sx, isCandidate = TRUE, stringsAsFactors = FALSE)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
# Centred VanRaden GRM: the uniform plan has ~0 kinship, so bound relative to
# the gain-greedy plan (all weight on the best male + best female) instead.
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(n); cg[iM] <- 0.5; cg[iF] <- 0.5
k_greedy <- as.numeric(t(cg) %*% sKin %*% cg)
ub <- 0.12 * k_greedy   # force spreading -> active constraint, modest support
cat(sprintf("k_greedy=%.5f  ub=%.5f\n", k_greedy, ub))
tt <- system.time(res <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE))
oc <- res$parent$oc
cat(sprintf("optiSel (real wheat GRM, n=%d): temps=%.2fs gain=%.5f support=%d cKc=%.5f ub=%.5f\n",
            n, tt["elapsed"], res$obj.fun, sum(oc > 1e-6), as.numeric(t(oc) %*% sKin %*% oc), ub))
write.table(sKin, "/tmp/bench_K_599.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = bv, oc = oc, sex = sx), "/tmp/bench_bc_599.csv", sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), "/tmp/bench_ub_599.txt")
cat("EXPORT OK\n")
