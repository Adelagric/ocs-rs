suppressMessages(library(optiSel))
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1]); m <- as.integer(args[2])
set.seed(1)
p <- runif(m, 0.05, 0.5)
# Structured population: families sharing a base genotype + ~5% mutations,
# so the kinship matrix has real off-diagonal structure (non-trivial OCS).
n_fam <- max(2L, n %/% 20L)
base  <- matrix(rbinom(n_fam * m, 2, rep(p, each = n_fam)), n_fam, m)
fam   <- sample(n_fam, n, replace = TRUE)
M     <- base[fam, ]
idx   <- which(matrix(runif(as.numeric(n) * m), n, m) < 0.05)
M[idx] <- rbinom(length(idx), 2, p[((idx - 1L) %/% n) + 1L])
Z <- sweep(M, 2, 2 * p)
s <- 2 * sum(p * (1 - p))
G <- (Z %*% t(Z)) / s            # VanRaden relationship (diag ~1)
sKin <- G / 2                    # kinship (diag ~0.5), optiSel convention
diag(sKin) <- diag(sKin) + 1e-5  # ridge
ids <- paste0("I", seq_len(n))
rownames(sKin) <- colnames(sKin) <- ids
phen <- data.frame(Indiv = ids, Born = 1L, Breed = "X", BV = rnorm(n),
                   Sex = sample(c("male", "female"), n, replace = TRUE),
                   isCandidate = TRUE, stringsAsFactors = FALSE)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
ub <- 1.04 * cand$mean$sKin
t <- system.time(res <- opticont("max.BV", cand, list(ub.sKin = ub),
                                 trace = FALSE, quiet = TRUE))
oc <- res$parent$oc
cKc <- as.numeric(t(oc) %*% sKin %*% oc)
cat(sprintf("OPTISEL n=%d m=%d : temps=%.2fs gain=%.5f support=%d cKc=%.5f ub=%.5f\n",
            n, m, t["elapsed"], res$obj.fun, sum(oc > 1e-6), cKc, ub))
write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = phen$BV, oc = oc, sex = phen$Sex),
            sprintf("/tmp/bench_bc_%d.csv", n), sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_%d.txt", n))
cat("EXPORT OK\n")
