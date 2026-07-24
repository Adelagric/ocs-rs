# support-first vs optiSel, both called from the same R session on the same instance.
#
# This is the claim the R binding exists to make testable: an optiSel user swaps one
# call and gets the same optimum. Conventions are aligned exactly rather than
# approximately — optiSel constrains c' sKin c <= ub with sKin = G/2 + 1e-5 I, so the
# identical problem for ocs_solve() is k = 2*ub on G + 2e-5 I.
#
# Usage: Rscript research/repro/r_binding_optisel.R [n] [m]

suppressMessages(library(optiSel))
suppressMessages(library(ocsrs))

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 1000L
m <- if (length(args) >= 2) as.integer(args[2]) else 500L

set.seed(1)
p <- runif(m, 0.05, 0.5)
# Structured population: families sharing a base genotype + ~5% mutations, so the
# kinship matrix has real off-diagonal structure (a non-trivial OCS instance).
n_fam <- max(2L, n %/% 20L)
base <- matrix(rbinom(n_fam * m, 2, rep(p, each = n_fam)), n_fam, m)
fam <- sample(n_fam, n, replace = TRUE)
M <- base[fam, ]
idx <- which(matrix(runif(as.numeric(n) * m), n, m) < 0.05)
M[idx] <- rbinom(length(idx), 2, p[((idx - 1L) %/% n) + 1L])
Z <- sweep(M, 2, 2 * p)
s <- 2 * sum(p * (1 - p))
G <- (Z %*% t(Z)) / s
sKin <- G / 2
diag(sKin) <- diag(sKin) + 1e-5
ids <- paste0("I", seq_len(n))
rownames(sKin) <- colnames(sKin) <- ids
phen <- data.frame(
  Indiv = ids, Born = 1L, Breed = "X", BV = rnorm(n),
  Sex = sample(c("male", "female"), n, replace = TRUE),
  isCandidate = TRUE, stringsAsFactors = FALSE
)

cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
ub <- 1.04 * cand$mean$sKin

t_os <- system.time(
  fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
)
oc <- fit$parent$oc

t_rs <- system.time(
  res <- ocs_solve(Z, phen$BV, k = 2 * ub, s = s, ridge = 2e-5,
                   male = phen$Sex == "male")
)

# Both solutions scored in optiSel's own metric, so the comparison is like for like.
kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)
gain_os <- sum(phen$BV * oc)
gain_rs <- sum(phen$BV * res$c)

cat(sprintf("instance      : n=%d m=%d, sexed, ub=%.8f\n", n, m, ub))
cat(sprintf("optiSel       : gain=%.6f  coancestry=%.8f (%.4f%% of ub)  support=%d  %.3fs\n",
            gain_os, kin_os, 100 * kin_os / ub, sum(oc > 1e-6), t_os[["elapsed"]]))
cat(sprintf("ocsrs         : gain=%.6f  coancestry=%.8f (%.4f%% of ub)  support=%d  %.3fs\n",
            gain_rs, kin_rs, 100 * kin_rs / ub, length(res$support), t_rs[["elapsed"]]))
cat(sprintf("speed-up      : %.0fx\n", t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4)))
cat(sprintf("gain gap      : %+.6f in favour of %s\n",
            gain_rs - gain_os, if (gain_rs > gain_os) "ocsrs" else "optiSel"))
cat(sprintf("sex split     : males=%.10f females=%.10f\n",
            sum(res$c[phen$Sex == "male"]), sum(res$c[phen$Sex == "female"])))
cat(sprintf("feasible      : sum=%.10f min=%.2e within ub=%s\n",
            sum(res$c), min(res$c), kin_rs <= ub * (1 + 1e-9)))
