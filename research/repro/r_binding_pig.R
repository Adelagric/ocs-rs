# The matrix-free solver vs optiSel on the PIC pig panel (n=3534, 52,843 SNP).
#
# The severest m-axis test: with 52k markers the matrix-free products and the Gram
# cache both scale with m, and the R binding must copy a 1.5 GB genotype matrix — so
# this is where the shipped end-to-end number is most exposed. Instance and cap
# reproduce pig_export.R exactly (ub = 0.15 * k_greedy on sKin = G/2 + 1e-5 I), so
# the identical problem for ocs_solve is k = 2*ub on G + 2e-5 I.
#
# Needs /tmp/pig/FileS1/{genotypes.txt,ebvs.txt}. Usage: Rscript r_binding_pig.R

suppressMessages({
  library(data.table)
  library(optiSel)
  library(ocsrs)
})

setwd("/tmp/pig/FileS1")
cat("loading genotypes (934 MB)...\n")
flush.console()
g_raw <- fread("genotypes.txt", showProgress = FALSE)
ids <- as.character(g_raw[[1]])
m_mat <- as.matrix(g_raw[, -1])
rm(g_raw)
gc()
cat(sprintf("genotypes: %d x %d\n", nrow(m_mat), ncol(m_mat)))
flush.console()

p <- colMeans(m_mat) / 2
z <- sweep(m_mat, 2, 2 * p)
rm(m_mat)
gc()
s <- 2 * sum(p * (1 - p))
g <- tcrossprod(z) / s
sKin <- g / 2
diag(sKin) <- diag(sKin) + 1e-5
rm(g)
gc()
rownames(sKin) <- colnames(sKin) <- ids

ebv <- fread("ebvs.txt")
bv <- ebv$ebv3[match(as.integer(ids), ebv$Id)]
set.seed(1)
sx <- sample(rep(c("male", "female"), length.out = length(ids)))

phen <- data.frame(
  Indiv = ids, Born = 1L, Breed = "pig", BV = bv, Sex = sx,
  isCandidate = TRUE, stringsAsFactors = FALSE
)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(length(ids))
cg[iM] <- 0.5
cg[iF] <- 0.5
ub <- 0.15 * as.numeric(t(cg) %*% sKin %*% cg)

cat("opticont...\n")
flush.console()
t_os <- system.time(
  fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
)
oc <- fit$parent$oc

t_marshal <- system.time(
  ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male", max_iter = 1)
)
t_rs <- system.time(
  res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
)

kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)

cat(sprintf("\n=== PIC pig: n=%d m=%d, sexed, ub=%.8f ===\n", length(ids), ncol(z), ub))
cat(sprintf("optiSel   : gain=%.6f  coancestry=%.8f (%.3f%%)  support=%d  %.3fs\n",
            sum(bv * oc), kin_os, 100 * kin_os / ub, sum(oc > 1e-6), t_os[["elapsed"]]))
cat(sprintf("ocsrs     : gain=%.6f  coancestry=%.8f (%.3f%%)  support=%d  %.3fs (copy %.3fs)\n",
            sum(bv * res$c), kin_rs, 100 * kin_rs / ub, length(res$support),
            t_rs[["elapsed"]], t_marshal[["elapsed"]]))
cat(sprintf("            iterations=%d  products=%d  status=%s\n",
            res$iterations, res$products, res$status))
cat(sprintf("speed-up  : %.0fx   gain gap: %+.6f   sex: M=%.6f F=%.6f\n",
            t_os[["elapsed"]] / max(t_rs[["elapsed"]], 1e-4),
            sum(bv * res$c) - sum(bv * oc),
            sum(res$c[sx == "male"]), sum(res$c[sx == "female"])))
cat(sprintf("feasible  : sum=%.8f  min=%.2e  within ub=%s\n",
            sum(res$c), min(res$c), kin_rs <= ub * (1 + 1e-9)))
