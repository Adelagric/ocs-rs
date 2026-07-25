# PIC pig, sexed sub-panel: the one instance with BOTH a real breeding value and a
# real recorded sex.
#
# The PIC download ships a pedigree that the other export scripts ignore. A pedigree
# assigns sex for free: any animal appearing in the SIRE column is male, in the DAM
# column female. Intersecting that with the genotyped animals recovers a real sex for
# 1194 of the 3534 (390 sires, 804 dams; no animal appears as both), and every one of
# them carries a real estimated breeding value from the accompanying genetic
# evaluation (ebvs.txt, with accuracies). So this sub-panel needs neither of the two
# stand-ins the other panels need — no phenotype standing in for a breeding value, no
# arbitrary balanced sex split.
#
# Honest caveat, to be carried into the text: these 1194 are exactly the animals that
# became parents, i.e. a post-selection subset rather than a random sample of
# selection candidates. That is a property of how sex is recoverable here, not a
# choice, and it does not affect what the benchmark measures (two solvers on one
# identical instance) — but it does mean the contribution vector is illustrative.
#
# Conventions match the other comparisons: optiSel constrains c' sKin c <= ub on
# sKin = G/2 + 1e-5 I, so ocs_solve is handed k = 2*ub on G + 2e-5 I. Cap is
# 0.15 * k_greedy, as in pig_export.R.
#
# Needs /tmp/pig/FileS1/{genotypes,ebvs,pedigree}.txt.
# Usage: Rscript research/repro/r_binding_pig_sexed.R [reps]

suppressMessages({
  library(data.table)
  library(optiSel)
  library(ocsrs)
})

reps <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  5L
}
fmt <- function(v) sprintf("median %.3fs [min %.3f, max %.3f]", median(v), min(v), max(v))

dir <- "/tmp/pig/FileS1"
cat("deriving sex from the pedigree...\n")
flush.console()
ped <- fread(file.path(dir, "pedigree.txt"), colClasses = "character")
sires <- unique(ped$SIRE[ped$SIRE != "0"])
dams <- unique(ped$DAM[ped$DAM != "0"])
stopifnot(length(intersect(sires, dams)) == 0) # no animal is both

cat("loading genotypes (934 MB)...\n")
flush.console()
g_raw <- fread(file.path(dir, "genotypes.txt"), showProgress = FALSE)
ids_all <- as.character(g_raw[[1]])
sex_all <- ifelse(ids_all %in% sires, "male",
                  ifelse(ids_all %in% dams, "female", NA_character_))
keep <- !is.na(sex_all)
cat(sprintf("genotyped %d; real sex recovered for %d (%d male, %d female)\n",
            length(ids_all), sum(keep), sum(sex_all == "male", na.rm = TRUE),
            sum(sex_all == "female", na.rm = TRUE)))
flush.console()

ids <- ids_all[keep]
sx <- sex_all[keep]
m_mat <- as.matrix(g_raw[keep, -1])
rm(g_raw)
gc()

# Allele frequencies are recomputed on the retained animals, so G describes this
# sub-panel rather than being sliced out of the full-panel matrix.
p <- colMeans(m_mat) / 2
z <- sweep(m_mat, 2, 2 * p)
rm(m_mat)
gc()
s <- 2 * sum(p * (1 - p))
sKin <- tcrossprod(z) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- ids

ebv <- fread(file.path(dir, "ebvs.txt"))
bv <- ebv$ebv3[match(as.integer(ids), ebv$Id)]
acc <- ebv$acc3[match(as.integer(ids), ebv$Id)]
stopifnot(!anyNA(bv))
cat(sprintf("real EBV (trait 3) present for all %d; mean accuracy %.3f\n", length(bv), mean(acc)))
flush.console()

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

t_os <- numeric(0)
for (r in seq_len(3L)) {
  t_os <- c(t_os, system.time(
    fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )[["elapsed"]])
}
oc <- fit$parent$oc

t_rs <- numeric(reps)
for (r in seq_len(reps)) {
  t_rs[r] <- system.time(
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  )[["elapsed"]]
}

kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)

cat(sprintf("\n=== PIC pig, sexed sub-panel: n=%d m=%d, ub=%.8f ===\n",
            length(ids), ncol(z), ub))
cat(sprintf("optiSel (3x) %s  gain=%.6f coancestry=%.8f (%.3f%%) support=%d\n",
            fmt(t_os), sum(bv * oc), kin_os, 100 * kin_os / ub, sum(oc > 1e-6)))
cat(sprintf("ocsrs   (%dx) %s  gain=%.6f coancestry=%.8f (%.3f%%) support=%d\n",
            reps, fmt(t_rs), sum(bv * res$c), kin_rs, 100 * kin_rs / ub,
            length(res$support)))
cat(sprintf("speed-up (medians) %.0fx   gain gap %+.6f\n",
            median(t_os) / median(t_rs), sum(bv * res$c) - sum(bv * oc)))
cat(sprintf("feasible: sum=%.10f min=%.2e males=%.10f females=%.10f within ub=%s\n",
            sum(res$c), min(res$c), sum(res$c[sx == "male"]),
            sum(res$c[sx == "female"]), kin_rs <= ub * (1 + 1e-9)))
cat(sprintf("selected: %d males, %d females\n",
            sum(res$c[sx == "male"] > 1e-9), sum(res$c[sx == "female"] > 1e-9)))

# Export the dense kinship at this exact cap so the dense-G prototype (sf_at_ub.py)
# can be timed on the identical instance, giving this row both Table 2 columns.
n <- length(ids)
write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = bv, oc = oc, sex = sx), sprintf("/tmp/bench_bc_%d.csv", n),
            sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_%d.txt", n))
cat(sprintf("exported /tmp/bench_{K,bc,ub}_%d for the algorithm column\n", n))
