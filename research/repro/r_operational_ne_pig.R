# Operational-point re-benchmark on the PIC pig sexed sub-panel (real EBV + real sex),
# the high-marker acid test: n=1194, m=52843. Here the dense-G memory wall does NOT
# bite (n is small), so optiSel runs fine and this measures pure SPEED at a large
# marker count — the worst case for the matrix-free path, which re-derives kinship
# from Z (O(n m) per iteration) while optiSel reuses a G formed once.
#
# Ne_c = 1/sum(c^2); ΔF ≈ 1/(2 Ne_c). Conventions as elsewhere: optiSel bounds
# c' sKin c <= ub on sKin = G/2 + 1e-5 I, so ocs_solve uses k = 2*ub, ridge = 2e-5.
#
# Needs /tmp/pig/FileS1/{genotypes,ebvs,pedigree}.txt.
# Usage: Rscript research/repro/r_operational_ne_pig.R

suppressMessages({
  library(data.table)
  library(optiSel)
  library(ocsrs)
})

dir <- "/tmp/pig/FileS1"
ped <- fread(file.path(dir, "pedigree.txt"), colClasses = "character")
sires <- unique(ped$SIRE[ped$SIRE != "0"])
dams <- unique(ped$DAM[ped$DAM != "0"])
stopifnot(length(intersect(sires, dams)) == 0)

cat("loading genotypes (934 MB)...\n")
flush.console()
g_raw <- fread(file.path(dir, "genotypes.txt"), showProgress = FALSE)
ids_all <- as.character(g_raw[[1]])
sex_all <- ifelse(ids_all %in% sires, "male", ifelse(ids_all %in% dams, "female", NA_character_))
keep <- !is.na(sex_all)
ids <- ids_all[keep]
sx <- sex_all[keep]
m_mat <- as.matrix(g_raw[keep, -1])
rm(g_raw)
gc()

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
stopifnot(!anyNA(bv))
cat(sprintf("pig sexed sub-panel: n=%d m=%d (real EBV trait3, real sex)\n", length(ids), ncol(z)))

phen <- data.frame(
  Indiv = ids, Born = 1L, Breed = "pig", BV = bv, Sex = sx,
  isCandidate = TRUE, stringsAsFactors = FALSE
)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)

n <- length(ids)
ne_of <- function(c) 1 / sum(c^2)
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(n)
cg[iM] <- 0.5
cg[iF] <- 0.5
ub_loose <- as.numeric(t(cg) %*% sKin %*% cg)
cu <- numeric(n)
cu[sx == "male"] <- 0.5 / sum(sx == "male")
cu[sx == "female"] <- 0.5 / sum(sx == "female")
ub_tight <- as.numeric(t(cu) %*% sKin %*% cu)

rows <- list()
for (t_ne in c(25, 50, 100)) {
  lo <- ub_tight
  hi <- ub_loose
  chosen <- 0.5 * (lo + hi)
  for (it in 1:35) {
    ub <- 0.5 * (lo + hi)
    chosen <- ub
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
    ne <- ne_of(res$c)
    if (ne > t_ne) lo <- ub else hi <- ub
    if (abs(ne - t_ne) / t_ne < 0.015) break
  }
  ub <- chosen
  # Capture results, then time separately (assignments inside replicate() would not
  # escape its evaluation frame).
  fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  oc <- fit$parent$oc
  t_os <- median(replicate(3, system.time(
    opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )[["elapsed"]]))
  res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  t_rs <- median(replicate(3, system.time(
    ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  )[["elapsed"]]))
  ne <- ne_of(res$c)
  gap <- sum(bv * res$c) - sum(bv * oc)
  kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)
  kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
  su <- t_os / max(t_rs, 1e-4)
  cat(sprintf(
    "Ne=%3.0f dF=%.2f%% | |S|=%3d (optiSel %3d, using %.0f%% of ub) | optiSel %6.3fs  ocsrs %6.3fs  speedup %5.2fx | gaingap %+.2e\n",
    ne, 50 / ne, length(res$support), sum(oc > 1e-6), 100 * kin_os / ub, t_os, t_rs, su, gap
  ))
  rows[[length(rows) + 1]] <- data.frame(
    panel = "pig (real EBV+sex)", n = n, m = ncol(z), target_ne = t_ne, achieved_ne = ne,
    delta_f_pct = 50 / ne, ub = ub, support_ocsrs = length(res$support),
    support_optisel = sum(oc > 1e-6), optisel_pct_of_ub = 100 * kin_os / ub,
    optisel_s = t_os, ocsrs_s = t_rs, speedup = su, gain_gap = gap
  )
}
out <- do.call(rbind, rows)
dir.create("artifacts", showWarnings = FALSE)
write.csv(out, "artifacts/operational_ne_pig.csv", row.names = FALSE)
cat("\n[operational_ne_pig] wrote artifacts/operational_ne_pig.csv\n")
