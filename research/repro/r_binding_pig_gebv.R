# PIC pig sexed sub-panel with a genuine GEBV as the selection criterion.
#
# The panel ships estimated breeding values, but those come from the breeding
# programme's own evaluation and are not genomic predictions. Here the criterion is
# computed the way genomic selection actually computes it: a GBLUP fitted on the
# panel's own phenotypes and its own genomic relationship matrix, and read back as
# the posterior mean of the genetic effect. Training uses every animal with a trait-3
# phenotype (3141 of 3534); the GEBV is then available for all of them, including the
# 119 animals of the sexed sub-panel that carry no phenotype of their own and whose
# value therefore comes purely from genomic relationships — the case genomic selection
# exists for.
#
# The instance is otherwise the sexed sub-panel of r_binding_pig_sexed.R: sex read off
# the pedigree (sire => male, dam => female), 1194 animals, G recomputed on them.
#
# Usage: Rscript research/repro/r_binding_pig_gebv.R [reps]

suppressMessages({
  library(data.table)
  library(BGLR)
  library(optiSel)
  library(ocsrs)
})

set.seed(20240617) # BGLR is a Gibbs sampler; fix the chain
reps <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  as.integer(commandArgs(trailingOnly = TRUE)[1])
} else {
  5L
}
fmt <- function(v) sprintf("median %.3fs [min %.3f, max %.3f]", median(v), min(v), max(v))
dir <- "/tmp/pig/FileS1"

ped <- fread(file.path(dir, "pedigree.txt"), colClasses = "character")
sires <- unique(ped$SIRE[ped$SIRE != "0"])
dams <- unique(ped$DAM[ped$DAM != "0"])

cat("loading genotypes (934 MB)...\n")
flush.console()
g_raw <- fread(file.path(dir, "genotypes.txt"), showProgress = FALSE)
ids_all <- as.character(g_raw[[1]])
m_all <- as.matrix(g_raw[, -1])
rm(g_raw)
gc()

# --- GBLUP on the whole panel -------------------------------------------------
# G here is the training kernel over all 3534 animals, so relationships to the
# phenotyped animals inform the prediction for the unphenotyped ones.
p_all <- colMeans(m_all) / 2
z_all <- sweep(m_all, 2, 2 * p_all)
s_all <- 2 * sum(p_all * (1 - p_all))
g_train <- tcrossprod(z_all) / s_all
diag(g_train) <- diag(g_train) + 1e-5
rm(z_all)
gc()

ph <- fread(file.path(dir, "phenotypes.txt"), na.strings = ".")
y <- ph$t3[match(as.integer(ids_all), ph$ID)]
cat(sprintf("GBLUP on trait 3: %d phenotyped of %d\n", sum(!is.na(y)), length(y)))
flush.console()

fit <- BGLR(
  y = y, ETA = list(list(K = g_train, model = "RKHS")),
  nIter = 12000, burnIn = 3000, verbose = FALSE,
  saveAt = file.path(tempdir(), "bglr_")
)
gebv_all <- as.numeric(fit$ETA[[1]]$u)
h2 <- fit$ETA[[1]]$varU / (fit$ETA[[1]]$varU + fit$varE)
cat(sprintf("estimated genomic heritability h2 = %.3f\n", h2))
rm(g_train)
gc()

# --- the sexed sub-panel ------------------------------------------------------
sex_all <- ifelse(ids_all %in% sires, "male",
                  ifelse(ids_all %in% dams, "female", NA_character_))
keep <- !is.na(sex_all)
ids <- ids_all[keep]
sx <- sex_all[keep]
bv <- gebv_all[keep]

ebv <- fread(file.path(dir, "ebvs.txt"))
ebv3 <- ebv$ebv3[match(as.integer(ids), ebv$Id)]
cat(sprintf("sub-panel n=%d; correlation GEBV vs the panel's own EBV = %.3f\n",
            length(ids), cor(bv, ebv3)))
cat(sprintf("  of these, %d have no trait-3 phenotype of their own\n",
            sum(is.na(ph$t3[match(as.integer(ids), ph$ID)]))))
flush.console()

m_mat <- m_all[keep, ]
rm(m_all)
gc()
p <- colMeans(m_mat) / 2
z <- sweep(m_mat, 2, 2 * p)
rm(m_mat)
gc()
s <- 2 * sum(p * (1 - p))
sKin <- tcrossprod(z) / s / 2
diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- ids

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
    fit_oc <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  )[["elapsed"]])
}
oc <- fit_oc$parent$oc

t_rs <- numeric(reps)
for (r in seq_len(reps)) {
  t_rs[r] <- system.time(
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  )[["elapsed"]]
}

kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)

cat(sprintf("\n=== PIC pig sexed sub-panel, GEBV criterion: n=%d m=%d, ub=%.8f ===\n",
            length(ids), ncol(z), ub))
cat(sprintf("optiSel (3x) %s  gain=%.6f coancestry=%.8f (%.3f%%) support=%d\n",
            fmt(t_os), sum(bv * oc), kin_os, 100 * kin_os / ub, sum(oc > 1e-6)))
cat(sprintf("ocsrs   (%dx) %s  gain=%.6f coancestry=%.8f (%.3f%%) support=%d\n",
            reps, fmt(t_rs), sum(bv * res$c), kin_rs, 100 * kin_rs / ub,
            length(res$support)))
cat(sprintf("speed-up (medians) %.0fx   gain gap %+.6f\n",
            median(t_os) / median(t_rs), sum(bv * res$c) - sum(bv * oc)))
cat(sprintf("feasible: sum=%.10f min=%.2e males=%.10f within ub=%s\n",
            sum(res$c), min(res$c), sum(res$c[sx == "male"]),
            kin_rs <= ub * (1 + 1e-9)))
cat(sprintf("selected: %d males, %d females\n",
            sum(res$c[sx == "male"] > 1e-9), sum(res$c[sx == "female"] > 1e-9)))

n <- length(ids)
write.table(sKin, sprintf("/tmp/bench_K_gebv_%d.csv", n), sep = ",",
            row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = bv, oc = oc, sex = sx),
            sprintf("/tmp/bench_bc_gebv_%d.csv", n), sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_gebv_%d.txt", n))
cat(sprintf("exported /tmp/bench_{K,bc,ub}_gebv_%d\n", n))
