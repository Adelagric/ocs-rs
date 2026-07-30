# §3.1 re-verified at operational caps (review point 6): does optiSel still halt
# just inside the constraint when the support is 89, not 19? Reports, for the mouse
# panel at ΔF ∈ {2,1,0.5%}, the coancestry each solver reaches as a % of the cap ub,
# the gain gap, and the support — so "support-first's edge = optiSel's unspent budget"
# can be checked at the operating point, not only at the loose one.
# Usage: Rscript research/repro/exactness_operational.R
suppressMessages({ library(BGLR); library(optiSel); library(ocsrs) })

data(mice, package = "BGLR")
x <- mice.X
gender <- as.character(mice.pheno$GENDER)
sx <- ifelse(gender %in% c("M", "Male", "male", "1"), "male", "female")
bv <- as.numeric(mice.pheno$Obesity.BMI)
ok <- !is.na(bv); x <- x[ok, ]; sx <- sx[ok]; bv <- bv[ok]
p <- colMeans(x) / 2; z <- sweep(x, 2, 2 * p); s <- 2 * sum(p * (1 - p))
sKin <- tcrossprod(z) / s / 2; diag(sKin) <- diag(sKin) + 1e-5
ids <- rownames(x); rownames(sKin) <- colnames(sKin) <- ids
n <- length(bv)
phen <- data.frame(Indiv = ids, Born = 1L, Breed = "m", BV = bv, Sex = sx,
                   isCandidate = TRUE, stringsAsFactors = FALSE)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
ne_of <- function(c) 1 / sum(c^2)

iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(n); cg[iM] <- 0.5; cg[iF] <- 0.5
ub_loose <- as.numeric(t(cg) %*% sKin %*% cg)
cu <- numeric(n); cu[sx == "male"] <- 0.5 / sum(sx == "male")
cu[sx == "female"] <- 0.5 / sum(sx == "female")
ub_tight <- as.numeric(t(cu) %*% sKin %*% cu)

cat(sprintf("mouse n=%d — exactness at operational caps\n", n))
cat(sprintf("%-6s %-5s %-24s %-24s %-10s\n", "Ne", "|S|", "ocsrs (coan, %ub)", "optiSel (coan, %ub)", "gain gap"))
for (t_ne in c(25, 50, 100)) {
  lo <- ub_tight; hi <- ub_loose; ub <- 0.5 * (lo + hi)
  for (it in 1:40) {
    ub <- 0.5 * (lo + hi)
    res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
    ne <- ne_of(res$c); if (ne > t_ne) lo <- ub else hi <- ub
    if (abs(ne - t_ne) / t_ne < 0.01) break
  }
  fit <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE)
  oc <- fit$parent$oc
  res <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = sx == "male")
  kin_rs <- as.numeric(t(res$c) %*% sKin %*% res$c)
  kin_os <- as.numeric(t(oc) %*% sKin %*% oc)
  cat(sprintf("%-6.0f %-5d (%.6f, %5.1f%%)      (%.6f, %5.1f%%)      %+.2e\n",
              ne_of(res$c), length(res$support),
              kin_rs, 100 * kin_rs / ub, kin_os, 100 * kin_os / ub,
              sum(bv * res$c) - sum(bv * oc)))
}
