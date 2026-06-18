suppressMessages({library(BGLR); library(optiSel)})
data(mice, package = "BGLR")
X <- mice.X
ids <- rownames(X)
g <- as.character(mice.pheno$GENDER)
cat("GENDER coding:", paste(unique(g), collapse = ","), "\n")
sx <- ifelse(g %in% c("M", "Male", "male", "1"), "male", "female")
cat("sex: males", sum(sx == "male"), "females", sum(sx == "female"), "\n")

bv <- as.numeric(mice.pheno$Obesity.BMI)        # real heritable phenotype (proxy EBV)
ok <- !is.na(bv)
X <- X[ok, ]; ids <- ids[ok]; sx <- sx[ok]; bv <- bv[ok]
cat("after dropping NA phenotypes: n =", length(ids), "\n")

p <- colMeans(X) / 2; Z <- sweep(X, 2, 2 * p); s <- 2 * sum(p * (1 - p))
G <- tcrossprod(Z) / s
sKin <- G / 2; diag(sKin) <- diag(sKin) + 1e-5
rownames(sKin) <- colnames(sKin) <- ids

phen <- data.frame(Indiv = ids, Born = 1L, Breed = "mouse", BV = bv, Sex = sx,
                   isCandidate = TRUE, stringsAsFactors = FALSE)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(length(ids)); cg[iM] <- 0.5; cg[iF] <- 0.5
k_greedy <- as.numeric(t(cg) %*% sKin %*% cg)
ub <- 0.15 * k_greedy
t <- system.time(res <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE))
oc <- res$parent$oc
cat(sprintf("OPTISEL mice REAL-SEX n=%d: temps=%.2fs gain=%.5f support=%d cKc=%.5f ub=%.5f\n",
            length(ids), t["elapsed"], res$obj.fun, sum(oc > 1e-6),
            as.numeric(t(oc) %*% sKin %*% oc), ub))
n <- length(ids)
write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = bv, oc = oc, sex = sx), sprintf("/tmp/bench_bc_%d.csv", n),
            sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_%d.txt", n))
cat("EXPORT OK n=", n, "\n")
