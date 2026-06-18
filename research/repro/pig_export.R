suppressMessages({library(data.table); library(optiSel)})
setwd("/tmp/pig/FileS1")
cat("loading genotypes (897MB)...\n"); flush.console()
G_raw <- fread("genotypes.txt", showProgress = FALSE)
ids <- as.character(G_raw[[1]])
M <- as.matrix(G_raw[, -1]); rm(G_raw); gc()
cat("genotypes:", nrow(M), "x", ncol(M), "\n"); flush.console()
p <- colMeans(M) / 2
Z <- sweep(M, 2, 2 * p); rm(M); gc()
s <- 2 * sum(p * (1 - p))
cat("computing VanRaden GRM (3534 x 3534)...\n"); flush.console()
G <- tcrossprod(Z) / s; rm(Z); gc()
sKin <- G / 2; diag(sKin) <- diag(sKin) + 1e-5; rm(G); gc()
rownames(sKin) <- colnames(sKin) <- ids

ebv <- fread("ebvs.txt")
bv <- ebv$ebv3[match(as.integer(ids), ebv$Id)]   # trait 3 (h2=0.38), real EBV
set.seed(1); sx <- sample(rep(c("male", "female"), length.out = length(ids)))  # arbitrary 50/50

phen <- data.frame(Indiv = ids, Born = 1L, Breed = "pig", BV = bv, Sex = sx,
                   isCandidate = TRUE, stringsAsFactors = FALSE)
cand <- candes(phen = phen, sKin = sKin, quiet = TRUE)
iM <- which(sx == "male")[which.max(bv[sx == "male"])]
iF <- which(sx == "female")[which.max(bv[sx == "female"])]
cg <- numeric(length(ids)); cg[iM] <- 0.5; cg[iF] <- 0.5
k_greedy <- as.numeric(t(cg) %*% sKin %*% cg)
ub <- 0.15 * k_greedy
cat(sprintf("opticont (ub=%.5f)...\n", ub)); flush.console()
t <- system.time(res <- opticont("max.BV", cand, list(ub.sKin = ub), trace = FALSE, quiet = TRUE))
oc <- res$parent$oc
cat(sprintf("OPTISEL pig REAL n=%d: temps=%.2fs gain=%.5f support=%d cKc=%.5f ub=%.5f\n",
            length(ids), t["elapsed"], res$obj.fun, sum(oc > 1e-6),
            as.numeric(t(oc) %*% sKin %*% oc), ub))
n <- length(ids)
write.table(sKin, sprintf("/tmp/bench_K_%d.csv", n), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(data.frame(bv = bv, oc = oc, sex = sx), sprintf("/tmp/bench_bc_%d.csv", n),
            sep = ",", row.names = FALSE)
writeLines(format(ub, digits = 15), sprintf("/tmp/bench_ub_%d.txt", n))
cat("EXPORT OK n=", n, "\n")
