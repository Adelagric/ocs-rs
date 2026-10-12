# Run with Rscript after installing ocsrs from the full repository checkout.
library(ocsrs)

set.seed(11)
n <- 12L
m <- 20L
frequencies <- runif(m, 0.1, 0.5)
dosages <- matrix(rbinom(n * m, 2, rep(frequencies, each = n)), n, m)
p <- colMeans(dosages) / 2
Z <- sweep(dosages, 2, 2 * p)
s <- vanraden_scale(p)
b <- rnorm(n)
male <- seq_len(n) %% 2L == 1L
caps <- rep(0.15, n)
ridge <- 1e-5 # G = tcrossprod(Z) / s + ridge * I

# A feasible plan contributes half from each sex, uniformly within that sex.
c0 <- ifelse(male, 0.5 / sum(male), 0.5 / sum(!male))
stopifnot(all(c0 <= caps))
q0 <- sum(crossprod(Z, c0)^2) / s + ridge * sum(c0^2)
k <- q0 + 0.05 # This cap admits c0; it is not a universal breeding threshold.

res <- ocs_solve(Z, b, k = k, s = s, ridge = ridge, male = male, caps = caps)
if (res$status != "Solved") stop("No solved optimum: ", res$status)
print(res)
print(res$support) # 1-based R indices
print(res$c[res$support])
