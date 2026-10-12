# ocsrs 0.4.0

R bindings for continuous optimum contribution selection. `ocs_solve` maximises
genetic gain under a quadratic relationship bound, non-negative contributions, optional individual
caps, and either a total contribution of one or half from each sex.

## Install from source

R 4.2+ and Rust 1.85+ with Cargo are required. Clone the **full repository** and
run from its root:

```sh
R CMD INSTALL bindings/r
```

The R crate depends on the core crate through a relative Cargo path. Copying only
`bindings/r` or building a standalone R source tarball does not include that
dependency and is not a supported installation path.

## Runnable example

```r
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
```

The same example is executable from this directory with
`Rscript examples/quickstart.R`. The positive ridge contributes
`ridge * sum(c^2)` to the quadratic relationship constraint; include it when choosing a feasible
cap. A `MaxIter` result has no guarantee of being an optimum.

## Data contract

The API reports `quad = t(c) %*% G %*% c`. With the additive relationship
convention, group coancestry is `quad / 2` using `K = G / 2`; a bound `ub`
on group coancestry corresponds to `k = 2 * ub`. This score does not specify
realised offspring inbreeding without a mating allocation and reference population.

The API accepts centred genotypes **`Z`, not a relationship matrix `G`**. It uses
`G = tcrossprod(Z) / s + ridge * I`. Passing `G` as the first argument instead
would square that matrix's effect. Arbitrary relationship matrices or matrices
modified by blending, shrinkage or pedigree adjustment need a separate
adapter/factorization that reproduces the desired operator.

`b`, `male` and `caps` must be in exactly the row order of `Z`. `read_plink(prefix)`
returns dense centred `Z`, allele frequencies `p`, scale `s` and `ids` in `.fam`
row order. `ids` contains **IID only**: FID and `.fam` sex are not returned. Keep
those separately, check uniqueness before matching by IID, and align every vector
to `panel$ids`. File dimension checks do not verify that the identities or order in
the three files match; use a trio produced together.

An already-double R matrix is borrowed for a solve without copying the genotype
buffer; integer matrices are first converted to double. Both the R reader and this
API use dense matrices. Packed 2-bit PLINK solving is currently available only in
Rust through `read_packed` and `solve_sexed_packed`.
