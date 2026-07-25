#' Exact, matrix-free optimum contribution selection
#'
#' Maximises genetic gain `t(b) %*% c` subject to a cap on the mean coancestry of
#' the next generation, `t(c) %*% G %*% c <= k`, with non-negative contributions
#' summing to one — or, when `male` is given, to one half per sex.
#'
#' The relationship matrix `G = Z Z' / s + ridge * I` is never formed: kinship
#' products come straight from `Z`, and only the handful of candidates in the
#' optimal support enter the linear algebra. The result is the exact optimum for
#' continuous contributions under a single coancestry constraint (KKT-certified),
#' not a heuristic.
#'
#' @param Z Numeric matrix of centred genotypes, `n` candidates by `m` markers
#'   (dosages minus twice the allele frequency). [read_plink()] returns one.
#' @param b Numeric vector of length `n`: the breeding values to maximise.
#' @param k Coancestry cap; the constraint is `t(c) %*% G %*% c <= k`.
#' @param s VanRaden scale `2 * sum(p * (1 - p))` — see [vanraden_scale()].
#' @param ridge Small ridge making `G` positive definite (default `1e-5`).
#' @param male Optional logical vector of length `n`. When supplied the sexed
#'   problem is solved, each sex contributing one half.
#' @param caps Optional numeric vector of length `n` of per-candidate upper
#'   bounds `c <= caps`.
#' @param max_iter,tol Active-set iteration cap and reduced-cost tolerance.
#'
#' @return An object of class `ocs_result`: a list with the contributions `c`,
#'   the active `support` (1-based indices), the `gain` `t(b) %*% c`, the realised
#'   coancestry `quad`, counts `iterations` and `products`, and `status`
#'   (`"Solved"` or `"MaxIter"`). Only a `"Solved"` result is an optimum.
#'
#' @examples
#' \dontrun{
#' panel <- read_plink("cattle")
#' res <- ocs_solve(panel$Z, b = ebv, k = 0.03, s = panel$s)
#' res$support
#' }
#' @export
ocs_solve <- function(Z, b, k, s, ridge = 1e-5, male = NULL, caps = NULL,
                      max_iter = 10000L, tol = 1e-9) {
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix of centred genotypes (n candidates x m markers).")
  }
  n <- nrow(Z)
  m <- ncol(Z)
  b <- as.numeric(b)
  if (length(b) != n) {
    stop(sprintf("`b` has length %d, expected nrow(Z) = %d.", length(b), n))
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) {
    stop("`k` must be a single finite number (the coancestry cap).")
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0) {
    stop("`s` must be a single positive number; see `vanraden_scale()`.")
  }

  # R stores matrices column-major and the Rust side borrows that buffer directly, so
  # neither `as.numeric` (which would flatten-copy) nor a transpose is needed. Coerce
  # to double only when Z is not already double, so the common case copies nothing.
  if (!is.double(Z)) storage.mode(Z) <- "double"
  n_i <- as.integer(n)
  m_i <- as.integer(m)
  it <- as.integer(max_iter)
  s <- as.numeric(s)
  ridge <- as.numeric(ridge)
  k <- as.numeric(k)
  tol <- as.numeric(tol)

  if (!is.null(male)) {
    male <- as.logical(male)
    if (length(male) != n) {
      stop(sprintf("`male` has length %d, expected nrow(Z) = %d.", length(male), n))
    }
    if (anyNA(male)) stop("`male` must not contain NA.")
    male <- as.integer(male)
  }
  if (!is.null(caps)) {
    caps <- as.numeric(caps)
    if (length(caps) != n) {
      stop(sprintf("`caps` has length %d, expected nrow(Z) = %d.", length(caps), n))
    }
    if (any(caps < 0)) stop("`caps` must be non-negative.")
  }

  out <- if (is.null(male) && is.null(caps)) {
    sf_solve(Z, n_i, m_i, s, ridge, b, k, it, tol)
  } else if (is.null(male)) {
    sf_solve_capped(Z, n_i, m_i, s, ridge, b, caps, k, it, tol)
  } else if (is.null(caps)) {
    sf_solve_sexed(Z, n_i, m_i, s, ridge, b, male, k, it, tol)
  } else {
    sf_solve_sexed_capped(Z, n_i, m_i, s, ridge, b, male, caps, k, it, tol)
  }

  structure(out, class = "ocs_result")
}

#' @export
print.ocs_result <- function(x, ...) {
  cat(sprintf(
    "<ocs_result> status=%s  |support|=%d  gain=%.6g  coancestry=%.6g  (%d iterations, %d products)\n",
    x$status, length(x$support), x$gain, x$quad, x$iterations, x$products
  ))
  invisible(x)
}

#' VanRaden scaling from allele frequencies
#'
#' @param p Numeric vector of allele frequencies.
#' @return `2 * sum(p * (1 - p))`, the scale in `G = Z Z' / s`.
#' @export
vanraden_scale <- function(p) {
  p <- as.numeric(p)
  2 * sum(p * (1 - p))
}

#' Read a PLINK 1 binary panel
#'
#' Reads `<prefix>.bed`, `<prefix>.bim` and `<prefix>.fam` into the centred
#' genotypes [ocs_solve()] takes. Genotypes are centred by twice the allele
#' frequency estimated from the non-missing calls, and missing calls are imputed
#' to the marker mean. Only SNP-major `.bed` files are read (PLINK >= 1.9 writes
#' those by default).
#'
#' @param prefix Path without extension, e.g. `"data/cattle"`.
#' @return An object of class `ocs_panel`: a list with the centred matrix `Z`,
#'   the allele frequencies `p`, the VanRaden scale `s`, and the individual `ids`.
#' @examples
#' \dontrun{
#' panel <- read_plink("cattle")
#' res <- ocs_solve(panel$Z, b = ebv, k = 0.03, s = panel$s)
#' }
#' @export
read_plink <- function(prefix) {
  out <- sf_read_plink(path.expand(as.character(prefix)))
  # The Rust side already emits column-major, so this just attaches dimensions.
  z <- matrix(out$z, nrow = out$n, ncol = out$m)
  rownames(z) <- out$ids
  structure(
    list(Z = z, p = out$p, s = out$s, ids = out$ids),
    class = "ocs_panel"
  )
}

#' @export
print.ocs_panel <- function(x, ...) {
  cat(sprintf(
    "<ocs_panel> %d individuals x %d markers  VanRaden s=%.6g\n",
    nrow(x$Z), ncol(x$Z), x$s
  ))
  invisible(x)
}
