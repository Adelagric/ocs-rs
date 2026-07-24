# The two conventions that fail silently across the R/Rust boundary get explicit
# tests: R's column-major matrix layout, and 1-based indices. Both are checked by
# recomputing the answer with R's own arithmetic rather than by trusting the solver.

make_panel <- function(n = 60, m = 90, seed = 11) {
  set.seed(seed)
  p <- runif(m, 0.05, 0.5)
  geno <- matrix(rbinom(n * m, 2, rep(p, each = n)), n, m)
  z <- geno - matrix(2 * p, n, m, byrow = TRUE)
  list(Z = z, p = p, s = vanraden_scale(p), b = rnorm(n))
}

# A cap partway between the uniform plan and a single candidate's own coancestry:
# binding, but not so tight that the support grows to a large fraction of n.
binding_cap <- function(g, frac = 0.3) {
  n <- nrow(g)
  kmin <- sum(g) / n^2
  kmin + frac * (mean(diag(g)) - kmin)
}

dense_g <- function(panel, ridge = 1e-5) {
  panel$Z %*% t(panel$Z) / panel$s + ridge * diag(nrow(panel$Z))
}

test_that("the simplex optimum is feasible and on the coancestry boundary", {
  panel <- make_panel()
  g <- dense_g(panel)
  k <- binding_cap(g)

  res <- ocs_solve(panel$Z, panel$b, k = k, s = panel$s)

  expect_equal(res$status, "Solved")
  expect_equal(sum(res$c), 1, tolerance = 1e-9)
  expect_gte(min(res$c), -1e-12)
  expect_lte(res$quad, k * (1 + 1e-9))
  # the reported gain is R's own b'c, not a solver-internal number
  expect_equal(res$gain, sum(panel$b * res$c), tolerance = 1e-10)
})

test_that("Z is read column-major", {
  # The coancestry recomputed through R's own G must match the solver's. Reading Z
  # transposed would build a different G and break this to many digits.
  panel <- make_panel(seed = 3)
  g <- dense_g(panel)
  k <- binding_cap(g)

  res <- ocs_solve(panel$Z, panel$b, k = k, s = panel$s)

  expect_equal(as.numeric(t(res$c) %*% g %*% res$c), res$quad, tolerance = 1e-10)
})

test_that("the support is 1-based and lists exactly the positive contributions", {
  panel <- make_panel(seed = 5)
  k <- binding_cap(dense_g(panel))

  res <- ocs_solve(panel$Z, panel$b, k = k, s = panel$s)

  expect_true(all(res$support >= 1L))
  expect_true(all(res$support <= nrow(panel$Z)))
  expect_equal(sort(res$support), which(res$c > 1e-12))
})

test_that("the sexed form splits the budget in half by sex", {
  panel <- make_panel(seed = 7)
  n <- nrow(panel$Z)
  male <- rep(c(TRUE, FALSE), length.out = n)
  k <- binding_cap(dense_g(panel))

  res <- ocs_solve(panel$Z, panel$b, k = k, s = panel$s, male = male)

  expect_equal(res$status, "Solved")
  expect_equal(sum(res$c[male]), 0.5, tolerance = 1e-9)
  expect_equal(sum(res$c[!male]), 0.5, tolerance = 1e-9)
  expect_lte(res$quad, k * (1 + 1e-9))
})

test_that("per-candidate caps are respected", {
  panel <- make_panel(seed = 9)
  n <- nrow(panel$Z)
  k <- binding_cap(dense_g(panel))

  res <- ocs_solve(panel$Z, panel$b, k = k, s = panel$s, caps = rep(0.08, n))

  expect_equal(res$status, "Solved")
  expect_lte(max(res$c), 0.08 + 1e-9)
  expect_equal(sum(res$c), 1, tolerance = 1e-9)
})

test_that("a tighter cap never yields a higher gain", {
  panel <- make_panel(seed = 13)
  g <- dense_g(panel)
  loose <- ocs_solve(panel$Z, panel$b, k = binding_cap(g, 0.5), s = panel$s)
  tight <- ocs_solve(panel$Z, panel$b, k = binding_cap(g, 0.2), s = panel$s)

  expect_gte(loose$gain, tight$gain - 1e-9)
})

test_that("malformed input is rejected with a message, not a crash", {
  panel <- make_panel(n = 30, m = 40, seed = 17)
  expect_error(ocs_solve(as.numeric(panel$Z), panel$b, k = 0.5, s = panel$s), "matrix")
  expect_error(ocs_solve(panel$Z, panel$b[-1], k = 0.5, s = panel$s), "length")
  expect_error(ocs_solve(panel$Z, panel$b, k = 0.5, s = -1), "positive")
  expect_error(
    ocs_solve(panel$Z, panel$b, k = 0.5, s = panel$s, male = rep(TRUE, 3)),
    "length"
  )
})

test_that("a PLINK trio written from R reads back correctly", {
  # A third independent encoding of the format (after the Rust and Python test
  # suites): the decoder would have to be wrong in the same way three times.
  n <- 37L # not a multiple of 4, so every row ends in padding bits
  m <- 8L
  set.seed(23)
  dosage <- matrix(rbinom(n * m, 2, 0.3), m, n) # markers x individuals
  dosage[1, 2] <- NA # a missing call

  prefix <- file.path(tempdir(), "ocsrs_test_panel")
  on.exit(unlink(paste0(prefix, c(".bed", ".bim", ".fam"))), add = TRUE)

  code_of <- function(g) if (is.na(g)) 1L else c(3L, 2L, 0L)[g + 1L]
  raw_bytes <- as.raw(c(0x6c, 0x1b, 0x01))
  for (j in seq_len(m)) {
    row <- integer(ceiling(n / 4))
    for (i in seq_len(n)) {
      slot <- (i - 1L) %/% 4L + 1L
      row[slot] <- bitwOr(row[slot], bitwShiftL(code_of(dosage[j, i]), 2L * ((i - 1L) %% 4L)))
    }
    raw_bytes <- c(raw_bytes, as.raw(row))
  }
  writeBin(raw_bytes, paste0(prefix, ".bed"))
  writeLines(sprintf("FAM%d IND%d 0 0 1 -9", seq_len(n), seq_len(n)), paste0(prefix, ".fam"))
  writeLines(sprintf("1 rs%d 0 %d A G", seq_len(m), seq_len(m) * 1000), paste0(prefix, ".bim"))

  panel <- read_plink(prefix)

  expect_equal(dim(panel$Z), c(n, m))
  expect_equal(panel$ids[1:2], c("IND1", "IND2"))
  # frequencies and centring recomputed in R from the same calls
  freq <- apply(dosage, 1, function(row) mean(row, na.rm = TRUE) / 2)
  expect_equal(panel$p, freq, tolerance = 1e-12)
  expected_z <- t(dosage) - matrix(2 * freq, n, m, byrow = TRUE)
  expected_z[is.na(expected_z)] <- 0
  expect_equal(unname(panel$Z), expected_z, tolerance = 1e-12)
  expect_equal(panel$s, 2 * sum(freq * (1 - freq)), tolerance = 1e-12)
})

test_that("a panel read from PLINK solves", {
  n <- 40L
  m <- 60L
  set.seed(29)
  dosage <- matrix(rbinom(n * m, 2, 0.3), m, n)
  prefix <- file.path(tempdir(), "ocsrs_test_solve")
  on.exit(unlink(paste0(prefix, c(".bed", ".bim", ".fam"))), add = TRUE)

  raw_bytes <- as.raw(c(0x6c, 0x1b, 0x01))
  for (j in seq_len(m)) {
    row <- integer(ceiling(n / 4))
    for (i in seq_len(n)) {
      slot <- (i - 1L) %/% 4L + 1L
      code <- c(3L, 2L, 0L)[dosage[j, i] + 1L]
      row[slot] <- bitwOr(row[slot], bitwShiftL(code, 2L * ((i - 1L) %% 4L)))
    }
    raw_bytes <- c(raw_bytes, as.raw(row))
  }
  writeBin(raw_bytes, paste0(prefix, ".bed"))
  writeLines(sprintf("FAM%d IND%d 0 0 1 -9", seq_len(n), seq_len(n)), paste0(prefix, ".fam"))
  writeLines(sprintf("1 rs%d 0 %d A G", seq_len(m), seq_len(m) * 1000), paste0(prefix, ".bim"))

  panel <- read_plink(prefix)
  g <- panel$Z %*% t(panel$Z) / panel$s + 1e-5 * diag(n)
  res <- ocs_solve(panel$Z, rnorm(n), k = binding_cap(g), s = panel$s)

  expect_equal(res$status, "Solved")
  expect_equal(sum(res$c), 1, tolerance = 1e-9)
})
