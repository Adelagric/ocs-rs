# A second *exact* solver, independently implemented and published in the field.
#
# COMA (Endelman, Genetics 229(2):iyae193, 2025; github.com/jendelman/COMA) solves, at
# ploidy 2,   max h'y   s.t.  y'Ky <= Ft1,  sum(y) = 1,  sum(female*y) = 0.5,
#                             min <= y <= max
# through CVXR and ECOS — a strict superset of the problem support-first solves. With
# K = sKin = G/2 (the optiSel convention used throughout these scripts) the cap maps as
# k = 2*Ft1, and Ft1 = dF + (1 - dF)*Ft0 with Ft0 = mean((2*diag(K) - 1)).
#
# So: bisect our own cap to a breeder operating point (target Ne), convert that cap to
# the dF COMA takes, and compare the two optima on the same matrix. The reviews asked
# for a second exact solver; ECOS is another conic interior-point method rather than a
# different algorithm class, but the implementation, the modelling layer and the group
# are independent of ours.
#
# Usage: Rscript research/repro/coma_crosscheck.R
suppressMessages({ library(BGLR); library(COMA); library(ocsrs) })

ne_of <- function(c) 1 / sum(c^2)

panel_wheat <- function() {
  data(wheat, package = "BGLR")
  x <- wheat.X
  bv <- as.numeric(wheat.Y[, 1])
  # Autogamous: no recorded sex, an arbitrary balanced split, as in the manuscript.
  sx <- rep(c("male", "female"), length.out = nrow(x))
  list(name = "CIMMYT wheat", x = x, bv = bv, sx = sx)
}

panel_mouse <- function() {
  data(mice, package = "BGLR")
  x <- mice.X
  gender <- as.character(mice.pheno$GENDER)
  sx <- ifelse(gender %in% c("M", "Male", "male", "1"), "male", "female")
  bv <- as.numeric(mice.pheno$Obesity.BMI)
  ok <- !is.na(bv)
  list(name = "HS mouse (real sex)", x = x[ok, ], bv = bv[ok], sx = sx[ok])
}

rows <- list()
for (make in list(panel_wheat, panel_mouse)) {
  pn <- make()
  x <- pn$x; bv <- pn$bv; sx <- pn$sx
  p <- colMeans(x) / 2
  z <- sweep(x, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))
  sKin <- tcrossprod(z) / s / 2
  diag(sKin) <- diag(sKin) + 1e-5          # the ridge, in sKin units
  ids <- rownames(x); if (is.null(ids)) ids <- paste0("I", seq_len(nrow(x)))
  rownames(sKin) <- colnames(sKin) <- ids
  n <- length(bv)
  male <- sx == "male"

  # Bracket for the bisection: uniform within-sex (tight) to the gain-greedy pair (loose).
  cu <- numeric(n)
  cu[male] <- 0.5 / sum(male); cu[!male] <- 0.5 / sum(!male)
  iM <- which(male)[which.max(bv[male])]; iF <- which(!male)[which.max(bv[!male])]
  cg <- numeric(n); cg[iM] <- 0.5; cg[iF] <- 0.5
  lo <- as.numeric(t(cu) %*% sKin %*% cu); hi <- as.numeric(t(cg) %*% sKin %*% cg)

  # COMA's inbreeding base, from its own definition.
  Ft0 <- mean((2 * diag(sKin) - 1) / (2 - 1))

  cat(sprintf("\n%s  n=%d  m=%d   (Ft0 = %.6f)\n", pn$name, n, ncol(x), Ft0))
  cat(sprintf("%-5s %-9s %-11s %-11s %-9s %-9s %-11s\n",
              "Ne", "cap ub", "|S| ocs-rs", "|S| COMA", "gain rs", "gain COMA", "max|c-y|"))

  for (t_ne in c(25, 50, 100)) {
    l <- lo; h <- hi; ub <- 0.5 * (l + h)
    for (it in 1:40) {
      ub <- 0.5 * (l + h)
      r <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = male)
      ne <- ne_of(r$c)
      if (ne > t_ne) l <- ub else h <- ub
      if (abs(ne - t_ne) / t_ne < 0.01) break
    }
    t0 <- proc.time()[["elapsed"]]
    rs <- ocs_solve(z, bv, k = 2 * ub, s = s, ridge = 2e-5, male = male)
    t_rs <- proc.time()[["elapsed"]] - t0

    dF <- (ub - Ft0) / (1 - Ft0)             # invert Ft1 = dF + (1 - dF) Ft0
    parents <- data.frame(id = ids, merit = bv, min = 0, max = 1,
                          female = !male, stringsAsFactors = FALSE)
    t0 <- proc.time()[["elapsed"]]
    cm <- try(COMA::ocs(dF = dF, parents = parents, ploidy = 2, K = sKin), silent = TRUE)
    t_cm <- proc.time()[["elapsed"]] - t0
    if (inherits(cm, "try-error")) {
      cat(sprintf("%-5.0f %-9.6f  COMA failed: %s\n", t_ne, ub,
                  trimws(strsplit(as.character(cm), ":")[[1]][2])))
      next
    }
    y <- cm$oc$value[match(ids, cm$oc$id)]
    y[is.na(y)] <- 0
    gain_rs <- sum(bv * rs$c); gain_cm <- sum(bv * y)
    kin_cm <- as.numeric(t(y) %*% sKin %*% y)
    cat(sprintf("%-5.0f %-9.6f %-11d %-11d %-9.5f %-9.5f %-11.2e  (%.2fs vs %.2fs, COMA coan %.6f)\n",
                t_ne, ub, sum(rs$c > 1e-6), sum(y > 1e-6), gain_rs, gain_cm,
                max(abs(rs$c - y)), t_rs, t_cm, kin_cm))
    rows[[length(rows) + 1]] <- data.frame(
      panel = pn$name, n = n, target_ne = t_ne, ub = ub, dF = dF,
      support_ocsrs = sum(rs$c > 1e-6), support_coma = sum(y > 1e-6),
      gain_ocsrs = gain_rs, gain_coma = gain_cm,
      coan_ocsrs = as.numeric(t(rs$c) %*% sKin %*% rs$c), coan_coma = kin_cm,
      max_abs_diff = max(abs(rs$c - y)), s_ocsrs = t_rs, s_coma = t_cm,
      stringsAsFactors = FALSE)
  }
}
if (length(rows)) {
  out <- do.call(rbind, rows)
  dir.create("artifacts", showWarnings = FALSE)
  write.csv(out, "artifacts/coma_crosscheck.csv", row.names = FALSE)
  cat("\nwrote artifacts/coma_crosscheck.csv\n")
}
