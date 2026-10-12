# Small repeated-seed smoke test for the simulation's actual inputs and solution.
suppressMessages({library(AlphaSimR); library(ocsrs)})
source("research/repro/repro_common.R")
state_dir <- tempfile("simulation-state-")
Sys.setenv(REPRO_SIM_STATE_DIR = state_dir)

simulate_once <- function() {
  set.seed(1)
  threads <- if (exists("repro_sim_threads")) repro_sim_threads else NULL
  founders <- repro_founders("smoke_founders", nInd = 60, nChr = 4,
                              segSites = 120, seed = 1)
  sp <- SimParam$new(founders)
  if (!is.null(threads)) sp$nThreads <- threads
  sp$setSexes("yes_rand")
  sp$addTraitA(nQtlPerChr = 40)
  sp$addSnpChip(nSnpPerChr = 40)
  parents <- newPop(founders, simParam = sp)
  candidates <- randCross(parents, nCrosses = 250, nProgeny = 1, simParam = sp)
  M <- pullSnpGeno(candidates, simParam = sp)
  bv <- as.numeric(gv(candidates)[, 1])
  male <- candidates@sex == "M"
  p <- colMeans(M) / 2
  z <- sweep(M, 2, 2 * p)
  s <- 2 * sum(p * (1 - p))
  ridge <- 2e-5
  cu <- ifelse(male, .5 / sum(male), .5 / sum(!male))
  cg <- numeric(nrow(M))
  cg[which(male)[which.max(bv[male])]] <- .5
  cg[which(!male)[which.max(bv[!male])]] <- .5
  q <- function(c) sum(as.numeric(crossprod(z, c))^2) / s + ridge * sum(c^2)
  lo <- q(cu); hi <- q(cg)
  for (i in 1:40) {
    k <- (lo + hi) / 2
    fit <- checked_solve(z, bv, k, s, ridge, male)
    ne <- 1 / sum(fit$c^2)
    if (ne > 50) lo <- k else hi <- k
    if (abs(ne - 50) / 50 < .01) break
  }
  require_target(ne, 50)
  list(founders = founders@geno, genotypes = M, sex = male, criterion = bv, support = fit$support,
       gain = sum(bv * fit$c), c = fit$c, threads = sp$nThreads)
}

first <- simulate_once()
second <- simulate_once()
if (!identical(first, second)) message("differing fields: ", paste(names(first)[!mapply(identical, first, second)], collapse = ", "))
stopifnot(identical(first, second))
cat("checkpoint replay identical: genotype/sex/criterion/support/gain/contributions; threads=",
    first$threads, "\n", sep = "")
unlink(state_dir, recursive = TRUE)
