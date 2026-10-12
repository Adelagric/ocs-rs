# Shared output paths and acceptance checks for automated reproduction.
repro_out <- Sys.getenv("REPRO_OUTPUT_DIR", "artifacts")
repro_inputs <- Sys.getenv("REPRO_INPUT_DIR", file.path(repro_out, "repro_inputs"))
dir.create(repro_out, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_inputs, recursive = TRUE, showWarnings = FALSE)

# Simulation RNG streams depend on thread scheduling. One thread is the default
# for repeatable seeded inputs; timings still use the native solver's own path.
repro_sim_threads <- suppressWarnings(as.integer(Sys.getenv("REPRO_SIM_THREADS", "1")))
if (length(repro_sim_threads) != 1L || is.na(repro_sim_threads) || repro_sim_threads < 1L)
  stop("REPRO_SIM_THREADS must be a positive integer")

# Some MaCS versions ignore the seed passed by runMacs. Replay the same founders
# and post-generation R RNG state instead of assuming a seeded rerun is identical.
repro_founders <- function(key, nInd, nChr, segSites, seed) {
  state_dir <- Sys.getenv("REPRO_SIM_STATE_DIR", "artifacts/simulation_state")
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(state_dir, paste0(key, ".rds"))
  params <- list(nInd = nInd, nChr = nChr, segSites = segSites, seed = seed,
                 threads = repro_sim_threads, species = "GENERIC",
                 alpha_version = as.character(packageVersion("AlphaSimR")))
  if (file.exists(path)) {
    state <- readRDS(path)
    if (!identical(state$params, params))
      stop("founder checkpoint parameters/version differ; select a separate REPRO_SIM_STATE_DIR")
  } else {
    set.seed(seed)
    pop <- runMacs(nInd = nInd, nChr = nChr, segSites = segSites,
                   species = "GENERIC", nThreads = repro_sim_threads)
    state <- list(params = params, founders = pop, rng_after = .Random.seed)
    temp <- paste0(path, ".", Sys.getpid(), ".tmp")
    saveRDS(state, temp)
    if (!file.rename(temp, path)) stop("could not save founder checkpoint")
  }
  assign(".Random.seed", state$rng_after, envir = .GlobalEnv)
  record <- data.frame(key = key, seed = seed, threads = repro_sim_threads,
                       alpha_version = params$alpha_version,
                       checkpoint = basename(path),
                       md5 = unname(tools::md5sum(path)))
  manifest <- file.path(repro_out, "simulation_state_manifest.csv")
  write.table(record, manifest, sep = ",", row.names = FALSE,
               col.names = !file.exists(manifest), append = file.exists(manifest))
  state$founders
}

checked_solve <- function(z, bv, k, s, ridge, male) {
  res <- ocs_solve(z, bv, k = k, s = s, ridge = ridge, male = male)
  if (!identical(res$status, "Solved")) stop("solver did not converge: ", res$status)
  c <- res$c
  q <- sum(as.numeric(crossprod(z, c))^2) / s + ridge * sum(c^2)
  if (any(!is.finite(c)) || min(c) < -1e-7 ||
      abs(sum(c[male]) - .5) > 1e-7 || abs(sum(c[!male]) - .5) > 1e-7 ||
      q > k + 1e-7 * max(1, abs(k))) stop("solver primal verification failed")
  res
}

require_target <- function(ne, target, tolerance = .015) {
  if (!is.finite(ne) || abs(ne - target) / target > tolerance)
    stop("cap calibration did not reach contribution Ne target: ", ne, " vs ", target)
}
