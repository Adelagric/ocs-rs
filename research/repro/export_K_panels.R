# Export the exact constructed R benchmark instance under a named panel key.
# This file is sourced by the operational scripts; direct execution validates
# their exports rather than rebuilding a second potentially different instance.
export_panel <- function(panel, K, bv, sx, m, ridge = 1e-5) {
  stopifnot(nrow(K) == length(bv), length(sx) == length(bv), all(is.finite(K)),
            all(is.finite(bv)), all(sx %in% c("male", "female")))
  kp <- file.path(repro_inputs, paste0(panel, "_K.csv"))
  bp <- file.path(repro_inputs, paste0(panel, "_bc.csv"))
  write.table(K, kp, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(bv = bv, oc = 0, sex = sx), bp, sep = ",", row.names = FALSE)
  hashes <- c(K_md5 = unname(tools::md5sum(kp)), bc_md5 = unname(tools::md5sum(bp)))
  write.csv(data.frame(panel = panel, n = length(bv), m = m,
                      kinship_definition = "ZZt/(2s)+ridge*I", kinship_ridge = ridge,
                      K_md5 = hashes[["K_md5"]], bc_md5 = hashes[["bc_md5"]]),
            file.path(repro_inputs, paste0(panel, "_metadata.csv")), row.names = FALSE)
  hashes
}

if (sys.nframe() == 0L) {
  source("research/repro/repro_common.R")
  panels <- commandArgs(trailingOnly = TRUE)
  if (!length(panels)) panels <- c("wheat", "mouse")
  for (panel in panels) {
    meta <- read.csv(file.path(repro_inputs, paste0(panel, "_metadata.csv")))
    for (kind in c("K", "bc")) {
      path <- file.path(repro_inputs, paste0(panel, "_", kind, ".csv"))
      stopifnot(unname(tools::md5sum(path)) == meta[[paste0(kind, "_md5")]])
    }
    cat("validated exact export:", panel, "\n")
  }
}
