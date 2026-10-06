#!/usr/bin/env Rscript
#
# Evaluate the golden scenarios (tests/testthat/helper-golden-scenarios.R from
# the harness checkout) against the installed hypervolume and save the raw
# values. lineage_report.R then compares them bit for bit across commits.
#
# Same scenarios as tools/make_golden.R, but writes to an arbitrary path and
# never compares against the committed fixture -- that fixture was recorded
# on a different machine, and the point here is commit-vs-commit on one.
#
#   Rscript lineage_golden.R <harness root> <out.rds>

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) != 2L) {
  stop("usage: Rscript lineage_golden.R <harness root> <out.rds>")
}
.out_path <- .args[[2]]

suppressPackageStartupMessages({
  library(testthat)
  library(hypervolume)
})

# The helpers resolve paths relative to tests/testthat, as in make_golden.R.
# Variables here are dot-prefixed because the helpers are sourced into the
# global environment.
setwd(file.path(normalizePath(.args[[1]], mustWork = TRUE), "tests", "testthat"))
source("helper-hv.R", local = FALSE)
source("helper-golden-scenarios.R", local = FALSE)

.scenarios <- golden_scenarios()
.values <- list()
.errors <- list()

for (.nm in names(.scenarios)) {
  message("  golden scenario: ", .nm)
  tryCatch(
    .values[[.nm]] <- .scenarios[[.nm]](),
    error = function(e) {
      message("    ERROR: ", conditionMessage(e))
      .errors[[.nm]] <<- conditionMessage(e)
    }
  )
}

saveRDS(list(
  values = .values,
  errors = .errors,
  meta = list(
    commit = Sys.getenv("HV_COMMIT", NA_character_),
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r_version = R.version.string,
    platform = R.version$platform,
    blas = sessionInfo()$BLAS,
    lapack = La_library(),
    seed = HV_SEED
  )
), .out_path)
message("wrote ", .out_path)
