#!/usr/bin/env Rscript
#
# Run a testthat directory against the *installed* hypervolume and write one
# row per test_that() block. Never exits non-zero for test failures -- a red
# test at an early commit is a data point in the lineage, not an error.
#
#   Rscript lineage_tests.R <tests/testthat dir> <out.csv> <suite label>

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) != 3L) {
  stop("usage: Rscript lineage_tests.R <tests/testthat dir> <out.csv> <suite label>")
}
.test_dir <- .args[[1]]
.out_path <- .args[[2]]
.suite    <- .args[[3]]
.commit   <- Sys.getenv("HV_COMMIT", NA_character_)

suppressPackageStartupMessages(library(testthat))

# testthat takes its edition from the installed package's DESCRIPTION, and
# 3.1.7 doesn't declare one (edition 2: expect_equal() is all.equal()). The
# tip suite was written for edition 3, so pin it there -- otherwise commits
# before the one that adds `Config/testthat/edition: 3` fail tip-suite tests
# for reasons that have nothing to do with their code. The `own` suite keeps
# whatever edition its commit declared.
if (.suite == "tip") Sys.setenv(TESTTHAT_EDITION = "3")

.cols <- c("commit", "suite", "file", "test", "status",
           "expectations", "failed", "skipped", "error", "warning", "seconds")

.rows <- tryCatch({
  res <- testthat::test_dir(
    .test_dir,
    package = "hypervolume",
    load_package = "installed",
    reporter = "summary",
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
  df <- as.data.frame(res)
  status <- ifelse(df$error, "error",
            ifelse(df$failed > 0L, "fail",
            ifelse(df$skipped, "skip", "pass")))
  data.frame(
    commit = .commit, suite = .suite,
    file = df$file, test = df$test, status = status,
    expectations = df$nb, failed = df$failed, skipped = df$skipped,
    error = df$error, warning = df$warning, seconds = round(df$real, 3),
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  # e.g. a helper file that doesn't load against this commit's package
  message("suite could not run: ", conditionMessage(e))
  data.frame(
    commit = .commit, suite = .suite, file = NA_character_,
    test = paste("SUITE ERROR:", conditionMessage(e)), status = "error",
    expectations = 0L, failed = 0L, skipped = FALSE, error = TRUE,
    warning = 0L, seconds = NA_real_, stringsAsFactors = FALSE
  )
})

utils::write.csv(.rows[, .cols], .out_path, row.names = FALSE)
message(sprintf("%s suite @ %s: %d pass, %d fail, %d error, %d skip -> %s",
                .suite, .commit,
                sum(.rows$status == "pass"), sum(.rows$status == "fail"),
                sum(.rows$status == "error"), sum(.rows$status == "skip"),
                .out_path))
