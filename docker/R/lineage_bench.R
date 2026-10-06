#!/usr/bin/env Rscript
#
# Time the tools/bench.R scenarios against the installed hypervolume and save
# every individual timing (not just the median), so the report can show spread.
#
# The scenario definitions are taken from the harness's tools/bench.R itself:
# its top-level expressions are evaluated up to, but not including, the final
# record/compare branch. The timed code is therefore exactly what bench.R
# times; only the bookkeeping differs.
#
#   Rscript lineage_bench.R <harness root> <out.rds>
#
# BENCH_TIMES (default 20, as in bench.R) sets iterations per scenario.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript lineage_bench.R <harness root> <out.rds>")
}
harness  <- normalizePath(args[[1]], mustWork = TRUE)
out_path <- args[[2]]
times    <- as.integer(Sys.getenv("BENCH_TIMES", "20"))

suppressPackageStartupMessages(library(microbenchmark))

is_mode_branch <- function(e) {
  is.call(e) && identical(e[[1]], as.name("if")) &&
    identical(e[[2]], as.name("compare_mode"))
}

bench_env <- new.env(parent = globalenv())
# bench.R defines `root` from the working directory; point it at a scratch dir
# so nothing can be written into the (read-only) harness.
setwd(tempdir())
for (e in parse(file.path(harness, "tools", "bench.R"), keep.source = FALSE)) {
  if (is_mode_branch(e)) break
  # bench.R shells out to git for its own bookkeeping; there is no git (or
  # .git) in the image, which only produces warnings we don't need.
  suppressWarnings(eval(e, bench_env))
}
scenarios <- bench_env$scenarios
if (!is.list(scenarios) || length(scenarios) == 0L) {
  stop("no scenarios found in ", file.path(harness, "tools", "bench.R"))
}

times_ns <- list()
errors   <- list()
for (nm in names(scenarios)) {
  message("  ", nm, "...", appendLF = FALSE)
  gc(verbose = FALSE, full = TRUE)
  tryCatch({
    mb <- microbenchmark(scenarios[[nm]](), times = times)
    times_ns[[nm]] <- as.numeric(mb$time)
    message(sprintf(" %.2f ms (median of %d)", median(mb$time) / 1e6, times))
  }, error = function(e) {
    message(" ERROR: ", conditionMessage(e))
    errors[[nm]] <<- conditionMessage(e)
  })
}

saveRDS(list(
  times_ns = times_ns,
  errors = errors,
  meta = list(
    commit = Sys.getenv("HV_COMMIT", NA_character_),
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r_version = R.version.string,
    platform = R.version$platform,
    blas = sessionInfo()$BLAS,
    lapack = La_library(),
    cores = parallel::detectCores(),
    omp_threads = Sys.getenv("OMP_NUM_THREADS", NA_character_),
    times = times
  )
), out_path)
message("wrote ", out_path)
