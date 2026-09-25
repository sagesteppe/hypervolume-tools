#!/usr/bin/env Rscript
#
# Benchmark harness for the performance work in hypervolume_speedup_notes.md.
# Companion to tools/make_golden.R: goldens catch correctness drift, bench
# catches (or confirms) the actual speedup.
#
# Usage:
#   Rscript tools/bench.R              # record baseline (refuses if exists)
#   Rscript tools/bench.R --force      # overwrite baseline
#   Rscript tools/bench.R --compare    # run against installed pkg, diff vs baseline
#
# Workflow:
#   1. On the untouched baseline: `Rscript tools/bench.R` writes
#      tools/bench-baseline.rds.
#   2. After a perf commit: `R CMD INSTALL .` then
#      `Rscript tools/bench.R --compare` shows a per-scenario ratio.
#
# Uses microbenchmark rather than bench so no new package dependency is
# added -- microbenchmark ships in most R installs and is already present
# on this machine.

suppressPackageStartupMessages({
  library(hypervolume)
  library(microbenchmark)
})

args <- commandArgs(trailingOnly = TRUE)
force_overwrite <- "--force" %in% args
compare_mode    <- "--compare" %in% args

# tiny helper so we don't need rlang
`%||%` <- function(a, b) if (is.null(a)) b else a

root <- normalizePath(".", mustWork = TRUE)
out_path <- file.path(root, "tools", "bench-baseline.rds")

# --- scenarios ---------------------------------------------------------------
#
# One scenario per proposal, sized so the total harness runs in ~1-2 minutes
# on the baseline. Numbers chosen large enough that a ~10x speedup is
# visible against microbenchmark's noise, small enough to be tolerable at
# baseline O(n^2) cost.

set.seed(20260925L)

scenarios <- list()

# proposal #1: kd-tree build (n scaling)
scenarios$kdtree_build_2k <- local({
  d <- matrix(rnorm(2000L * 3L), ncol = 3L)
  function() invisible(hypervolume:::kdtree_build(d, verbose = FALSE))
})
scenarios$kdtree_build_10k <- local({
  d <- matrix(rnorm(10000L * 3L), ncol = 3L)
  function() invisible(hypervolume:::kdtree_build(d, verbose = FALSE))
})

# proposals #5 / #6: ball-query throughput
scenarios$evalfspherical_10k_500q <- local({
  d <- matrix(rnorm(10000L * 3L), ncol = 3L)
  q <- matrix(rnorm(500L * 3L), ncol = 3L)
  colnames(d) <- colnames(q) <- c("X1", "X2", "X3")
  function() invisible(hypervolume:::evalfspherical(d, 0.75, q, verbose = FALSE))
})

# proposal #2: calculate_density (dmvnorm loop -> BLAS)
scenarios$calculate_density_200x1000 <- local({
  means <- matrix(rnorm(200L * 3L), ncol = 3L)
  s     <- matrix(rnorm(1000L * 3L), ncol = 3L)
  colnames(means) <- colnames(s) <- c("X1", "X2", "X3")
  bw <- c(0.5, 0.5, 0.5)
  w  <- rep(1 / nrow(means), nrow(means))
  function() invisible(hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w, verbose = FALSE
  ))
})

# proposal #3: threshold sweep (the 1000-iteration hot path)
scenarios$hypervolume_threshold_1000 <- local({
  n <- 5000L
  pts <- matrix(runif(n * 2L), ncol = 2L)
  colnames(pts) <- c("X1", "X2")
  vals <- runif(n)
  hv <- methods::new(
    "Hypervolume",
    Name = "bench", Method = "bench fixture",
    Data = pts, Dimensionality = 2,
    Volume = 1, PointDensity = n, Parameters = list(),
    RandomPoints = pts, ValueAtRandomPoints = vals
  )
  function() invisible(suppressMessages(hypervolume_threshold(
    hv, num.thresholds = 1000L, plot = FALSE, verbose = FALSE
  )))
})

# proposal #4: estimate_probability
scenarios$estimate_probability_100q <- local({
  n <- 1000L
  pts <- matrix(runif(n * 2L), ncol = 2L)
  colnames(pts) <- c("X1", "X2")
  hv <- methods::new(
    "Hypervolume",
    Name = "bench", Method = "bench fixture",
    Data = pts, Dimensionality = 2,
    Volume = 1, PointDensity = n, Parameters = list(),
    RandomPoints = pts, ValueAtRandomPoints = runif(n)
  )
  q <- matrix(runif(100L * 2L), ncol = 2L)
  colnames(q) <- c("X1", "X2")
  function() invisible(suppressMessages(hypervolume_estimate_probability(
    hv, q, verbose = FALSE
  )))
})

# proposal #9: fastPdist2
scenarios$fastPdist2_500x500 <- local({
  A <- matrix(rnorm(500L * 5L), ncol = 5L)
  B <- matrix(rnorm(500L * 5L), ncol = 5L)
  function() invisible(hypervolume:::fastPdist2(A, B))
})

# Flagship end-to-end: hypervolume_gaussian on a Quercus subset. Uses
# only one species and a slice of samples.per.point so the baseline runs
# in a few seconds -- the wall-clock here is what a downstream user
# actually sees.
scenarios$hypervolume_gaussian_quercus <- local({
  data(quercus, package = "hypervolume", envir = environment())
  q <- as.matrix(quercus[quercus$Species == "Quercus rubra", 2:3])
  q <- q[!duplicated(q), ]
  q <- q[seq_len(min(nrow(q), 300L)), ]  # cap for bench tractability
  function() {
    set.seed(20260925L)
    invisible(suppressMessages(hypervolume_gaussian(
      q, samples.per.point = 20, verbose = FALSE
    )))
  }
})

# --- benchmarking loop -------------------------------------------------------

run_bench <- function() {
  times_ms <- list()
  for (nm in names(scenarios)) {
    message("  ", nm, "...", appendLF = FALSE)
    gc(verbose = FALSE, full = TRUE)
    # 5 iterations is enough for scenarios that take >10ms; microbenchmark's
    # setup overhead dominates for anything faster, but nothing here is that
    # fast on the baseline.
    mb <- microbenchmark(scenarios[[nm]](), times = 5L, unit = "ms")
    med_ms <- median(mb$time) / 1e6
    times_ms[[nm]] <- med_ms
    message(sprintf(" %.1f ms (median of 5)", med_ms))
  }
  times_ms
}

git_sha <- tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
git_branch <- tryCatch(
  system2("git", c("rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

# --- compare mode ------------------------------------------------------------

if (compare_mode) {
  if (!file.exists(out_path)) {
    stop(
      "No baseline at ", out_path, "\n",
      "Record one first with: Rscript tools/bench.R"
    )
  }
  baseline <- readRDS(out_path)

  message("Running current build (installed hypervolume ",
          utils::packageVersion("hypervolume"), ") ...")
  current <- run_bench()

  nms <- union(names(baseline$times_ms), names(current))
  base_v <- vapply(nms, function(n) baseline$times_ms[[n]] %||% NA_real_, numeric(1))
  curr_v <- vapply(nms, function(n) current[[n]] %||% NA_real_, numeric(1))
  ratio  <- curr_v / base_v
  pct    <- (curr_v - base_v) / base_v * 100

  tbl <- data.frame(
    scenario  = nms,
    baseline  = sprintf("%.1f ms", base_v),
    current   = sprintf("%.1f ms", curr_v),
    ratio     = sprintf("%.2fx", ratio),
    change    = sprintf("%+.0f%%", pct),
    stringsAsFactors = FALSE
  )
  message("")
  message("Baseline recorded ", baseline$meta$date,
          " on ", baseline$meta$git_branch, "@", baseline$meta$git_sha)
  message("Current build:    ", git_branch, "@", git_sha)
  message("")
  print(tbl, row.names = FALSE, right = FALSE)
  invisible()
} else {
  # --- record mode ---------------------------------------------------------
  if (file.exists(out_path) && !force_overwrite) {
    stop(
      "Baseline already exists at ", out_path, "\n",
      "Re-recording silently destroys the speedup measurement. ",
      "Pass --force only if you deliberately want to re-baseline."
    )
  }

  if (!identical(git_branch, "main") && !identical(git_branch, "perf/baseline-tests")) {
    message(
      "NOTE: current branch is '", git_branch,
      "'. Bench baselines should be recorded on the pre-perf baseline."
    )
  }

  message("Recording bench baseline...")
  times_ms <- run_bench()

  baseline <- list(
    times_ms = times_ms,
    meta = list(
      date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      package_version = as.character(utils::packageVersion("hypervolume")),
      r_version = R.version.string,
      platform = R.version$platform,
      blas = tryCatch(sessionInfo()$BLAS, error = function(e) NA_character_),
      git_sha = git_sha,
      git_branch = git_branch
    )
  )
  saveRDS(baseline, out_path)
  message("")
  message("wrote ", out_path)
  message("  package ", baseline$meta$package_version,
          "  |  ", baseline$meta$git_branch, "@", baseline$meta$git_sha)
  message("")
  message("Compare against this baseline after each perf commit with:")
  message("  R CMD INSTALL . && Rscript tools/bench.R --compare")
}
