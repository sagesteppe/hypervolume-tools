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
#   2. After a perf commit: force a clean rebuild with
#        rm -f src/*.o src/*.so && R CMD INSTALL .
#      then `Rscript tools/bench.R --compare` shows a per-scenario ratio.
#
# The clean rebuild is not optional -- R's default Makefile does not track
# header dependencies, so editing src/KDTree.h will not by itself trigger
# a rebuild of the .cpp files that include it, and --compare would silently
# read stale numbers.
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

# --report="label" appends a section to tools/bench-report.md and a set of
# rows to tools/bench-report.csv, tagged with the label. Both artifacts are
# gitignored -- they accumulate across proposals and get pasted into the PR
# body at the end.
report_label <- {
  hit <- grepl("^--report=", args)
  if (any(hit)) sub("^--report=", "", args[hit][1]) else NA_character_
}

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

# proposal #1: kd-tree build. Sizes chosen to span the actual range of
# ecological field data hypervolume() sees in the wild -- roughly the
# 5th percentile (n=500), median (~n=5000), and 99th percentile (n=10000)
# of user datasets. Anything larger is algorithmically interesting but
# doesn't reflect what published users actually run.
scenarios$kdtree_build_500 <- local({
  d <- matrix(rnorm(500L * 3L), ncol = 3L)
  function() invisible(hypervolume:::kdtree_build(d, verbose = FALSE))
})
scenarios$kdtree_build_5000 <- local({
  d <- matrix(rnorm(5000L * 3L), ncol = 3L)
  function() invisible(hypervolume:::kdtree_build(d, verbose = FALSE))
})
scenarios$kdtree_build_10000 <- local({
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

# proposal #8: hypervolume_set on two 2D overlapping boxes at a realistic
# niche-analysis scale (~5000 RandomPoints per hypervolume, ~0.5 Jaccard
# overlap). The hot path here used to be unique.data.frame after
# rbind()ing two 5000-row data.frames.
scenarios$hypervolume_set_5000 <- local({
  set.seed(20260925L)
  make_box <- function(lo, hi, seed, name) {
    set.seed(seed)
    n <- 5000L
    pts <- cbind(runif(n, lo[1], hi[1]), runif(n, lo[2], hi[2]))
    colnames(pts) <- c("X1", "X2")
    methods::new(
      "Hypervolume",
      Name = name, Method = "bench fixture",
      Data = pts, Dimensionality = 2L,
      Volume = prod(hi - lo), PointDensity = n / prod(hi - lo),
      Parameters = list(),
      RandomPoints = pts, ValueAtRandomPoints = rep(1, n)
    )
  }
  hv1 <- make_box(c(0, 0), c(1, 1), seed = 1L, name = "A")
  hv2 <- make_box(c(0.5, 0), c(1.5, 1), seed = 2L, name = "B")
  function() {
    invisible(suppressMessages(hypervolume_set(
      hv1, hv2, check.memory = FALSE, verbose = FALSE
    )))
  }
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
    # 20 iterations because the cheap scenarios (kdtree_build_500 at
    # sub-millisecond, fastPdist2 at 700us) need enough samples for the
    # median to be stable to a few percent. Cost is ~30s total on the
    # baseline, dominated by the expensive scenarios.
    mb <- microbenchmark(scenarios[[nm]](), times = 20L, unit = "ms")
    med_ms <- median(mb$time) / 1e6
    times_ms[[nm]] <- med_ms
    message(sprintf(" %.2f ms (median of 20)", med_ms))
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

  if (!is.na(report_label)) {
    md_path  <- file.path(root, "tools", "bench-report.md")
    csv_path <- file.path(root, "tools", "bench-report.csv")

    md <- c(
      sprintf("## %s", report_label),
      "",
      sprintf("Recorded %s. Current build: `%s@%s`. Baseline: `%s@%s`.",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
              git_branch, git_sha,
              baseline$meta$git_branch, baseline$meta$git_sha),
      "",
      "| scenario | baseline | current | ratio | change |",
      "|---|---|---|---|---|",
      sprintf("| %s | %s | %s | %s | %s |",
              tbl$scenario, tbl$baseline, tbl$current, tbl$ratio, tbl$change),
      ""
    )
    if (!file.exists(md_path)) {
      writeLines(c(
        "# hypervolume perf report",
        "",
        sprintf("Deltas against the untouched 3.1.7 baseline at `%s@%s`, ",
                baseline$meta$git_branch, baseline$meta$git_sha),
        "recorded by `tools/bench.R --compare --report=\"label\"`.",
        "Median of 20 microbenchmark runs.",
        ""
      ), md_path)
    }
    cat(md, sep = "\n", file = md_path, append = TRUE)
    cat("\n", file = md_path, append = TRUE)

    csv_new <- !file.exists(csv_path)
    csv_rows <- data.frame(
      label      = report_label,
      timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      git_sha    = git_sha,
      scenario   = nms,
      baseline_ms = base_v,
      current_ms  = curr_v,
      ratio      = ratio,
      pct_change = pct,
      stringsAsFactors = FALSE
    )
    utils::write.table(
      csv_rows, csv_path,
      sep = ",", row.names = FALSE, col.names = csv_new,
      append = !csv_new, quote = TRUE
    )
    message("")
    message("Appended report section '", report_label, "' to ", md_path)
    message("Appended ", nrow(csv_rows), " rows to ", csv_path)
  }

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
