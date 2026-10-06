#!/usr/bin/env Rscript
#
# Collect everything docker/lineage.sh wrote into <out dir> and produce
#
#   lineage-report.md   speed, numerical identity and test status per commit
#   bench.csv           one row per commit x scenario
#   golden.csv          one row per commit x scenario x output element
#   tests.csv           one row per commit x suite x test
#
#   Rscript lineage_report.R <out dir>
#
# Needs only base R; runs in the deps image.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript lineage_report.R <out dir>")
out_dir <- normalizePath(args[[1]], mustWork = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

commits <- utils::read.delim(file.path(out_dir, "commits.tsv"),
                             colClasses = "character", quote = "")
commits$label <- sprintf("%s %s", commits$idx, commits$sha)

result_path <- function(kind, i, ext) {
  file.path(out_dir, kind, sprintf("%s-%s.%s", commits$idx[i], commits$sha[i], ext))
}
read_if <- function(path) if (file.exists(path)) readRDS(path) else NULL

bench  <- lapply(seq_len(nrow(commits)), function(i) read_if(result_path("bench", i, "rds")))
golden <- lapply(seq_len(nrow(commits)), function(i) read_if(result_path("golden", i, "rds")))

md <- character()
add <- function(...) md <<- c(md, ...)
fmt_ms <- function(x) ifelse(is.na(x), "–",
                      ifelse(x >= 100, sprintf("%.0f", x),
                      ifelse(x >= 1, sprintf("%.1f", x), sprintf("%.2f", x))))
md_table <- function(df) {
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
subject_short <- function(s, n = 60L) ifelse(nchar(s) > n, paste0(substr(s, 1L, n - 1L), "…"), s)

# --- header ------------------------------------------------------------------

any_meta <- Filter(Negate(is.null), c(bench, golden))
meta <- if (length(any_meta)) any_meta[[1]]$meta else list()
run_info <- if (file.exists(file.path(out_dir, "run-info.txt"))) {
  readLines(file.path(out_dir, "run-info.txt"))
} else character()

add("# hypervolume perf lineage", "",
    "Every commit below was installed into its own Docker image built on one",
    "shared base (same R, BLAS/LAPACK, compiler and dependency versions), then",
    "run against the same test/bench harness. Differences between rows come",
    "from the commit alone.", "",
    sprintf("- R: %s (%s)", meta$r_version %||% "?", meta$platform %||% "?"),
    sprintf("- BLAS: `%s`", meta$blas %||% "?"),
    sprintf("- LAPACK: `%s`", meta$lapack %||% "?"))
if (!is.null(meta$cores)) {
  add(sprintf("- cores visible: %s, OMP_NUM_THREADS=%s, %s iterations per bench scenario",
              meta$cores, meta$omp_threads, meta$times))
}
add(paste0("- ", run_info), "")

# --- commits -----------------------------------------------------------------

add("## Commits", "",
    "`code` marks commits that change what gets installed (R/, src/,",
    "DESCRIPTION, NAMESPACE). Only those, plus the baseline, are timed and",
    "have their numbers captured; the others install a byte-identical package.", "")
add(md_table(data.frame(
  `#` = commits$idx, commit = paste0("`", commits$sha, "`"),
  code = ifelse(commits$code_changed == "baseline", "baseline",
         ifelse(commits$code_changed == "yes", "code", "")),
  subject = subject_short(commits$subject),
  check.names = FALSE
)), "")

# --- speed -------------------------------------------------------------------

have_bench <- which(!vapply(bench, is.null, logical(1)))
bench_rows <- list()
if (length(have_bench)) {
  scen <- unique(unlist(lapply(bench[have_bench], function(b) names(b$times_ns))))
  base_i <- have_bench[[1]]
  med <- function(i, s) {
    t <- bench[[i]]$times_ns[[s]]
    if (is.null(t)) NA_real_ else median(t) / 1e6
  }
  prev_i <- NA_integer_
  for (i in have_bench) {
    for (s in scen) {
      t <- bench[[i]]$times_ns[[s]]
      m <- med(i, s)
      bench_rows[[length(bench_rows) + 1L]] <- data.frame(
        idx = commits$idx[i], commit = commits$sha[i], scenario = s,
        median_ms = m,
        q25_ms = if (is.null(t)) NA_real_ else unname(quantile(t, 0.25)) / 1e6,
        q75_ms = if (is.null(t)) NA_real_ else unname(quantile(t, 0.75)) / 1e6,
        speedup_vs_baseline = med(base_i, s) / m,
        speedup_vs_previous = if (is.na(prev_i)) NA_real_ else med(prev_i, s) / m,
        error = bench[[i]]$errors[[s]] %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }
    prev_i <- i
  }
  bench_df <- do.call(rbind, bench_rows)
  utils::write.csv(bench_df, file.path(out_dir, "bench.csv"), row.names = FALSE)

  add("## Speed", "",
      "Median ms per call, with speed-up against the baseline (first row) in",
      "brackets. Per-iteration spread is in `bench.csv` (q25/q75).", "")
  wide <- data.frame(commit = paste0("`", commits$sha[have_bench], "`"), check.names = FALSE)
  for (s in scen) {
    sub <- bench_df[bench_df$scenario == s, ]
    wide[[s]] <- ifelse(is.na(sub$median_ms), "error",
                 ifelse(sub$commit == commits$sha[base_i], fmt_ms(sub$median_ms),
                        sprintf("%s (%.1f×)", fmt_ms(sub$median_ms), sub$speedup_vs_baseline)))
  }
  # scenarios as rows reads better than 10+ columns of scenario names
  tw <- as.data.frame(t(wide[, -1L, drop = FALSE]), stringsAsFactors = FALSE)
  names(tw) <- wide$commit
  tw <- cbind(scenario = rownames(tw), tw)
  add(md_table(tw), "")
}

# --- numerical identity ------------------------------------------------------

compare_values <- function(a, b) {
  # a, b: one element of a golden scenario's output
  if (is.null(a) || is.null(b)) return(list(identical = NA, max_abs = NA_real_, max_rel = NA_real_, note = "missing"))
  if (identical(a, b)) return(list(identical = TRUE, max_abs = 0, max_rel = 0, note = ""))
  a <- unclass(a); b <- unclass(b)
  if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b)) {
    return(list(identical = FALSE, max_abs = NA_real_, max_rel = NA_real_,
                note = sprintf("shape changed (%d -> %d values)", length(b), length(a))))
  }
  d <- abs(as.numeric(a) - as.numeric(b))
  scale <- pmax(abs(as.numeric(b)), .Machine$double.xmin)
  list(identical = FALSE, max_abs = max(d, na.rm = TRUE),
       max_rel = max(d / scale, na.rm = TRUE),
       note = if (all(d == 0, na.rm = TRUE)) "same values, different attributes" else "")
}

have_golden <- which(!vapply(golden, is.null, logical(1)))
if (length(have_golden)) {
  base_i <- have_golden[[1]]
  scen <- unique(unlist(lapply(golden[have_golden], function(g) c(names(g$values), names(g$errors)))))
  rows <- list()
  prev_i <- NA_integer_
  for (i in have_golden) {
    for (s in scen) {
      cur <- golden[[i]]$values[[s]]
      els <- unique(c(names(cur), names(golden[[base_i]]$values[[s]])))
      if (!is.null(golden[[i]]$errors[[s]])) els <- "(scenario)"
      for (el in els) {
        vs_base <- compare_values(cur[[el]], golden[[base_i]]$values[[s]][[el]])
        vs_prev <- if (is.na(prev_i)) list(identical = NA, max_abs = NA_real_, max_rel = NA_real_, note = "")
                   else compare_values(cur[[el]], golden[[prev_i]]$values[[s]][[el]])
        err <- golden[[i]]$errors[[s]]
        notes <- c(if (!is.null(err)) paste("error:", err), vs_prev$note, vs_base$note)
        rows[[length(rows) + 1L]] <- data.frame(
          idx = commits$idx[i], commit = commits$sha[i], scenario = s, element = el,
          identical_to_baseline = vs_base$identical,
          max_abs_diff_baseline = vs_base$max_abs, max_rel_diff_baseline = vs_base$max_rel,
          identical_to_previous = vs_prev$identical,
          max_abs_diff_previous = vs_prev$max_abs, max_rel_diff_previous = vs_prev$max_rel,
          note = paste(unique(notes[nzchar(notes)]), collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
    prev_i <- i
  }
  golden_df <- do.call(rbind, rows)
  utils::write.csv(golden_df, file.path(out_dir, "golden.csv"), row.names = FALSE)

  cell <- function(sub, which) {
    ident <- sub[[paste0("identical_to_", which)]]
    rel   <- sub[[paste0("max_rel_diff_", which)]]
    if (any(startsWith(sub$note, "error:"))) return("error")
    if (all(is.na(ident))) return("–")
    if (all(ident, na.rm = TRUE)) return("=")
    worst <- suppressWarnings(max(rel, na.rm = TRUE))
    if (!is.finite(worst)) return("changed")
    sprintf("%.1e", worst)
  }
  identity_table <- function(which) {
    wide <- data.frame(commit = paste0("`", commits$sha[have_golden], "`"), check.names = FALSE)
    for (s in scen) {
      wide[[s]] <- vapply(commits$sha[have_golden], function(sha) {
        cell(golden_df[golden_df$commit == sha & golden_df$scenario == s, ], which)
      }, character(1))
    }
    tw <- as.data.frame(t(wide[, -1L, drop = FALSE]), stringsAsFactors = FALSE)
    names(tw) <- wide$commit
    cbind(scenario = rownames(tw), tw)
  }

  add("## Numerical identity", "",
      "Outputs of the golden scenarios (`tests/testthat/helper-golden-scenarios.R`),",
      "compared bit for bit. `=` means every output element is `identical()`;",
      "otherwise the cell is the largest relative difference. Element-level",
      "detail is in `golden.csv`.", "",
      "### Against the previous timed commit", "",
      "Shows which commit moved which numbers.", "")
  add(md_table(identity_table("previous")), "")
  add("### Against the baseline", "")
  add(md_table(identity_table("baseline")), "")
}

# --- tests -------------------------------------------------------------------

test_files <- list.files(file.path(out_dir, "tests"), pattern = "\\.csv$", full.names = TRUE)
if (length(test_files)) {
  # commit must stay character: a short sha like 5e17080 otherwise parses as 5e+17080
  tests_df <- do.call(rbind, lapply(test_files, utils::read.csv, stringsAsFactors = FALSE,
                                    colClasses = c(commit = "character")))
  tests_df$idx <- commits$idx[match(tests_df$commit, commits$sha)]
  tests_df <- tests_df[order(tests_df$idx, tests_df$suite, tests_df$file), ]
  utils::write.csv(tests_df, file.path(out_dir, "tests.csv"), row.names = FALSE)

  counts <- function(sha, suite) {
    sub <- tests_df[tests_df$commit == sha & tests_df$suite == suite, ]
    if (!nrow(sub)) return("–")
    n <- function(st) sum(sub$status == st)
    sprintf("%d / %d / %d / %d", n("pass"), n("fail"), n("error"), n("skip"))
  }
  add("## Tests", "",
      "`own`: the commit's own test suite, as it stood at that commit.",
      "`tip`: the final harness's suite run against the commit's package, so",
      "every row is judged by the same tests.", "",
      "Counts are test_that() blocks: pass / fail / error / skip.", "")
  add(md_table(data.frame(
    commit = paste0("`", commits$sha, "`"),
    own = vapply(commits$sha, counts, character(1), suite = "own"),
    tip = vapply(commits$sha, counts, character(1), suite = "tip"),
    subject = subject_short(commits$subject, 50L),
    check.names = FALSE
  )), "")

  # Where does each tip-suite test change status along the lineage?
  tip <- tests_df[tests_df$suite == "tip", ]
  if (nrow(tip)) {
    tip$key <- paste(tip$file, tip$test, sep = " :: ")
    changes <- character()
    for (k in unique(tip$key)) {
      sub <- tip[tip$key == k, ]
      sub <- sub[order(sub$idx), ]
      flips <- which(sub$status[-1L] != sub$status[-nrow(sub)]) + 1L
      if (length(flips)) {
        steps <- sprintf("%s → %s at `%s`", sub$status[flips - 1L], sub$status[flips], sub$commit[flips])
        changes <- c(changes, sprintf("| %s | %s | %s | %s |",
                                      sub$file[1], sub$test[1], sub$status[1], paste(steps, collapse = "; ")))
      }
    }
    add("### Tip-suite tests that change status along the lineage", "")
    if (length(changes)) {
      add("| file | test | at baseline | changes |", "|---|---|---|---|", changes, "")
    } else {
      add("None: every tip-suite test has the same status at every commit.", "")
    }
  }
}

writeLines(md, file.path(out_dir, "lineage-report.md"))
message("wrote ", file.path(out_dir, "lineage-report.md"))
