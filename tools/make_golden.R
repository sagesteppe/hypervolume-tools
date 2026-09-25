#!/usr/bin/env Rscript
#
# Generate the golden-value baseline used by tests/testthat/test-10-golden.R.
#
# RUN THIS ON THE BASELINE COMMIT (branch `main`), BEFORE ANY PERFORMANCE WORK:
#
#     git checkout main
#     R CMD INSTALL .
#     Rscript tools/make_golden.R
#
# Then switch back to the working branch. The fixture records what the
# unmodified package produces; the tests assert that the optimised package
# still produces it.
#
# Re-running this on a modified package will happily overwrite the baseline
# with the modified numbers, which defeats the point -- so the script refuses
# unless the fixture is absent or --force is passed.

suppressPackageStartupMessages({
  library(testthat)
  library(hypervolume)
})

args <- commandArgs(trailingOnly = TRUE)
# Not named `force`: helper-hv.R is sourced into this same environment and its
# with_seed()/quietly() both call base::force().
force_overwrite <- "--force" %in% args

root <- normalizePath(".", mustWork = TRUE)
helper_dir <- file.path(root, "tests", "testthat")
if (!dir.exists(helper_dir)) {
  stop(
    "tests/testthat not found. Run this from the package root:\n",
    "    Rscript tools/make_golden.R"
  )
}

# testthat::test_path() only resolves inside a running test, so point it here.
# No on.exit(): at top level under Rscript there is no return context to
# attach to, so the handler would be silently discarded. Restored explicitly
# at the end instead.
old_wd <- setwd(helper_dir)

source("helper-hv.R", local = FALSE)
source("helper-golden-scenarios.R", local = FALSE)

fixture_dir <- file.path(helper_dir, "fixtures")
dir.create(fixture_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(fixture_dir, "golden-baseline.rds")

if (file.exists(out_path) && !force_overwrite) {
  stop(
    "golden-baseline.rds already exists.\n",
    "Overwriting it re-baselines against whatever the package currently does, ",
    "which silently destroys the regression check.\n",
    "Pass --force only if you are deliberately re-baselining."
  )
}

git_sha <- tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
git_branch <- tryCatch(
  system2("git", c("rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)

if (!identical(git_branch, "main")) {
  message(
    "NOTE: current branch is '", git_branch, "', not 'main'. ",
    "Make sure this is really the unmodified baseline."
  )
}

scenarios <- golden_scenarios()
values <- list()

for (nm in names(scenarios)) {
  message("  running scenario: ", nm)
  t0 <- proc.time()[["elapsed"]]
  values[[nm]] <- scenarios[[nm]]()
  message(sprintf("    done in %.2fs", proc.time()[["elapsed"]] - t0))
}

golden <- list(
  values = values,
  meta = list(
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    package_version = as.character(utils::packageVersion("hypervolume")),
    r_version = R.version.string,
    platform = R.version$platform,
    blas = tryCatch(sessionInfo()$BLAS, error = function(e) NA_character_),
    git_sha = git_sha,
    git_branch = git_branch,
    seed = HV_SEED
  )
)

saveRDS(golden, out_path)

message("")
message("wrote ", out_path)
message("  package ", golden$meta$package_version,
        "  |  ", golden$meta$git_branch, "@", golden$meta$git_sha)
message("  ", golden$meta$r_version)
message("")
message("Commit this file so the comparison survives across checkouts:")
message("  git add -f tests/testthat/fixtures/golden-baseline.rds")

setwd(old_wd)
