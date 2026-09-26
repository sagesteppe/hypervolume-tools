# Testing this fork

Baseline test suite for the performance work described in
`hypervolume_speedup_notes.md`. Nothing in the package source has been changed
yet — `main` is the untouched 3.1.7 download, and this branch adds tests only.

## Layout

```
tests/
  testthat.R                       # entry point
  testthat/
    helper-hv.R                    # fixtures + independent reference implementations
    helper-golden-scenarios.R      # scenario definitions shared with tools/make_golden.R
    test-01-fastpdist2.R           # proposal #9
    test-02-kdtree.R               # proposals #1, #5, #6
    test-04-calculate-density.R    # proposal #2
    test-05-volume-and-threshold.R # proposal #3
    test-06-samplers.R             # proposal #7
    test-07-gaussian-integration.R # proposals #2 + #3 end to end
    test-08-set-operations.R       # proposal #8
    test-09-estimate-probability.R # proposal #4
    test-10-golden.R               # cross-cutting regression snapshot
    fixtures/                      # golden-baseline.rds lives here (generated)
tools/
  make_golden.R                    # generates the snapshot; run on `main`
```

Numbering has a gap at 03 on purpose, so a file can be inserted next to the
kd-tree tests later without renumbering everything.

## Running

```r
# install once, then iterate
R CMD INSTALL .
testthat::test_local()

# or from a shell
Rscript -e 'devtools::test()'
```

### C++ changes: force a clean rebuild

R's default `Makevars` **does not track `.h` dependencies**. Editing
`src/KDTree.h` (or `MyHeaps.h`) will not by itself trigger a rebuild of
the `.cpp` files that include it — `make` compares source-file mtimes,
sees nothing changed, and re-uses the stale `.o`. `R CMD INSTALL .`
will then happily re-link the old object into a "new" shared library
and both tests and bench numbers will silently reflect the previous
build.

**When editing anything under `src/`, always:**

```sh
rm -f src/*.o src/*.so && R CMD INSTALL .
```

before running the test suite or `tools/bench.R --compare`. Without
this, `--compare` will show phantom "no improvement" numbers. This
cost about an hour of debugging during the KDTree work.

A handful of tests are gated behind environment variables so a fresh checkout
runs green:

```r
Sys.setenv(HYPERVOLUME_TEST_SLOW = "true")   # slow (>1s) tests
Sys.setenv(HYPERVOLUME_TEST_OMP  = "true")   # OpenMP determinism check
testthat::test_local()
```

`HYPERVOLUME_TEST_SLOW` covers the kd-tree build-scaling guard and the
"recognisable 2D Gaussian" sanity check. `HYPERVOLUME_TEST_OMP` covers the
thread-count determinism test in `test-02` — see the acceptance-criteria
checklist below.

## The golden baseline

`test-10-golden.R` compares against numbers recorded from the *unmodified*
package. Generate them before starting any optimisation work:

```sh
git checkout main
R CMD INSTALL .
Rscript tools/make_golden.R
git checkout perf/baseline-tests
git add -f tests/testthat/fixtures/golden-baseline.rds
git commit -m "test: golden baseline from unmodified 3.1.7"
```

If the fixture is missing the golden tests skip rather than fail, so a fresh
checkout still runs green. The generator refuses to overwrite an existing
fixture unless passed `--force`, because re-baselining against modified code
silently destroys the regression check.

The fixture records the R version, platform and BLAS alongside the values.
`fastPdist2` goes through BLAS `dgemm`, so a different BLAS can move the last
couple of bits — if the golden tests fail only in the `fastpdist` and
`calculate_density` scenarios and only at ~1e-12, check the metadata before
assuming a real regression.

## What is and is not covered

Covered, because a proposal touches it:

| Area | Functions |
|---|---|
| kd-tree | `kdtree_build`, `kdtree_ball_query_multiple`, `kdtree_ball_query_id_multiple`, `kdtree_range_query_multiple`, `evalfspherical`, `evalfrectangular`, `repmat` |
| density | `calculate_density`, `hypervolume_gaussian` |
| thresholding | `hypervolume_threshold`, `volume_calculation` |
| sampling | `sample_ellipsoid`, `ellipsoid_volume`, `ellipsoid_inverse_weight`, `nball_volume`, `sample_model_ellipsoid`, `sample_model_rejection` |
| set operations | `hypervolume_set`, `hypervolume_overlap_statistics`, `hypervolume_join`, `hypervolume_distance` |
| inference | `hypervolume_estimate_probability`, `hypervolume_inclusion_test` |
| distance | `fastPdist2` |

Deliberately **not** covered — no proposal touches them, and testing them would
have meant a much larger suite for no safety gain:

- the whole occupancy family (`hypervolume_n_occupancy*`, `occupancy_*`,
  `get_occupancy_*`)
- resampling and bootstrap (`hypervolume_resample`, `hypervolume_permute`,
  `hypervolume_n_resample`, `bootstrap*`)
- other constructors (`hypervolume_box`, `hypervolume_svm`,
  `hypervolume_general_model`) — note these share `sample_model_ellipsoid`
  and `sample_model_rejection`, which *are* covered
- plotting, `hypervolume_save_animated_gif`, projection, segmentation, holes,
  `expectation_*`

Measure the real number rather than trusting this table:

```r
covr::report(covr::package_coverage(type = "tests"))
```

## Tests that pin behaviour rather than endorse it

A few tests assert things that look like bugs. They are there so a
performance refactor cannot change them *by accident* — if any is fixed
deliberately, update the test and say so in the commit.

0. **`hypervolume_threshold(uniform.density = TRUE)` returns more values than
   points** (`test-05`). The line is
   `hv_new@ValueAtRandomPoints <- rep(1, length(hv_new@ValueAtRandomPoints))`,
   but at that point `hv_new` is still a copy of `hv` with only `@RandomPoints`
   subset, so `length()` is the *original* n. Every `uniform.density = TRUE`
   result therefore violates the `length(@ValueAtRandomPoints) ==
   nrow(@RandomPoints)` invariant. Found while reviewing the suite, not in the
   original performance pass. Cheap to fix, unrelated to performance — worth
   its own commit.

1. **`calculate_density(weight = NULL)` returns `numeric(0)`**
   (`test-04`). The signature defaults `weight` to `NULL`, but the body does
   `weight[i] * dmvnorm(...)`; `NULL[i]` is `NULL`, so the accumulator collapses.
   `hypervolume_gaussian` always fills `weight` first, so this never fires in
   practice.

2. **`hypervolume_estimate_probability` returns `NaN` when a query point
   coincides exactly with an observation** (`test-09`). With the default
   `weight.exponent = -1`, a distance of zero produces an `Inf` weight, and
   `sum(prob * Inf) / sum(Inf)` collapses to `NaN` rather than the limit
   (which would be `prob[coincident_idx]`). `hypervolume_inclusion_test`
   forwards to this and propagates the `NaN`. Almost impossible to hit in
   real use — query points are generally continuous, and coincidence with
   a Monte-Carlo `RandomPoints` sample is measure-zero — so this is
   deliberately left in place. The test dodges it by offsetting the query
   by `1e-4`.

### Fixed and removed from this list

- **`calculate_density` integrated to `1/n` instead of `1`.** The trailing
  `density / n_points` was a double-normalisation (weight is already
  sum-to-1 in the `hypervolume_gaussian` path). Fixed; `test-04` now
  pins integral == 1 as a proper probability density. Volumes are
  unaffected (thresholding is quantile-based). Downstream consumers of
  `@ValueAtRandomPoints` in absolute units — chiefly
  `hypervolume_estimate_probability` — now return values with proper
  probability-density units instead of being off by 1/n. Golden fixture
  re-baselined on the corrected values.

- **`fastPdist2` returned `NaN` for near-coincident points.** The
  `|a|² + |b|² − 2a·b` expansion can go slightly negative from
  rounding. Fixed by clamping negatives to zero before `sqrt`; `test-01`
  now asserts no `NaN` on the diagonal of `fastPdist2(A, A)`.

## Acceptance-criteria checklist for the perf proposals

Three tests are wired up but currently skipped. They are the pass/fail gates
for the proposals they name — un-skipping them (and passing) is part of the
commit that lands the proposal, not a follow-up.

| Test | Location | Un-skip in the commit that lands |
|---|---|---|
| `fastPdist2 returns real distances for near-coincident points` | `test-01-fastpdist2.R` | proposal #9 (fastPdist2 clamp) |
| `tree build cost is not quadratic in n` | `test-02-kdtree.R` | proposal #1 (KDTree workarray) |
| `ball-query results are independent of OpenMP thread count` | `test-02-kdtree.R` | proposal #6 (OpenMP over ball queries) — also add `hypervolume:::omp_set_threads()` and wire `HYPERVOLUME_TEST_OMP=true` into CI |

If a proposal lands without flipping its skip, the regression net for that
proposal isn't active — the CI would stay green whether or not the fix worked.

## Notes for whoever lands the optimisations

- `test-07` has a scale-equivariance check (`volume` scales as `c^d` when data
  and bandwidth are both multiplied by `c`). Because acceptance is
  deterministic and thresholding is quantile-based, it should hold to
  floating-point precision, not just Monte Carlo precision. It is the sharpest
  single check on the `calculate_density` rewrite getting the normalising
  constant right.
- `test-04` asserts `calculate_density` returns a bare numeric vector with no
  `dim`. A `%*%`-based rewrite returns an n×1 matrix unless it is unwrapped,
  and that would survive several lines downstream before failing confusingly.
- `test-06` pins `sample_model_ellipsoid` on three cases where the volume is
  *analytic* rather than noisy (one ellipsoid; two coincident; two far apart).
  Those are the checks that hoisting the kd-tree build out of the `while` loop
  must not disturb.
