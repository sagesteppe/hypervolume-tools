# hypervolume performance work

Proposed speed-ups to [hypervolume](https://github.com/bblonder/hypervolume)
3.1.7, one commit per change, with the tests that check each one and a way to
reproduce every measurement.

`main` is the untouched 3.1.7 source. `perf/baseline-tests` first adds a test
suite, golden values and benchmarks recorded against that unmodified code, then
lands the changes one at a time. Each commit says what it changes and why.

## Results

Measured with [`docker/lineage.sh`](docker/README.md): R 4.5.0, reference
BLAS/LAPACK, single thread, x86_64 Linux, 20 iterations per scenario. "Before"
is the median over the four builds whose installed code is identical to 3.1.7.
"After" is the last perf commit (`008415f`).

| scenario (`tools/bench.R`) | before (ms) | after (ms) | speed-up |
|---|---:|---:|---:|
| `hypervolume_estimate_probability`, 100 queries | 18.3 | 2.1 | 8.9× |
| `calculate_density`, 200 × 1000 | 18.5 | 2.3 | 8.0× |
| `kdtree_build`, n = 10 000 | 21.8 | 6.7 | 3.3× |
| `evalfspherical`, 10 000 points × 500 queries | 32.6 | 13.2 | 2.5× |
| `kdtree_build`, n = 5000 | 6.7 | 2.7 | 2.4× |
| `hypervolume_set`, 2 × 5000 points | 27.8 | 11.8 | 2.4× |
| `hypervolume_gaussian`, *Quercus rubra* subset (end to end) | 651 | 451 | 1.4× |
| `kdtree_build`, n = 500 | 0.27 | 0.20 | 1.4× |
| `hypervolume_threshold`, 1000 thresholds | 212 | 181 | 1.2× |
| `fastPdist2`, 500 × 500 | 1.5 | 1.6 | 0.9× |

Timings of identical code vary by up to 1.2× between builds on most scenarios,
and up to 1.5× on the sub-2 ms ones. Treat anything under about 1.3× as noise.
The `fastPdist2` change is a correctness fix (proposal #9), not a speed-up; it
costs nothing measurable.

Numerically, six of the nine golden scenarios are bit-identical to 3.1.7 at
every commit: kd-tree ball counts, threshold statistics, sampler volume, 3D
Gaussian hypervolume, set operations and `fastPdist2`. The rest move only where
expected:

- **`17400e8`** (proposal #2): `calculate_density` and the 2D Gaussian
  hypervolume differ by at most 1.9 × 10⁻¹⁵ relative. That's summation order
  in the matrix product.
- **`bc33d9a`** (proposal #4): `hypervolume_estimate_probability` differs by
  at most 4.2 × 10⁻¹⁰ relative, because distances now come from `fastPdist2`'s
  |a|² + |b|² − 2a·b expansion rather than `pdist`. That's within the golden
  test's 10⁻⁸ tolerance.
- **`a26da65`**: the deliberate normalisation fix described below.

Tests run against every commit flip from fail to pass only at the commit that
fixes them: the `fastPdist2` near-coincident test at `56a4c85`, and the three
`calculate_density` normalisation tests plus the golden test at `a26da65`.

## Reproducing the results

Every commit can be installed into its own Docker image, and the tests, golden
values and benchmarks run in each. See [`docker/README.md`](docker/README.md).

```sh
docker/lineage.sh        # writes lineage-out/lineage-report.md
```

## The changes

"Proposal #N" in commit messages and test comments refers to this list.

| # | Where | Change | Commit |
|---|---|---|---|
| 1 | `src/KDTree.h` | Tree build cleared an O(n) scratch array at every node, making it O(n²). Now clears only the entries it wrote. | `e1c835f` |
| 2 | `R/hypervolume_gaussian.R` | `calculate_density` looped `dmvnorm` over every observation. Now one blocked matrix product. | `17400e8` |
| 3 | `R/hypervolume_threshold.R` | Built up to 1000 full hypervolumes to sweep thresholds. Now one sort plus cumulative sums. | `187dbc1` |
| 4 | `R/hypervolume_estimate_probability.R` | One `pdist()` call per query point. Now blocked through `fastPdist2`. | `bc33d9a` |
| 5 | `src/KDTree.h`, `src/convert.cpp` | Points were copied three times into nested vectors. Now one flat row-major array. | `008415f` |
| 6 | kd-tree ball queries | OpenMP over independent queries. **Not done**; its test is scaffolded and skipped. | – |
| 7 | `R/samplers.R` | `ellipsoid_inverse_weight` rebuilt the same kd-tree on every sampler iteration. Built once now. | `4f92bcd` |
| 8 | `R/hypervolume_set.R` | `unique()` and `rbind()` on data.frames. Now stays in matrices. | `c78c4ee` |
| 9 | `src/fastPdist2.cpp` | Rounding could make squared distances slightly negative, giving `NaN` after `sqrt`. Clamped at zero. | `56a4c85` |

Considered and rejected: truncating the Gaussian density sum to nearby kernels.
Blonder et al. (2018, §2.3) identify truncating the density estimate as the
source of bias that the Gaussian method was designed to avoid.

## Changes to results

All perf commits are intended to reproduce the 3.1.7 numbers. The golden tests
hold them to that, and the lineage report shows bit for bit where any output
moves. Two commits change results on purpose:

**`a26da65`: `calculate_density` integrated to 1/m instead of 1.** Blonder et
al. (2014, Box 2) define the estimator with normalising constant
α = 1/(m·∏h), so exactly one factor of 1/m. The code applied it twice: once
through `weight`, which `hypervolume_gaussian` normalises to sum to 1, and again
through a trailing `density / n_points`.

- Volumes are unchanged. Thresholding is a ratio of sums, and the sampler's
  acceptance test only checks density > 0, so neither depends on the density's
  scale.
- `@ValueAtRandomPoints` is now a probability density, m times larger than
  before. Before the fix it could not be compared between hypervolumes built
  from different numbers of observations.
- `hypervolume_estimate_probability`, which reads those values directly,
  changes by the same factor.

The fix is its own commit, after the BLAS rewrite (proposal #2) had already
reproduced the old scaling exactly. A numbers change is then never hidden
inside a speed-up commit.

**`56a4c85`: `fastPdist2` returned `NaN` for near-coincident points.** It now
returns 0.

## Tests

```
tests/testthat/
  helper-hv.R                    fixtures and independent reference implementations
  helper-golden-scenarios.R      scenarios shared with tools/make_golden.R
  test-01 … test-09              one file per area a proposal touches
  test-10-golden.R               regression snapshot against the unmodified package
  fixtures/golden-baseline.rds
tools/
  make_golden.R                  regenerates the golden snapshot
  bench.R                        benchmark scenarios, one or more per proposal
```

Without Docker:

```r
R CMD INSTALL .
testthat::test_local()
```

After any change under `src/`, run `rm -f src/*.o src/*.so` before installing.
R's default `Makevars` doesn't track header dependencies, so editing
`KDTree.h` alone can leave a stale object file in the build. The Docker images
always build from a clean copy.

Some tests are skipped unless an environment variable is set:

- `HYPERVOLUME_TEST_SLOW=true`: the kd-tree build-scaling guard and the
  "recognisable 2D Gaussian" sanity check.
- `HYPERVOLUME_TEST_OMP=true`: thread-count determinism of ball queries. Only
  meaningful once proposal #6 lands.

The golden fixture records the R version, platform and BLAS it was made with.
`fastPdist2` goes through BLAS `dgemm`, so a different BLAS can move the last
few bits. If only the `fastpdist` and `calculate_density` scenarios fail, and
only at around 1e-12, compare the BLAS in the fixture's metadata with yours
before assuming a regression.

Coverage is limited to the functions the proposals touch. The occupancy
family, resampling and bootstrap, the other constructors (`hypervolume_box`,
`hypervolume_svm`, `hypervolume_general_model`), plotting and projection are
not tested here.

### Tests that pin existing behaviour

A few tests assert behaviour that looks like a bug. They exist so that a
performance change can't alter it by accident. None of these is fixed here,
because each is a separate change from the speed work.

0. **`hypervolume_threshold(uniform.density = TRUE)` returns more values than
   points** (`test-05`). `@ValueAtRandomPoints` is filled with
   `rep(1, length(hv_new@ValueAtRandomPoints))` while `hv_new` still holds the
   original, unthresholded vector. So `length(@ValueAtRandomPoints)` no longer
   equals `nrow(@RandomPoints)`.
1. **`calculate_density(weight = NULL)` returns `numeric(0)`** (`test-04`).
   `NULL` is the default, but the body indexes `weight[i]`.
   `hypervolume_gaussian` always supplies `weight`, so users can't reach this.
2. **`hypervolume_estimate_probability` returns `NaN` when a query point
   exactly equals an observation** (`test-09`). A zero distance gives an `Inf`
   weight, and `sum(prob * Inf) / sum(Inf)` is `NaN`. This needs the query to
   hit a Monte Carlo sample point exactly, which almost never happens with
   continuous data.
