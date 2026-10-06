# Per-commit Docker images

Every commit in the perf lineage (`main` .. `perf/baseline-tests`) gets its own
image with hypervolume installed from that commit's source. All of them share
one base image, so R, BLAS/LAPACK, compilers and dependency versions are
identical. Any difference in speed, numbers or test results between two images
comes from the commits.

This replaces the checkout → `rm src/*.o` → `R CMD INSTALL` → test → checkout
loop described in the [main README](../README.md#tests). Each image is installed from a clean
`git archive` of its commit, so there are never stale object files.

## Requirements

Docker (or Podman: `DOCKER=podman`), about 3 GB of disk, and an internet
connection for the first build.

## Run the whole lineage

```sh
docker/lineage.sh
```

This builds the images (around 15 minutes the first time, mostly the shared
dependency image), runs everything in each one, and writes
`lineage-out/lineage-report.md`. Steps can be run separately:

```sh
docker/lineage.sh build     # images only
docker/lineage.sh run       # tests, golden values, benchmarks
docker/lineage.sh report    # rebuild the report from existing results
```

Settings are environment variables, documented at the top of `lineage.sh`. The
ones you are most likely to want:

```sh
BENCH_TIMES=50 docker/lineage.sh run       # more iterations per bench scenario
HYPERVOLUME_TEST_SLOW=true docker/lineage.sh run   # include the slow tests
```

## Step into one commit

Images are tagged by short commit hash. The commit's full source, tests
included, is at `/pkg`.

```sh
docker images hypervolume-tools                      # list them
docker run --rm -it hypervolume-tools:008415f R      # R with that commit installed
docker run --rm hypervolume-tools:008415f \
  Rscript -e 'testthat::test_dir("/pkg/tests/testthat", package = "hypervolume", load_package = "installed")'
```

Two images can be compared side by side, for example by running the same R
snippet in the image before and after a perf commit.

## What gets run in each image

| | What | Why |
|---|---|---|
| **own tests** | the commit's own `tests/testthat`, as it was at that commit | the test lineage as it actually happened |
| **tip tests** | the final test suite, run against the commit's package | every commit judged by the same tests; shows which tests turn green where |
| **golden** | the scenarios in `helper-golden-scenarios.R`; raw outputs saved | bit-for-bit comparison between commits |
| **bench** | the scenarios in `tools/bench.R`, every iteration timed | speed-up per commit |

The harness (tests, `tools/bench.R` and `docker/R/`) is taken from the checkout
you run `lineage.sh` in, and mounted read-only into each image.

Golden values and benchmarks are only produced for commits that change what
gets installed (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`, ...). Other commits
install a byte-identical package. A few non-perf commits are flagged as code
changes anyway (one adds `testthat` to `DESCRIPTION`, two add and remove a stray
`R/.DS_Store`). They act as controls: they should show `=` for every number and
roughly 1.0× for every timing.

## Reading the report

- **Speed**: median ms per call, with the speed-up against the baseline.
  Interquartile ranges are in `bench.csv`.
- **Numerical identity**: `=` means every output is `identical()` to the
  compared commit. Otherwise the cell shows the largest relative difference.
  The "previous commit" table shows which commit moved which numbers.
  Deliberate changes, like the `calculate_density` normalisation fix, should
  show up only in the commit that made them.
- **Tests**: pass / fail / error / skip counts per suite, plus every tip-suite
  test whose status changes along the lineage, and the commit where it changes.

`bench.csv`, `golden.csv` and `tests.csv` in `lineage-out/` hold the detail
behind each table.

## Environment choices

- **R 4.5.0**, the version the fixtures in `tests/` and `tools/` were recorded
  with. Dependencies come from the dated Posit Package Manager snapshot that
  `rocker/r-ver:4.5.0` pins, so a rebuild months later gets the same versions.
- **Reference BLAS/LAPACK**, not the OpenBLAS that rocker ships with. CRAN's
  macOS and Windows builds use reference BLAS, so these timings match what most
  users will see. OpenBLAS would flatter the rewrites that move work into
  `dgemm`. To measure under OpenBLAS, remove the marked block in
  `Dockerfile.deps`.
- **Single-threaded** (`OMP_NUM_THREADS=1`), for deterministic numbers and
  stable timings.

Timings are only comparable within one run on one machine. Under Docker Desktop
on macOS everything runs inside a Linux VM, so absolute times will differ from
native R, but the ratios between commits hold. Bit-for-bit results can also
differ between CPU architectures (x86_64 vs arm64), so compare commits within
one run rather than against results from another machine.
