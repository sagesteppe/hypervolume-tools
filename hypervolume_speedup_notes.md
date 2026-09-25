# `hypervolume` — low-hanging performance fruit

Reviewed `bblonder/hypervolume` @ master (DESCRIPTION v3.1.7, 2026-09-17).
Note: GitHub clone was blocked by the sandbox proxy, so files were read individually
from `raw.githubusercontent.com`. Everything below cites file + function; line numbers
will be close but not exact.

Ranked by (expected speedup) × (1 / risk). The top three are, in my view, the whole game.

---

## 1. `src/KDTree.h` — `build_recursively()` is accidentally O(n²)

**This is the big one.** In `build_recursively`, every internal node does:

```cpp
// initialize the "partition" array
for (int i = 0; i < npoints; i++)
    workarray[i] = -1;                                  // <-- O(npoints), EVERY node
for (unsigned int i = 0; i < pidx.size(); i++)
    workarray[ sortidx[dim][pidx[i]] ] = pidx[i];
```

`npoints` is the size of the *whole* dataset, not the current subset, and there are
~`npoints` internal nodes. So tree construction spends exactly `n(n-1)` writes just
re-zeroing a scratch array.

I simulated the exact op counts (the split is deterministic: `ceil(m/2)` / `floor(m/2)`):

| n | wasted clear-writes | actual useful work | ratio |
|---|---|---|---|
| 1,000 | 999,000 | 44,754 | 22× |
| 10,000 | 99,990,000 | 786,351 | 127× |
| 100,000 | 9,999,900,000 | 11,770,059 | **850×** |
| 1,000,000 | 999,999,000,000 | 169,888,045 | **5,886×** |

This matters because `kdtree_build()` is called on `hv@RandomPoints`, which is routinely
10⁵–10⁶ rows (`samples.per.point` defaults to `10^(3+sqrt(d))/n`). Every
`hypervolume_set()`, `hypervolume_inclusion_test()`, and every iteration of the
`sample_model_ellipsoid()` accept/reject loop pays it.

**Fix — two lines, no behaviour change.** Only clear what you wrote, and skip the whole
thing in the branch that doesn't use `workarray`:

```cpp
int pivot = -1;
if ((double)pidx.size()*log((double)pidx.size()) < npoints) {
    // heapsort branch — workarray is never read here, so don't touch it at all
    Larray.resize(pidx.size()/2 + (pidx.size()%2==0 ? 0:1), -1);
    Rarray.resize(pidx.size()/2, -1);
    heapsort(dim, pidx, pidx.size());
    std::copy(pidx.begin(), pidx.begin()+Larray.size(), Larray.begin());
    std::copy(pidx.begin()+Larray.size(), pidx.end(),   Rarray.begin());
    pivot = pidx[(pidx.size()-1)/2];
} else {
    for (size_t i = 0; i < pidx.size(); i++)
        workarray[ sortidx[dim][pidx[i]] ] = pidx[i];
    ... existing scan ...
    for (size_t i = 0; i < pidx.size(); i++)          // restore the invariant
        workarray[ sortidx[dim][pidx[i]] ] = -1;
}
```

The `else`-branch scan (`for i < npoints`) is also O(npoints) per node, but it only fires
near the root (the heapsort branch takes over once `m log m < n`), so it totals ~O(n log n)
and is fine. Only the unconditional clear needs to go.

**Expected:** tree build goes from quadratic to ~n log n. Order-of-magnitude-plus on any
realistic hypervolume; the difference between "slow" and "appears to hang" at n ≥ 10⁶.

Two smaller things in the same file while you're there:

- `check_border_distance(int, int, vector<double> pnt, ...)` and `closest_point(vector<double> p)`
  take the query point **by value**, and `check_border_distance` recurses — so the point is
  deep-copied at every level of the tree. Make them `const vector<double>&`.
- `distance_squared()` declares `double N = a.size();` then loops `for (int i=0; i<N; i++)`.
  The int↔double compare each iteration blocks vectorisation. Use `size_t N`.
- `k_closest_points()` calls `std::reverse(idxs.begin(), idxs.end())` twice in a row — a
  no-op. Looks like a bug rather than a perf issue, but it means results come back in the
  opposite order from what the comment claims.

---

## 2. `R/hypervolume_gaussian.R` — `calculate_density()` is the wrong algorithm

```r
for (i in index_vals) {
    density = density + weight[i] * dmvnorm(s, mean = means[i, ], sigma = kernel_sd^2 * diag(d))
}
```

Per data point, per call, this:

- allocates `diag(d)` and a fresh d×d `sigma`,
- hands it to `mvtnorm::dmvnorm`, which Cholesky-factorises that (diagonal!) matrix,
- and does it all in an interpreted R loop.

The covariance is diagonal and *constant across the loop*. The whole thing is a weighted
sum of spherical Gaussians, i.e. one matrix product:

```r
calculate_density <- function(s, means, kernel_sd, weight=NULL, chunk=2e4, ...) {
  d <- ncol(means); n <- nrow(means)
  if (is.null(weight)) weight <- rep(1/n, n)
  Z  <- sweep(as.matrix(s),     2, kernel_sd, "/")   # scale once
  M  <- sweep(means,            2, kernel_sd, "/")
  logconst <- -0.5*d*log(2*pi) - sum(log(kernel_sd))
  out <- numeric(nrow(Z))
  for (k in split(seq_len(nrow(Z)), ceiling(seq_len(nrow(Z))/chunk))) {
    D2 <- fastPdist2(Z[k, , drop=FALSE], M)^2        # already in the package!
    out[k] <- exp(logconst) * (exp(-0.5*D2) %*% weight)
  }
  out / n
}
```

That's one BLAS `dgemm` per chunk instead of `nrow(data)` `dmvnorm` calls. Expect 1–2
orders of magnitude on `hypervolume_gaussian()`, which is the package's flagship entry point.

Two related notes:

- The chunking machinery is quadratic for no benefit: `chunks = ceiling((1:n_points)/chunksize)`
  then `which(chunks == chunk_this)` inside the loop scans all `n_points` every chunk —
  `n²/chunksize` comparisons purely to drive a progress bar. Use `seq()` arithmetic.
- **Possible correctness issue, worth a look:** `weight` is normalised to sum to 1, and then
  the result is divided by `n_points` again — so the returned "density" is off by a factor of
  `n`. It cancels out downstream because everything is quantile-thresholded, but
  `ValueAtRandomPoints` won't be a probability density in the units you'd expect.
- Bigger algorithmic win if you want it: sampling is truncated at `sd.count * bandwidth`
  (default 3 SD), so kernels further than that contribute ~1e-4 of a kernel peak. Restricting
  the sum to neighbours within `sd.count * bandwidth` — you already have the kd-tree — turns
  an O(n_data × n_samples) sum into roughly O(n_samples × k). More invasive, so I'd land #1
  and the BLAS rewrite first and re-measure.

---

## 3. `R/hypervolume_threshold.R` — 1000 full passes + 1000 materialised hypervolumes

`hypervolume_gaussian()` calls this with `num.thresholds = 1000` and then keeps **one**
result. The loop body does, per threshold:

```r
hv_new@RandomPoints      <- hv@RandomPoints[hv@ValueAtRandomPoints > thresholds[i], , drop=F]
hv_new@Volume            <- length(which(hv@ValueAtRandomPoints > thresholds[i])) / hv@PointDensity
hv_new@ValueAtRandomPoints <- hv@ValueAtRandomPoints[hv@ValueAtRandomPoints > thresholds[i]]
probability_sums[i]      <- sum(hv@ValueAtRandomPoints[hv@ValueAtRandomPoints > thresholds[i]])
```

Four independent O(n) comparisons of the same predicate, plus a full subset copy of the
n×d point matrix, ×1000. At n = 5×10⁵, d = 4 that's ~16 GB of transient allocation to
produce one object.

**Fix:** sort once, get every threshold's statistics from cumulative sums, and only
materialise the hypervolume actually returned.

```r
o  <- order(hv@ValueAtRandomPoints)
v  <- hv@ValueAtRandomPoints[o]
# counts above each threshold
idx    <- findInterval(thresholds, v)          # number of values <= threshold
counts <- length(v) - idx
psums  <- sum(v) - c(0, cumsum(v))[idx + 1]
volumes <- counts / hv@PointDensity
```

Then build the single `Hypervolume` for `quantile.obtained.index` at the end. When
`quantile.requested` is `NULL` and the caller really does want all 1000, keep the loop but
compute the predicate once per iteration and reuse it.

**Expected:** effectively removes `hypervolume_threshold` from the profile, and removes the
memory spike that makes `hypervolume_gaussian` fail on large inputs.

---

## 4. `R/hypervolume_estimate_probability.R` — one `pdist()` call per test point

```r
probabilities <- sapply(1:nrow(points), function(i) {
    distances <- pdist::pdist(points[i, , drop=FALSE], hv_points_ss)@dist
    weights   <- distances^weight.exponent
    sum(prob_ss * weights) / sum(weights)
})
```

`nrow(points)` separate R-level calls, each with its own S4 construction and allocation, to
do what is one blocked matrix operation. `fastPdist2()` is already compiled into the package.
Chunk over test points (say 2000 at a time) and do it in one shot; also, since
`weight.exponent` defaults to `-1`, special-case it to `1/d` rather than `d^-1` (`pow` is
~10× a divide).

The `parallel = TRUE` path spins up a cluster to parallelise these tiny calls, which mostly
buys back the overhead it's paying. Vectorise first, then decide if parallelism is still needed.

---

## 5. `src/convert.cpp` + `R/kdtree_build.R` — the data is copied three times

`kdtree_build()` does `t(data)` (copy 1) → `convertMatrixToVector()` builds a
`vector<vector<double>>` (copy 2, with a `/` and `%` per element) → the `KDTree` constructor
does `this->points = points` (copy 3). For 10⁶ × 5 doubles that's ~120 MB of avoidable
traffic, and `vector<vector<double>>` scatters every point onto its own heap allocation,
which is why the ball queries are cache-hostile.

Fixes, in increasing order of effort:

- `KDTree(vector<vector<double>>&& points) : points(std::move(points))` — kills copy 3 for free.
- Drop the `t()` in R and index the column-major buffer directly in C++ — kills copy 1.
- Store points as one flat `vector<double>` with stride `ndim` and index `pts[i*ndim + j]`.
  Also replaces the `vector<Node*>` + per-node `new` with a `vector<Node>` you can `reserve()`.
  This is the change that makes queries fast, not just build.

Also: the tree currently has **one leaf per point**, so a dataset of 10⁶ points allocates
~2×10⁶ individually-`new`'d 24-byte nodes. Bucketing leaves at ~16–32 points and
brute-forcing inside the bucket is the standard fix and typically gives 2–4× on queries
plus a large drop in build time and memory.

---

## 6. Free parallelism in the ball queries

`src/Makevars` is just:

```make
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

No OpenMP. The loops in `kdtree_ball_query_multiple`, `kdtree_ball_query_id_multiple` and
`kdtree_range_query_multiple` are embarrassingly parallel — each query point is independent
and the tree is read-only. Add:

```make
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS     = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

and `#pragma omp parallel for` over `i`, writing into a pre-sized `finalCounts` (currently
`push_back` into an unreserved vector) instead of appending. Near-linear in cores. Keep it
off by default on CRAN-checked builds if you'd rather not deal with the thread-limit policy —
`omp_set_num_threads(...)` from an R option is the usual dodge.

While in there: `kdtree_ball_query_multiple` only needs a **count**, but `ball_query()`
builds and fills `vector<int> thisIndices` and `vector<double> thisDistances` for every query.
A count-only overload that increments an integer removes two heap allocations per query point
— and this is the single hottest call in the package.
Also `vector<double> thisPoint = dataMatrix[i];` copies each query point; use a reference.

---

## 7. `R/samplers.R` — small stuff, real money

In `sample_model_ellipsoid()`:

- `ellipsoid_inverse_weight()` rebuilds the kd-tree over `centers` on **every iteration of the
  `while (total_accepted < N.samples)` loop**, but `centers` (= your data) never changes.
  Hoist the build out of the loop. Combined with fix #1 this is a large win on its own.
- `full_samples = lapply(1:nrow(data), function(i) data.frame(sample_ellipsoid(...)))` builds
  one `data.frame` per data point just to feed `rbindlist`. `data.frame()` is one of the more
  expensive constructors in R. Pre-allocate a single `nrow(data)*chunk.size × d` matrix and
  fill it, or at minimum keep them as matrices and `do.call(rbind, ...)`.
- `samples <- c(samples, list(data.frame(samples_final_this)))` — same conversion again per
  iteration.

In `sample_model_rejection()`: `random_points` is built as a `data.frame` and then rescaled
column-by-column in a loop; a matrix + `sweep()` is cheaper, though `predict()` methods may
want the data.frame — worth checking per-model.

---

## 8. `R/hypervolume_set.R` — `unique()` on data.frames

```r
final_points_intersection = unique(rbind(p1_in_2, p2_in_1))
final_points_union        = unique(rbind(p1_not_in_2, p2_not_in_1, p_in_1_and_2))
```

`unique.data.frame` pastes rows into strings and calls `duplicated` on them. On 10⁵–10⁶ rows
of doubles this is brutal and allocates a character vector the size of the dataset. Since
everything here starts as a numeric matrix, stay in matrix-land and use
`m[!duplicated(m), , drop=FALSE]` (which dispatches to `duplicated.matrix` / `duplicated.array`)
— or, better, note that these point clouds come from disjoint samples, so the `unique()` calls
may be defensive rather than necessary.

Related: `matrix(as.matrix(as.data.frame(final_points_union)), ncol=dim)` appears four times —
three conversions where zero are needed once the inputs stay matrices.

The two `evalfspherical()` calls (`p2_in_1_all`, `p1_in_2_all`) are independent and each
builds its own tree — an easy `foreach` pair if you want the 2× (`doParallel` is already a
dependency).

---

## 9. `src/fastPdist2.cpp` — correctness, not speed

```cpp
arma::mat C = -2 * (A * B.t());
C.each_col() += An;
C.each_row() += Bn.t();
return wrap(sqrt(C));
```

The `|a|² + |b|² − 2a·b` identity is fine for speed but is numerically unstable for nearby
points: `C` can come out slightly negative, and `sqrt` of that gives `NaN`. `hypervolume_distance(type="minimum")`
then does `min(..., na.rm=TRUE)`, which silently hides it — but the true minimum distance is
exactly the case most at risk. Clamp before the sqrt:

```cpp
C.transform([](double v){ return v < 0.0 ? 0.0 : v; });
```

Also worth exposing a squared-distance variant (no `sqrt` over m×n elements) since most
callers — including the `calculate_density` rewrite in #2 — square it right back.

---

## Suggested order

1. **#1** (kd-tree `workarray`) — two lines, no API change, biggest single win. Do this first
   and re-profile; it may change the ranking of everything below it.
2. **#3** (threshold) — self-contained, removes the memory ceiling on `hypervolume_gaussian`.
3. **#2** (`calculate_density`) — biggest win on the flagship function, but needs a numerical
   regression test against the current output before you trust it.
4. **#7** (hoist the tree build out of the sampler loop) — one line.
5. **#6** (OpenMP) and **#5** (flat storage / bucketed leaves) — the "make queries genuinely
   fast" work, more involved.
6. **#4**, **#8**, **#9** — worthwhile, lower stakes.

### Before you change anything

There's no test suite or benchmark harness in the repo. Given #2 and #3 change numerical
output paths, I'd pin down expected volumes for the `penguins` / `quercus` demo datasets
under a fixed seed first, so the refactors are verifiable rather than hopeful.
