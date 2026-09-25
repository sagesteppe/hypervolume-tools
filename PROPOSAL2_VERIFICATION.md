# Verifying proposal #2 against the published method

Source: Blonder, Lamanna, Violle & Enquist (2014), "The n-dimensional
hypervolume", *Global Ecology and Biogeography* 23:595–609. The mathematics is
in **Box 2, p. 599**.

**Caveat on scope, stated up front.** This paper describes the *box-kernel*
method (`hypervolume_box`, and the v1 package). `hypervolume_gaussian` and
`calculate_density` come from the later Blonder et al. (2018) *Methods in
Ecology and Evolution* paper, "New approaches for delineating n-dimensional
hypervolumes", which is not in the folder. So this document verifies the parts
of proposal #2 that rest on the shared estimator scaffolding — the normalising
constant, the thresholding rule, and the importance-sampling correction — all
of which carry over unchanged. The specifically Gaussian choices (the `sd.count`
truncation radius, the diagonal bandwidth matrix) are *not* adjudicated here.
Getting the 2018 paper would close that gap.

---

## 1. The `weight` note is confirmed. There is one 1/m too many.

Box 2 defines the kernel density estimate as

> Ẑ(pᵢ) = α · Σⱼ₌₁..m k((pᵢ − wⱼ)/h)

and fixes the normalising constant as

> α = 1 / (m · ∏ᵢ₌₁..ⁿ h·eᵢ)

So the published estimator contains **exactly one factor of 1/m**, where m is
the number of observations. The `∏ h` part is the per-kernel normalisation.

Now compare `R/hypervolume_gaussian.R`:

```r
for (i in index_vals) {
    density = density + weight[i] * dmvnorm(s, mean = means[i, ], sigma = kernel_sd^2 * diag(d))
}
density = density / n_points
```

`dmvnorm` is already a normalised density, so it supplies the `k(...)/∏h` part
itself. That leaves the 1/m to come from somewhere, and it arrives **twice**:

- `weight` is normalised by the caller to sum to 1 — for the default
  `rep(1/nrow(data), nrow(data))`, that is a factor of 1/m; and then
- `density / n_points` divides by m again.

Result: `Ẑ_code = Ẑ_paper / m`. The tests pin this as
`calculate_density` integrating to `1/n` rather than `1`
(`tests/testthat/test-04-calculate-density.R`).

### Does it matter?

For volume, no. Two independent scale-invariances protect it:

1. **The threshold rule is a ratio.** Box 2 flags indices where Ẑ(pᵢ) > q for
   the q satisfying `Σ_flagged Ẑ(pᵢ) ≥ τ · Σ_all Ẑ(pᵢ)`. Both sides scale
   together, so q moves with the scale and the same points are selected.
   `hypervolume_threshold` implements this as
   `probability_sums / max(probability_sums)`, equally scale-free.
2. **Volume never reads the density magnitude.** It comes out of
   `sample_model_ellipsoid` as `ellipsoid_volume × n × mean_weight ×
   total_accepted/total_tried`. The only place the density enters is the
   acceptance test `predicted_values > min.value` with `min.value = 0`, and a
   Gaussian density is strictly positive, so acceptance is unconditional and
   scale cannot reach it.

For `ValueAtRandomPoints`, yes — mildly. The slot is off by 1/m, so it is not a
probability density in the units a reader would assume, and it is **not
comparable across hypervolumes built from different numbers of observations**.
Anything that reads it in absolute terms inherits that:
`hypervolume_estimate_probability` (a weighted mean of `prob_ss`, so linear in
the scale), `occupancy_bootstrap_gof`, and the plot colour scales.

### Recommendation

**Do not bundle this into proposal #2.** It is a pre-existing correctness
question, not a performance one, and fixing it silently inside a "pure speedup"
commit is how a refactor gets blamed for a numbers change six months later.
Land the BLAS rewrite reproducing the current scaling exactly — the tests
enforce that — then fix the normalisation as its own commit with its own
before/after on the functions listed above.

---

## 2. The BLAS rewrite itself is safe

Box 2's estimator is a plain weighted sum of kernels over the data points. With
a diagonal bandwidth the Gaussian kernel factorises across dimensions, so

    Ẑ(p) ∝ Σⱼ wⱼ · exp(−½ Σᵢ ((pᵢ − wⱼᵢ)/hᵢ)²)

is a matrix of scaled squared distances, exponentiated, multiplied by a weight
vector. That is exactly what the rewrite computes. Nothing in the method
constrains the *order* of the sum, which is the only thing that changes.

Two caveats the tests already cover:

- **Return shape.** `exp(-0.5*D2) %*% weight` is an n×1 matrix, not a vector.
  `sample_model_ellipsoid` does `as.numeric(predicted_values) > min.value` and
  then row-indexes with the result, so a stray `dim` attribute survives several
  lines before failing somewhere unhelpful. Pinned by "returns a bare numeric
  vector, not a matrix".
- **The normalising constant.** Pinned by the scale-equivariance test in
  `test-07`: multiplying data and bandwidth by c must scale volume by exactly
  `c^d`. Because acceptance is deterministic and thresholding is quantile-based,
  that holds to floating point, not just Monte Carlo precision — so a wrong
  constant shows up as a hard failure rather than a noisy one.

---

## 3. Proposal #7 is not just safe, it is what the paper describes

Box 2, on the importance-sampling correction:

> "We build an n-dimensional recursive partitioning tree T **from the data
> points w**. Then for each pᵢ, we perform a range-box query on T […] and count
> the number of non-zero returns, which is proportional to Ẑ(pᵢ)."

The tree is built once from the observations and queried per sample point. The
current code rebuilds it on every iteration of the `while (total_accepted <
N.samples)` loop in `sample_model_ellipsoid`, even though `centers` is the
(fixed) input data. Hoisting the build out of the loop restores the published
algorithm rather than departing from it.

The oversampling correction itself — `1/Ẑ(pᵢ)` weighting, effective sample size
`ρ = Σ 1/Ẑ(pᵢ)`, volume `π/ρ` — maps onto `ellipsoid_inverse_weight` and the
`mean_weight` arithmetic. `test-06-samplers.R` pins it on three cases where the
answer is analytic rather than Monte Carlo: one ellipsoid returns exactly
`ellipsoid_volume(scales)`; two coincident ellipsoids return the same, not
double; two far-apart ellipsoids return exactly twice.

---

## 4. Proposals #1, #5 and #6 are unconstrained by the method

The tree is a spatial index. Box 2 is explicit that the box kernel was chosen so
that "the evaluation of the kernel density estimate [is] reduced to a counting
problem" — the method specifies *what* is counted, never how the index is
built. Fixing the O(n²) `workarray` clear, flattening point storage, bucketing
leaves, and parallelising the query loop all leave query results bit-identical.
`test-02-kdtree.R` checks them against brute force, including the degenerate
inputs (duplicates, collinear points, n = 2) that a rewritten split can break.

---

## 5. One proposal I now think is NOT safe: neighbour truncation

In the notes I floated restricting the density sum to kernels within
`sd.count × bandwidth`, on the grounds that sampling is already truncated there.
Reading Box 2 makes me less comfortable with it.

The box kernel was chosen specifically because it "reaches zero in a finite
distance" — the truncation is a *property of the kernel*, so restricting the sum
is exact. A Gaussian kernel has infinite support, so the same restriction is an
approximation, and it is not the same approximation as truncating the sampling
region. A sample point near the edge of observation A's ellipsoid still receives
genuine contributions from observations B, C… further away; dropping them
changes Ẑ, which changes the τ-quantile threshold, which changes which points
are retained, which changes the volume.

The error is probably small — beyond 3 SD a Gaussian contributes ~1e-4 of its
peak — but "probably small" is a claim that needs measuring, not asserting.

**Reclassify it**: not low-hanging fruit, but a method change needing a
before/after volume comparison across bandwidths and dimensionalities, and
ideally a word from Blonder. The BLAS rewrite alone should give 1–2 orders of
magnitude, which is likely enough that the truncation is not worth the argument.

---

## Summary

| Proposal | Verdict against the paper |
|---|---|
| #2 BLAS rewrite of `calculate_density` | **Safe.** Same estimator, different summation order. |
| #2 `weight` / `1/n` normalisation | **Confirmed off by 1/m** (α = 1/(m·∏h), Box 2). Affects `ValueAtRandomPoints` only, not volume. Fix separately, not inside the perf commit. |
| #2 neighbour truncation | **Downgraded.** Exact for a box kernel, approximate for a Gaussian. Needs measurement first. |
| #7 hoist the tree build | **Safe, and restores the published algorithm.** |
| #1, #5, #6 kd-tree internals | **Unconstrained by the method.** Pure implementation. |
| #3 threshold rewrite | **Safe.** Box 2's rule is a ratio of sums; cumulative sums compute it exactly. |
| #9 `fastPdist2` clamp | **Safe.** Repairs a numerical defect; no method content. |

**Open item:** get Blonder et al. (2018) *MEE* for the Gaussian-specific
choices. Nothing above depends on it, but the `sd.count` truncation radius and
the diagonal-bandwidth assumption are defined there, not here.
