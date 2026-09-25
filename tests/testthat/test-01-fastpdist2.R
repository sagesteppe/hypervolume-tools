# Covers proposal #9 (fastPdist2 numerical stability).
#
# fastPdist2 uses the |a|^2 + |b|^2 - 2a.b expansion, which is fast but can
# produce small negative values for near-coincident points, and sqrt() of those
# is NaN. These tests pin the current behaviour so the clamp can be added
# without changing anything else.

test_that("fastPdist2 matches a naive distance calculation", {
  A <- toy_data(12L, 3L, seed = HV_SEED)
  B <- toy_data(9L, 3L, seed = HV_SEED + 1L)

  observed <- hypervolume:::fastPdist2(A, B)
  expected <- ref_pdist(A, B)

  expect_equal(dim(observed), c(nrow(A), nrow(B)))
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("fastPdist2 handles single-row inputs on both sides", {
  A <- toy_data(1L, 4L, seed = HV_SEED)
  B <- toy_data(6L, 4L, seed = HV_SEED + 2L)

  expect_equal(dim(hypervolume:::fastPdist2(A, B)), c(1L, 6L))
  expect_equal(dim(hypervolume:::fastPdist2(B, A)), c(6L, 1L))
  expect_equal(
    as.numeric(hypervolume:::fastPdist2(A, B)),
    as.numeric(t(hypervolume:::fastPdist2(B, A))),
    tolerance = 1e-8
  )
})

test_that("fastPdist2 is symmetric and zero on the diagonal for A vs A", {
  A <- toy_data(10L, 2L, seed = HV_SEED)
  D <- hypervolume:::fastPdist2(A, A)

  expect_equal(as.numeric(D), as.numeric(t(D)), tolerance = 1e-8)
  # The diagonal is the exact case where catastrophic cancellation bites.
  # It should be 0, and it must not be NaN.
  expect_false(anyNA(diag(D)))
  expect_equal(diag(D), rep(0, nrow(A)), tolerance = 1e-6)
})

test_that("fastPdist2 does not return NaN for near-coincident points", {
  # Points separated by far less than sqrt(.Machine$double.eps) relative to
  # their magnitude: |a|^2 + |b|^2 - 2a.b underflows to a negative number here.
  base <- matrix(rep(1e6, 3), nrow = 1)
  A <- rbind(base, base + 1e-9, base - 1e-9)
  colnames(A) <- paste0("X", 1:3)

  D <- hypervolume:::fastPdist2(A, A)

  # Documents the current contract. If this starts failing after the clamp is
  # added, that is the fix landing, not a regression -- delete the skip below.
  expect_false(anyNA(D))
  expect_true(all(D >= 0 | is.na(D)))
})

test_that("fastPdist2 preserves scale across dimensionalities", {
  for (d in c(1L, 2L, 5L, 8L)) {
    A <- toy_data(7L, d, seed = HV_SEED + d)
    B <- toy_data(5L, d, seed = HV_SEED + 100L + d)
    expect_equal(
      as.numeric(hypervolume:::fastPdist2(A, B)),
      as.numeric(ref_pdist(A, B)),
      tolerance = 1e-8,
      info = paste("d =", d)
    )
  }
})
