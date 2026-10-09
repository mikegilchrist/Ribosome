# A.4 invariants. Each was confirmed on the reference solver first (max
# violation ~1e-13 or smaller on every grid point). Tested for imax in
# {4, 39, 194}, both which.interfere settings, all 10 mu x 25 kappa points
# (kappa 1e-4 .. 0.5 log grid), unless noted. Class i is at R index i+1.

for (nm in names(INTERFERE_SETS)) for (im in IMAX_SET) local({
  nm <- nm; im <- im; wi <- INTERFERE_SETS[[nm]]
  lab <- sprintf("imax=%d, which.interfere=%s", im, paste(wi, collapse = "+"))
  g <- make_grid()
  U <- solve_grid(g, im, wi, "unmarked"); M <- solve_grid(g, im, wi, "marked")
  C <- solve_grid(g, im, wi, "classes")
  cls <- 0:im

  test_that(paste("capped total sum(unmarked)==lambda/mu;", lab), {
    expect_equal(colSums(U), g$lambda / g$mu, tolerance = 1e-9)
  })
  test_that(paste("decapped class 0 marked[1]==lambda/delta;", lab), {
    expect_equal(M[1, ], g$lambda / g$delta, tolerance = 1e-9)
  })
  test_that(paste("decapped recurrence;", lab), {
    for (r in seq_len(nrow(g))) {
      i <- 1:im
      lhs <- g$tau[r] * i / (9 * im) * M[i + 1, r]
      rhs <- g$mu[r] * sapply(i, function(k) sum(U[(k + 1):(im + 1), r]))
      expect_equal(lhs, rhs, tolerance = 1e-9)
    }
  })
  test_that(paste("steady-state residual A u + b == 0 (reference tridiagonal);", lab), {
    for (r in seq_len(nrow(g))) {
      ref <- ref_for(g[r, ], im, wi)
      expect_lt(max(abs(ref$A %*% U[, r] + ref$b)), 1e-9)
    }
  })
  test_that(paste("probabilities: Ribo_percent sums to 1 (vector and columns);", lab), {
    expect_equal(sum(Ribo_percent(C[, 1])), 1, tolerance = 1e-12)
    expect_equal(colSums(Ribo_percent(C)), rep(1, ncol(C)), tolerance = 1e-12)
  })
  test_that(paste("decapped share of protein output identity;", lab), {
    pU <- colSums(cls * U); pM <- colSums(cls * M)
    rr <- 9 * g$mu * im / g$tau
    expect_equal(pM / (pU + pM), rr / (1 + rr), tolerance = 1e-8)
  })
  test_that(paste("decapped marked[-1] non-increasing in i;", lab), {
    expect_true(all(apply(M[-1, , drop = FALSE], 2, function(m) all(diff(m) <= 0))))
  })
  test_that(paste("non-negativity of all classes;", lab), {
    expect_true(all(U >= 0)); expect_true(all(M >= 0)); expect_true(all(C >= 0))
  })
})
