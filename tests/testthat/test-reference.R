# Package solvers (first set) vs the independent base-R oracle, plus a check
# that the helper's generalised reference reduces to the oracle itself.

test_that("helper ref_solve (kappa interference) equals the oracle solve_model", {
  for (r in c(1, 77, 250)) {
    row <- make_grid()[r, ]
    a <- solve_model(row$kappa, row$mu, row$tau, 39)
    b <- ref_solve(row$kappa, row$mu, row$tau, 39)
    expect_equal(b$m, a$m, tolerance = 1e-12)
    expect_equal(b$mstar, a$mstar, tolerance = 1e-12)
  }
})

for (nm in names(INTERFERE_SETS)) for (im in IMAX_SET) local({
  nm <- nm; im <- im; wi <- INTERFERE_SETS[[nm]]
  test_that(sprintf("CalcUnmarked/Marked/Classes match reference: imax=%d, which.interfere=%s",
                    im, paste(wi, collapse = "+")), {
    g <- make_grid()
    U <- solve_grid(g, im, wi, "unmarked"); M <- solve_grid(g, im, wi, "marked")
    C <- solve_grid(g, im, wi, "classes")
    for (r in seq_len(nrow(g))) {
      ref <- ref_for(g[r, ], im, wi)
      expect_equal(U[, r], ref$m, tolerance = 1e-8)
      expect_equal(M[, r], ref$mstar, tolerance = 1e-8)
      expect_equal(C[, r], ref$m + ref$mstar, tolerance = 1e-8)
    }
  })
})
