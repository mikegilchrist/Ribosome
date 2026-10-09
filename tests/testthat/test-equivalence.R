# A.5 equivalence between the duplicate solver sets.
# Result (measured, see KNOWN_FAILURES.md):
#  * _newtau set  : AGREED with the first set / reference (max rel diff ~1e-13)
#                   at 03be8d5; removed from the package in Part C.3, so its
#                   tests were removed with it (see KNOWN_FAILURES.md #1).
#  * *NP set      : does NOT agree at the same input vector, because it uses
#                   tau_i = i (rescaled) instead of tau*i/(9*imax). It DOES equal
#                   the first set (kappa+tau interference) evaluated at
#                   tau_first = 9*imax*tau_NP; that relation is tested below.

test_that("NP set equals first set at tau_first = 9*imax*tau_NP (kappa+tau interference)", {
  g <- make_grid()
  for (im in IMAX_SET) for (r in c(1, 100, 250)) {
    p <- par_vec(g[r, ]); p2 <- p; p2["tau"] <- 9 * im * p["tau"]
    expect_equal(CalcUnmarkedClassNP(p, im),
                 CalcUnmarkedClass(p2, im, which.interfere = c("kappa", "tau")), tolerance = 1e-10)
    expect_equal(CalcMarkedClassNP(p, im),
                 CalcMarkedClass(p2, im, which.interfere = c("kappa", "tau")), tolerance = 1e-10)
    expect_equal(ClassesNP(p, im),
                 Classes(p2, im, which.interfere = c("kappa", "tau")), tolerance = 1e-10)
  }
})

test_that("NP set equals the reference at the SAME input vector", {
  skip("KNOWN FAILURE (KNOWN_FAILURES.md #2): NP uses tau_i=i, not tau*i/(9*imax); differs from reference by orders of magnitude")
  p <- par_vec(make_grid()[100, ]); ref <- ref_for(make_grid()[100, ], 39, c("kappa", "tau"))
  expect_equal(CalcUnmarkedClassNP(p, 39), ref$m, tolerance = 1e-8)
  expect_equal(CalcMarkedClassNP(p, 39), ref$mstar, tolerance = 1e-8)
})

test_that("NP set equals the first set at the SAME input vector", {
  skip("KNOWN FAILURE (KNOWN_FAILURES.md #2): same cause as NP-vs-reference")
  p <- par_vec(make_grid()[100, ])
  expect_equal(ClassesNP(p, 39), Classes(p, 39), tolerance = 1e-8)
})
