# A.5 equivalence between the duplicate solver sets.
# Result (measured, see KNOWN_FAILURES.md):
#  * _newtau set  : AGREES with the first set / reference (max rel diff ~1e-14).
#                   It is textually identical apart from print() calls, and
#                   CalcMarkedClass_newtau falls back to CalcUnmarkedClass.
#  * *NP set      : does NOT agree at the same input vector, because it uses
#                   tau_i = i (rescaled) instead of tau*i/(9*imax). It DOES equal
#                   the first set (kappa+tau interference) evaluated at
#                   tau_first = 9*imax*tau_NP; that relation is tested below.

test_that("_newtau set equals first set and reference (both interference settings, 3 imax)", {
  g <- make_grid()
  rows <- c(1, 25, 26, 100, 149, 250)
  for (nm in names(INTERFERE_SETS)) for (im in IMAX_SET) {
    wi <- INTERFERE_SETS[[nm]]
    for (r in rows) {
      p <- par_vec(g[r, ]); ref <- ref_for(g[r, ], im, wi)
      capture.output({
        nu <- CalcUnmarkedClass_newtau(p, im, which.interfere = wi)
        nm_ <- CalcMarkedClass_newtau(p, im, which.interfere = wi)
        nm2 <- CalcMarkedClass_newtau(p, im, unmarked = nu, which.interfere = wi)
        nc <- Classes_newtau(p, im, which.interfere = wi)
      })
      # outputs carry parameter names as attributes in some paths; compare values
      nu <- as.numeric(nu); nm_ <- as.numeric(nm_); nm2 <- as.numeric(nm2); nc <- as.numeric(nc)
      expect_equal(nu, as.numeric(CalcUnmarkedClass(p, im, which.interfere = wi)), tolerance = 1e-10)
      expect_equal(nm_, as.numeric(CalcMarkedClass(p, im, which.interfere = wi)), tolerance = 1e-10)
      expect_equal(nm2, nm_, tolerance = 1e-10)
      expect_equal(nc, as.numeric(Classes(p, im, which.interfere = wi)), tolerance = 1e-10)
      expect_equal(nu, ref$m, tolerance = 1e-8)
      expect_equal(nm_, ref$mstar, tolerance = 1e-8)
    }
  }
})

test_that("_newtau functions print to stdout (debug calls still present)", {
  p <- par_vec(make_grid()[100, ])
  out <- capture.output(CalcUnmarkedClass_newtau(p, 4))
  expect_gt(length(out), 0)
})

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
