# CHECK_STEADY_STATE (line 158) was examined as a possible reusable residual
# test (brief A.4). Findings (KNOWN_FAILURES.md #4, #5):
#  - it uses tau_i = i (the NP convention), so it FAILS on first-set output;
#  - its `return` sits inside the for loop, so only class 0 is ever checked.
p <- par_vec(make_grid()[100, ])

test_that("CHECK_STEADY_STATE on first-set output (should pass; does not)", {
  skip("KNOWN FAILURE (KNOWN_FAILURES.md #4): CHECK_STEADY_STATE assumes tau_i=i, not tau*i/(9*imax)")
  U <- CalcUnmarkedClass(p, 39); M <- CalcMarkedClass(p, 39, U)
  expect_equal(CHECK_STEADY_STATE(U, M, p, verbose = FALSE), c(0, 0))
})
test_that("CHECK_STEADY_STATE currently flags first-set output as failing (negative bookkept)", {
  U <- CalcUnmarkedClass(p, 39); M <- CalcMarkedClass(p, 39, U)
  expect_equal(CHECK_STEADY_STATE(U, M, p, verbose = FALSE), c(1, 1))
})
test_that("CHECK_STEADY_STATE blind spot: only class 0 is checked (corrupt class 5, still passes)", {
  U <- CalcUnmarkedClassNP(p, 39); M <- CalcMarkedClassNP(p, 39, U)  # NP matches its tau_i=i convention
  expect_equal(CHECK_STEADY_STATE(U, M, p, verbose = FALSE), c(0, 0))
  U2 <- U; U2[6] <- U2[6] * 10
  expect_equal(CHECK_STEADY_STATE(U2, M, p, verbose = FALSE), c(0, 0))
  res <- CHECK_STEADY_STATE(U2, M, p, verbose = FALSE, post = TRUE)
  expect_true(all(res[[2]][, -1] == 0))
})
