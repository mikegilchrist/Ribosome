# A.7 figure regression: the three Fig 10 curves (Rmd: blocks 45001:50000,
# 40001:45000, 20001:25000 -> mu 2.2216e-4, 3.3007e-4, 2.3105e-3), imax=39,
# tau=lambda=1, delta=1e5, which.interfere=c("kappa","tau").
# Protein production = colSums(classes) * Ribo_mean(classes) (Rmd convention).
# Production rises monotonically in kappa, so the peak is at kappa = 0.5, the
# last point of the 25-point grid; no finer grid is needed.
fig10_mu <- c(0.0002221626, 0.0003300701, 0.0023104906)

test_that("Fig 10 normalised peaks 1.000, 0.692, 0.143 and decapped shares 0.072, 0.104, 0.448", {
  nk <- length(KAPPA_GRID)
  g <- make_grid(fig10_mu, KAPPA_GRID)
  wi <- c("kappa", "tau")
  C <- solve_grid(g, 39, wi, "classes"); U <- solve_grid(g, 39, wi, "unmarked")
  prod <- colSums(C) * as.numeric(Ribo_mean(C))
  blk <- rep(1:3, each = nk)
  peaks <- tapply(prod, blk, max)
  expect_equal(as.numeric(round(peaks / max(peaks), 3)), c(1.000, 0.692, 0.143))
  expect_equal(as.numeric(tapply(prod, blk, which.max)), rep(nk, 3))  # peak at kappa=0.5
  M <- C - U
  pU <- colSums((0:39) * U); pM <- colSums((0:39) * M)
  share <- pM / (pU + pM)
  for (b in 1:3) expect_equal(range(share[blk == b]), rep(c(0.072, 0.104, 0.448)[b], 2), tolerance = 1e-8 + 0.0006)
  expect_equal(as.numeric(round(tapply(share, blk, mean), 3)), c(0.072, 0.104, 0.448))
  # production equals total bound ribosomes
  expect_equal(prod, colSums((0:39) * C), tolerance = 1e-10)
})

test_that("Fig 10 peaks unchanged on the full Rmd kappa grid seq(0.0001,0.5,0.0001)", {
  g <- make_grid(fig10_mu, seq(0.0001, 0.5, 0.0001))
  C <- solve_grid(g, 39, c("kappa", "tau"), "classes")
  prod <- colSums(C) * as.numeric(Ribo_mean(C))
  peaks <- tapply(prod, rep(1:3, each = 5000), max)
  expect_equal(as.numeric(round(peaks / max(peaks), 3)), c(1.000, 0.692, 0.143))
})
