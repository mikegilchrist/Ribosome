# A.6 Ribo_* summary statistics against hand-computed values.
# Toy input: classes 0..3 with weights 1,2,3,4 -> p = .1,.2,.3,.4;
#   mean = .2+.6+1.2 = 2 ; var = .1*4+.2*1+.3*0+.4*1 = 1 -> sd 1 ; mode = class 3.
# Second column (reverse): p = .4,.3,.2,.1 ; mean = 1 ; var = .4+0+.2+.4 = 1 ; mode = class 0.
x <- c(1, 2, 3, 4)
M <- cbind(a = x, b = rev(x))

test_that("Ribo_percent: vector and matrix", {
  expect_equal(Ribo_percent(x), c(.1, .2, .3, .4))
  expect_equal(unname(Ribo_percent(M)), cbind(c(.1, .2, .3, .4), c(.4, .3, .2, .1)))
})
test_that("Ribo_mean: vector (1x1 matrix) and matrix (ncol x 1)", {
  expect_equal(as.numeric(Ribo_mean(x)), 2)
  expect_equal(dim(Ribo_mean(x)), c(1, 1))
  expect_equal(as.numeric(Ribo_mean(M)), c(2, 1))
  expect_equal(dim(Ribo_mean(M)), c(2, 1))
})
test_that("Ribo_mode: returns class index (0-based)", {
  expect_equal(as.numeric(Ribo_mode(x)), 3)
  expect_equal(as.numeric(Ribo_mode(M)), c(3, 0))
})
test_that("Ribo_sd: matrix branch", {
  expect_equal(as.numeric(Ribo_sd(M)), c(1, 1))
})
test_that("Ribo_sd: vector branch (fixed in Part C.5, KNOWN_FAILURES.md #3)", {
  expect_equal(as.numeric(Ribo_sd(x)), 1)
  expect_equal(as.numeric(Ribo_sd(rev(x))), 1)
  # vector and one-column matrix input agree
  expect_equal(as.numeric(Ribo_sd(x)), as.numeric(Ribo_sd(cbind(x))))
})
test_that("Ribo_frac and Ribo_frac_load", {
  expect_equal(Ribo_frac(M, 2 * M), c(.5, .5))
  out <- capture.output(r <- Ribo_frac_load(c(.5, .5), c(1, 2), c(4, 8)))  # prints length(Denominator)
  expect_equal(r, c(.125, .125))
})
test_that("Ribo_* on a scalar are not supported (apply on non-matrix errors)", {
  expect_error(Ribo_percent(1))
})
