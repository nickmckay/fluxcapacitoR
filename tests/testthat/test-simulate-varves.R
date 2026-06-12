test_that("simulateVarves AR(1) returns positive thicknesses with correct dims", {
  set.seed(42)
  v <- simulateVarves(n = 50, n.ens = 5, ar1 = 0.5, mean = 2)
  expect_equal(dim(v), c(50, 5))
  expect_true(all(v > 0))
  expect_equal(mean(v), 2, tolerance = 0.1)
})

test_that("simulateVarves long-memory (H) path works", {
  set.seed(42)
  v <- simulateVarves(n = 30, n.ens = 2, H = 0.8)
  expect_equal(dim(v), c(30, 2))
  expect_true(all(v > 0))
})

test_that("fgnEnsemble matches fractional Gaussian noise statistics", {
  set.seed(42)
  H <- 0.8
  x <- fgnEnsemble(n = 2000, n.ens = 200, H = H)
  expect_equal(dim(x), c(2000, 200))

  # unit variance and mean ~ 0
  expect_equal(mean(apply(x, 2, var)), 1, tolerance = 0.05)
  expect_equal(mean(colMeans(x)), 0, tolerance = 0.1)

  # empirical autocorrelation tracks the theoretical fGn function. The sample
  # ACF is biased slightly low for long memory (the sample mean absorbs
  # low-frequency variance), so the tolerance allows for that downward bias --
  # it still easily separates fGn from white noise or a different H.
  rho <- function(k, H) 0.5 * (abs(k + 1)^(2*H) - 2*abs(k)^(2*H) + abs(k - 1)^(2*H))
  empAcf <- rowMeans(apply(x, 2, \(col) stats::acf(col, lag.max = 3, plot = FALSE)$acf[2:4]))
  expect_true(all(abs(empAcf - rho(1:3, H)) < 0.06))
  expect_true(all(empAcf > 0.2))

  # stronger persistence at higher H
  set.seed(7)
  lowH  <- mean(apply(fgnEnsemble(2000, 100, 0.55), 2, \(c) stats::acf(c, 1, plot = FALSE)$acf[2]))
  highH <- mean(apply(fgnEnsemble(2000, 100, 0.95), 2, \(c) stats::acf(c, 1, plot = FALSE)$acf[2]))
  expect_gt(highH, lowH)
})

test_that("simulateVarves pads with NA to length.out", {
  set.seed(42)
  v <- simulateVarves(n = 20, n.ens = 3, ar1 = 0.3, length.out = 25)
  expect_equal(dim(v), c(25, 3))
  expect_true(all(is.na(v[21:25, ])))
  expect_true(all(is.finite(v[1:20, ])))
})

test_that("simulateVarves requires exactly one of ar1 and H", {
  expect_error(simulateVarves(n = 10, ar1 = 0.5, H = 0.8), "not both")
  expect_error(simulateVarves(n = 10), "must specify")
})

test_that("gammify preserves rank order and matches target mean", {
  set.seed(11)
  x <- matrix(rnorm(200), ncol = 2)
  g <- gammify(x, shape = 2, mean = 5)
  expect_equal(dim(g), dim(x))
  expect_true(all(g > 0))
  expect_equal(order(x[, 1]), order(g[, 1]))
  expect_equal(mean(g), 5, tolerance = 0.5)
})
