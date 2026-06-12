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
