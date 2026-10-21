# Unit tests for argument validation and informative messages.

test_that("an unknown wavelet family raises an error", {
  x <- sort(runif(16))
  expect_error(
    PHI(x, J = 3, family = "not-a-family"),
    "Unknown family"
  )
})

test_that("a negative coarsest level is rejected", {
  x <- sort(runif(16))
  expect_error(
    wavedec(x, j0 = -1, family = "daublets", filter.size = 8),
    "non-negative"
  )
})

test_that("family 'Own' without a filter is rejected", {
  x <- sort(runif(16))
  expect_error(
    PHI(x, J = 3, family = "Own"),
    "Provide your own filter"
  )
})

test_that("complex input is coerced to its real part with a warning", {
  set.seed(123)
  x <- complex(real = sort(runif(8)), imaginary = rnorm(8))
  expect_warning(
    wavedec(x, j0 = 0, family = "daublets", filter.size = 8),
    "complex"
  )
})

test_that("wavedec and waverec reject a matrix, an array and a data frame", {
  set.seed(4)
  M <- matrix(rnorm(8 * 16), 8, 16)
  expect_error(wavedec(M), "not a matrix")
  expect_error(waverec(M), "not a matrix")
  expect_error(
    wavedec(M, j0 = 4, family = "d", filter.size = 8, boundary = "interval"),
    "not a matrix"
  )
  expect_error(wavedec(array(rnorm(64), c(4, 4, 4))), "not a matrix")
  expect_error(wavedec(data.frame(a = rnorm(8), b = rnorm(8))),
               "numeric vector")
  expect_error(waverec(data.frame(a = rnorm(8))), "numeric vector")

  # a single row or column is still a vector
  x <- rnorm(16)
  expect_equal(wavedec(matrix(x, ncol = 1)), wavedec(x))
  expect_equal(wavedec(matrix(x, nrow = 1)), wavedec(x))
  expect_equal(waverec(matrix(x, ncol = 1)), waverec(x))
})
