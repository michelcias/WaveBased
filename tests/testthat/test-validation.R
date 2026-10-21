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

test_that("wavedec and waverec reject an array and a data frame", {
  set.seed(4)
  A <- array(rnorm(64), c(4, 4, 4))
  expect_error(wavedec(A), "not an array")
  expect_error(waverec(A), "not an array")
  expect_error(wavedec(data.frame(a = rnorm(8), b = rnorm(8))),
               "numeric vector or matrix")
  expect_error(waverec(data.frame(a = rnorm(8))), "numeric vector or matrix")

  # a matrix with a single column, or an array with a single dimension of
  # extent greater than one, is still a vector
  x <- rnorm(16)
  expect_equal(wavedec(matrix(x, ncol = 1)), wavedec(x))
  expect_equal(waverec(matrix(x, ncol = 1)), waverec(x))
  expect_equal(wavedec(array(x, c(1, 16, 1))), wavedec(x))
})

test_that("an unavailable filter size lists the available ones", {
  x <- sort(runif(16))
  expect_error(PHI(x, J = 3, family = "Coiflets", filter.size = 8),
               "6, 12, 18, ..., 96 and 102", fixed = TRUE)
  expect_error(PHI(x, J = 3, family = "Symmlets", filter.size = 2),
               "4, 6, 8, ..., 38 and 40", fixed = TRUE)
  expect_error(PHI(x, J = 3, family = "Daublets", filter.size = 3),
               "2, 4, 6, ..., 74 and 76", fixed = TRUE)
})
