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
