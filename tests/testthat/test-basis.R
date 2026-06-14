# Unit tests for the wavelet basis constructors PHI, PSI and wbasis.

test_that("PHI returns a periodic scaling-function basis with the right shape", {
  set.seed(123)
  x <- sort(runif(50))
  J <- 3

  mat <- PHI(x, J = J, family = "daublets", filter.size = 18,
             prec.wavelet = 30, periodic = TRUE)

  expect_true(is.matrix(mat))
  # With periodic = TRUE the columns correspond to k = 0, ..., 2^J - 1.
  expect_equal(nrow(mat), length(x))
  expect_equal(ncol(mat), 2^J)
  expect_true(all(is.finite(mat)))
})

test_that("PSI returns a periodic wavelet basis with the right shape", {
  set.seed(123)
  x <- sort(runif(40))
  J <- 3

  mat <- PSI(x, J = J, family = "daublets", filter.size = 18,
             prec.wavelet = 30, periodic = TRUE)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), length(x))
  expect_equal(ncol(mat), 2^J)
  expect_true(all(is.finite(mat)))
})

test_that("wbasis returns a 2^J-column basis when periodic", {
  set.seed(123)
  x <- sort(runif(60))
  j0 <- 0
  J <- 3

  mat <- wbasis(x, j0 = j0, J = J, family = "Daublets", filter.size = 18,
                prec.wavelet = 30, periodic = TRUE)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), length(x))
  # 2^j0 scaling columns + sum_{j=j0}^{J-1} 2^j wavelet columns = 2^J.
  expect_equal(ncol(mat), 2^J)
  expect_true(all(is.finite(mat)))
})
