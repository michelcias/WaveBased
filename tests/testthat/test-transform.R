# Unit tests for the discrete wavelet transform (wavedec) and its inverse
# (waverec).

test_that("waverec inverts wavedec (perfect reconstruction)", {
  set.seed(123)
  n <- 64                 # length must be a power of 2
  x <- sort(runif(n))

  wdx <- wavedec(x, j0 = 3, family = "daublets", filter.size = 18)
  wrx <- waverec(wdx, j0 = 3, family = "daublets", filter.size = 18)

  expect_length(wdx, n)
  expect_length(wrx, n)
  expect_equal(wrx, x, tolerance = 1e-8)
})

test_that("reconstruction works across families and coarsest levels", {
  set.seed(1)
  n <- 128
  x <- cumsum(rnorm(n))

  for (fam in c("daublets", "symmlets")) {
    for (j0 in c(0, 2, 4)) {
      wdx <- wavedec(x, j0 = j0, family = fam, filter.size = 8)
      wrx <- waverec(wdx, j0 = j0, family = fam, filter.size = 8)
      expect_equal(wrx, x, tolerance = 1e-8,
                   info = paste("family =", fam, "j0 =", j0))
    }
  }
})
