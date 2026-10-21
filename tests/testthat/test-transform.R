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

test_that("wavedec and waverec stop on a length that is not a power of two", {
  set.seed(2)
  for (n in c(60, 100, 1000, 1)) {
    x <- rnorm(n)
    expect_error(wavedec(x), "power of two", info = paste("n =", n))
    expect_error(waverec(x), "power of two", info = paste("n =", n))
    expect_error(wavedec(x, j0 = 4, family = "d", filter.size = 8,
                         boundary = "interval"),
                 "power of two", info = paste("interval, n =", n))
    expect_error(waverec(x, j0 = 4, family = "d", filter.size = 8,
                         boundary = "interval"),
                 "power of two", info = paste("interval, n =", n))
  }

  # the C entry points check the length themselves, for direct callers
  x <- rnorm(100)
  expect_error(.Call("_WaveBased_C_WaveDec", x, 1L, 20L, 0L, 0, 1L, NULL,
                     PACKAGE = "WaveBased"), "power of two")
  expect_error(.Call("_WaveBased_C_WaveRec", x, 1L, 20L, 0L, 0, 1L, NULL,
                     PACKAGE = "WaveBased"), "power of two")

  # a corrupted heap only shows up at the next garbage collection
  invisible(gc())
})

test_that("dyadic lengths still reconstruct exactly", {
  set.seed(3)
  for (n in c(64, 1024)) {
    x <- rnorm(n)
    expect_lt(max(abs(waverec(wavedec(x)) - x)), 1e-10)
    wdx <- wavedec(x, j0 = 4, family = "d", filter.size = 8,
                   boundary = "interval")
    wrx <- waverec(wdx, j0 = 4, family = "d", filter.size = 8,
                   boundary = "interval")
    expect_lt(max(abs(wrx - x)), 1e-10)
  }
})

test_that("a matrix is transformed by rows", {
  set.seed(5)
  X <- matrix(rnorm(6 * 64), 6, 64, dimnames = list(letters[1:6], NULL))

  for (fam in c("daublets", "symmlets", "coiflets")) {
    WX <- wavedec(X, j0 = 2, family = fam, filter.size = 6)
    expect_identical(WX, t(apply(X, 1, wavedec, j0 = 2, family = fam,
                                 filter.size = 6)))
    XR <- waverec(WX, j0 = 2, family = fam, filter.size = 6)
    expect_identical(XR, t(apply(WX, 1, waverec, j0 = 2, family = fam,
                                 filter.size = 6)))
    expect_lt(max(abs(XR - X)), 1e-10)
  }

  WX <- wavedec(X, j0 = 4, family = "d", filter.size = 8,
                boundary = "interval")
  expect_identical(WX, t(apply(X, 1, wavedec, j0 = 4, family = "d",
                               filter.size = 8, boundary = "interval")))
  XR <- waverec(WX, j0 = 4, family = "d", filter.size = 8,
                boundary = "interval")
  expect_lt(max(abs(XR - X)), 1e-10)

  # a single row stays a matrix
  W1 <- wavedec(X[1, , drop = FALSE], j0 = 2)
  expect_equal(dim(W1), c(1L, 64L))
  expect_identical(W1[1, ], wavedec(X[1, ], j0 = 2))
  expect_equal(dim(wavedec(matrix(0, 0, 16))), c(0L, 16L))

  # no decomposition when 2^j0 is the number of columns
  X8 <- matrix(rnorm(24), 3, 8)
  expect_warning(W8 <- wavedec(X8, j0 = 3), "no decomposition")
  expect_equal(W8, X8)

  # the number of columns must be a power of two
  expect_error(wavedec(matrix(0, 3, 60)), "power of two")
  expect_error(waverec(matrix(0, 3, 60)), "power of two")
  expect_error(.Call("_WaveBased_C_WaveDec", matrix(0, 3, 60), 1L, 20L, 0L,
                     0, 1L, NULL, PACKAGE = "WaveBased"), "power of two")
  invisible(gc())
})
