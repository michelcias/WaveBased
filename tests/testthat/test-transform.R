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
