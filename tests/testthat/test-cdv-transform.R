# Unit tests for the boundary-corrected (Cohen-Daubechies-Vial) wavelet
# transform, wavedec/waverec with boundary = "interval", and for the
# precomputed high-precision boundary block tables.

test_that("the interval transform is orthogonal and inverts exactly", {
  set.seed(21)
  for (case in list(c(1, 8), c(2, 12))) {
    blk <- WaveBased:::.cdv_lookup(case[1], case[2])
    jmin <- WaveBased:::.cdv_min_level(case[2], blk$uwidth, "decompose")
    fam <- c("d", "s")[case[1]]
    n <- 2^(jmin + 2)
    x <- rnorm(n)

    w <- wavedec(x, j0 = jmin, family = fam, filter.size = case[2],
                 boundary = "interval")
    # Parseval: the transform is orthogonal
    expect_equal(sum(w^2), sum(x^2), tolerance = 1e-10)
    # exact round trip
    xr <- waverec(w, j0 = jmin, family = fam, filter.size = case[2],
                  boundary = "interval")
    expect_lt(max(abs(xr - x)), 1e-10)
  }
})

test_that("the interval transform matches the wbasis filter bank", {
  set.seed(22)
  x <- c(0, 1, sort(runif(50)))
  jmin <- WaveBased:::.cdv_min_level(8, WaveBased:::.cdv_lookup(1, 8)$uwidth,
                                     "decompose")
  Jf <- jmin + 1

  # decomposing the rows of the finest-level scaling basis must reproduce
  # the decomposed basis returned by wbasis
  V <- wbasis(x, j0 = Jf, J = Jf, family = "d", filter.size = 8,
              boundary = "interval")
  W <- wbasis(x, j0 = jmin, J = Jf, family = "d", filter.size = 8,
              boundary = "interval")
  Wt <- t(apply(V, 1, wavedec, j0 = jmin, family = "d", filter.size = 8,
                boundary = "interval"))
  expect_lt(max(abs(Wt - W)), 1e-12)
})

test_that("wavedec validates the interval boundary requirements", {
  x <- rnorm(32)
  expect_error(wavedec(x, j0 = 2, family = "d", filter.size = 8,
                       boundary = "interval"), "j0 >=")
  expect_error(wavedec(x, j0 = 3, family = "d", filter.size = 8,
                       boundary = "none"), "periodic")
  # default remains the periodic transform
  expect_equal(wavedec(x, j0 = 3, family = "d", filter.size = 8),
               wavedec(x, j0 = 3, family = "d", filter.size = 8,
                       boundary = "periodic"))
})

test_that("precomputed tables cover sizes beyond the runtime construction", {
  set.seed(23)
  x <- c(0, 1, sort(runif(120)))
  # Daublets 20 and 24 are impossible to derive at run time in double
  # precision, but are tabulated in high precision
  for (case in list(c(1, 20), c(1, 24), c(2, 24))) {
    blk <- WaveBased:::.cdv_lookup(case[1], case[2])
    expect_false(is.null(blk))
    h <- WaveBased:::.wb_filter(case[1], case[2])
    Nv <- case[2] / 2

    expect_lt(max(abs(blk$BL %*% t(blk$BL) + blk$bL %*% t(blk$bL) -
                        diag(Nv))), 1e-11)
    jmin <- WaveBased:::.cdv_min_level(case[2], blk$uwidth, "decompose")
    O <- WaveBased:::.cdv_step_matrix(h, blk, 2^(jmin + 1))
    expect_lt(max(abs(O %*% t(O) - diag(2^(jmin + 1)))), 1e-11)

    # end-to-end polynomial reproduction with the tabulated blocks
    fam <- c("d", "s")[case[1]]
    W <- wbasis(x, j0 = jmin, J = jmin + 1, family = fam,
                filter.size = case[2], boundary = "interval")
    for (a in c(0, 3, Nv - 1)) {
      fit <- lm.fit(W, x^a)
      expect_lt(max(abs(fit$residuals)), 1e-6)
    }
  }
})

test_that("tabulated blocks agree with the runtime construction", {
  # for a filter that both paths can handle, the two constructions must
  # produce the same basis up to the double-precision derivation error
  h <- WaveBased:::.wb_filter(1, 8)
  tab <- WaveBased:::.cdv_lookup(1, 8)
  run <- WaveBased:::.cdv_blocks(h)
  expect_equal(tab$uwidth, run$uwidth)
  for (nm in c("BL", "bL", "UL", "uL", "phi0L", "BR", "bR", "UR", "uR",
               "phi0R"))
    expect_lt(max(abs(tab[[nm]] - run[[nm]])), 1e-9)
})

test_that("wdensity works with the interval boundary", {
  data(bac, envir = environment())
  # suppressWarnings: approx() warns about tied values in 'bac', which is
  # unrelated to the boundary treatment being tested here
  m <- suppressWarnings(
    wdensity(data = bac, wf = function(x) 0.1 + 0.9 * x, power.dens = 0.5,
             J1 = 5, family = "d", filter.size = 8, warped = TRUE,
             thresh = "linear", plot = FALSE, boundary = "interval"))
  expect_true(all(is.finite(m$dens)))
  expect_true(all(m$dens >= 0))
})
