# Tests for wdensity

test_that("the thresholding of wdensity uses the filter of the projection", {
  set.seed(1)
  d <- rbeta(200, 2, 5)
  h <- WaveBased:::.wb_filter(1L, 8L)
  w <- function(x) 0.1 + 0.9 * x
  est <- function(...) wdensity(d, wf = w, power.dens = 0.5, J1 = 5,
                                warped = FALSE, plot = FALSE, ...)

  for (th in c("hard", "soft")) {
    own <- est(family = "Own", wavelet.filter = h, thresh = th)
    # the same coefficients as a built-in family: the same estimate
    expect_equal(own, est(family = "Daublets", filter.size = 8, thresh = th))
    # a filter given alongside a built-in family is the one used throughout
    expect_warning(mixed <- est(family = "Coiflets", filter.size = 6,
                                wavelet.filter = h, thresh = th),
                   "not 'Own'")
    expect_equal(mixed, own)
  }
})

test_that("wdensity removes the non-finite elements of data", {
  set.seed(1)
  d <- rbeta(200, 2, 5)
  w <- function(x) 0.1 + 0.9 * x
  est <- function(x) wdensity(x, wf = w, power.dens = 0.5, J1 = 5,
                              warped = FALSE, plot = FALSE)

  expect_warning(dna <- est(c(d[1:100], NA, Inf, d[101:200], -Inf, NaN)),
                 "not finite")
  expect_equal(dna, est(d))
  expect_error(suppressWarnings(est(c(NA, Inf, NaN))), "no finite")
})

test_that("wdensity does not warn about tied data with the default warping", {
  set.seed(1)
  d <- round(rbeta(200, 2, 5), 2)
  expect_gt(anyDuplicated(d), 0)
  w <- function(x) 0.1 + 0.9 * x
  expect_no_warning(wdensity(d, wf = w, power.dens = 1, J1 = 6, warped = TRUE,
                             plot = FALSE))
})

test_that("wdensity takes the power of the positive part of the estimate", {
  set.seed(1)
  d <- rbeta(200, 2, 5)
  w <- function(x) 0.1 + 0.9 * x

  # a negative estimate of f^a has no real power 1/a, and gave NaN
  for (a in c(0.75, 2)) for (wp in c(FALSE, TRUE)) {
    e <- wdensity(d, wf = w, power.dens = a, J1 = if (wp) 7 else 5,
                  warped = wp, plot = FALSE)
    expect_false(anyNA(e$dens))
    expect_true(all(e$dens >= 0))
  }
})

test_that("the fine grid of the correction agrees with PHI", {
  set.seed(3)
  cfgs <- list(list("d", 1L, 10, 5, "periodic"),
               list("s", 2L, 20, 3, "periodic"),   # filter longer than 2^J
               list("d", 1L, 2, 4, "periodic"),    # Haar
               list("c", 3L, 12, 6, "periodic"),
               list("s", 2L, 8, 5, "interval"),
               list("d", 1L, 20, 6, "interval"))
  for (k in cfgs) {
    J <- k[[4]]
    cf <- rnorm(2^J)
    u <- (0:2^(J + 3)) / 2^(J + 3)
    g <- WaveBased:::.wb_fine_grid(cf, J = J, r = 3, fam = k[[2]],
                                   filter.size = k[[3]], prec.wavelet = 30,
                                   boundary = k[[5]])
    expect_equal(g, drop(PHI(u, J = J, family = k[[1]], filter.size = k[[3]],
                             boundary = k[[5]]) %*% cf),
                 tolerance = 1e-10)
  }

  # a filter given as "Own" gives the same grid as its built-in family
  h <- WaveBased:::.wb_filter(1L, 8L)
  cf <- rnorm(32)
  expect_identical(
    WaveBased:::.wb_fine_grid(cf, J = 5, r = 3, fam = 4L, filter.size = 8,
                              prec.wavelet = 30, wavelet.filter = h),
    WaveBased:::.wb_fine_grid(cf, J = 5, r = 3, fam = 1L, filter.size = 8,
                              prec.wavelet = 30))
})

test_that("the constant of the correction solves its equation", {
  xi <- WaveBased:::.wb_bonafide_xi
  f <- c(-0.2, 0.1, 0.5, 2, 1.3, -0.4, 0.8)
  w <- rep(0.25, 7)

  expect_equal(xi(f, w), 0.15)                       # mass above one
  expect_equal(sum(w * pmax(f - xi(f, w), 0)), 1)
  expect_equal(xi(c(0.1, 0.2, 0.3), rep(1, 3)), -2/15)  # mass below one
  expect_equal(xi(c(1, 5, 2), c(0, 0.5, 0.5)), 3)    # zero weights ignored
  expect_error(xi(c(1, NA), c(1, 1)), "not finite")
})

test_that("the corrected estimate of wdensity is a density on its domain", {
  set.seed(1)
  d <- rbeta(300, 2, 5)
  w <- function(x) 0.1 + 0.9 * x
  trap <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

  # without warping, the domain is [0, 1] in the scale of the basis: with
  # rescale = TRUE, [min - a, max + a] in the scale of the data
  J1 <- 5
  eps <- 1.9^(-J1)
  a <- eps * diff(range(d)) / (1 - 2 * eps)
  for (pw in c(0.75, 1, 2)) {
    e <- wdensity(d, wf = w, power.dens = pw, J1 = J1, warped = FALSE,
                  from = min(d) - a, to = max(d) + a, length.obs = 2^11 + 1,
                  correct = TRUE, plot = FALSE)
    expect_true(all(e$dens >= 0))
    expect_equal(trap(e$y, e$dens), 1, tolerance = 1e-3)
  }
  # (the interval basis refuses the points that rounding puts above 1, so it
  # is checked without rescaling, where the domain is [0, 1] itself; its
  # thresholding would need J0 >= 5)
  for (pw in c(0.75, 1)) {
    e <- wdensity(d, wf = w, power.dens = pw, J1 = J1, warped = FALSE,
                  rescale = FALSE, boundary = "interval", thresh = "linear",
                  from = 0, to = 1,
                  length.obs = 2^11 + 1, correct = TRUE, plot = FALSE)
    expect_true(all(e$dens >= 0))
    expect_equal(trap(e$y, e$dens), 1, tolerance = 1e-3)
  }

  # with warping, the domain is the range of the data
  for (pw in c(0.5, 1)) {
    e <- wdensity(d, wf = w, power.dens = pw, J1 = 6, warped = TRUE,
                  length.obs = 2^11, correct = TRUE, plot = FALSE)
    expect_true(all(e$dens >= 0))
    expect_equal(trap(e$y, e$dens), 1, tolerance = 1e-3)
  }
  e <- wdensity(d, wf = w, power.dens = 1, J1 = 6, warped = TRUE,
                warp.fun = function(x) pbeta(x, 2, 4),
                dwarp.fun = function(x) dbeta(x, 2, 4),
                length.obs = 2^11, correct = TRUE, plot = FALSE)
  expect_equal(trap(e$y, e$dens), 1, tolerance = 1e-3)
})

test_that("the correction of wdensity does not depend on the evaluation grid", {
  set.seed(1)
  d <- rbeta(300, 2, 5)
  w <- function(x) 0.1 + 0.9 * x

  for (wp in c(FALSE, TRUE)) {
    est <- function(...) wdensity(d, wf = w, power.dens = 1,
                                  J1 = if (wp) 6 else 5, warped = wp,
                                  correct = TRUE, plot = FALSE, ...)
    expect_identical(est()$xi, est(from = quantile(d, 0.2),
                                   to = quantile(d, 0.8),
                                   length.obs = 50)$xi)
  }
})

test_that("with power.dens = 0.5 and no warping the estimate needs no correction", {
  set.seed(1)
  d <- rbeta(300, 2, 5)
  w <- function(x) 0.1 + 0.9 * x
  est <- function(...) wdensity(d, wf = w, power.dens = 0.5, J1 = 5,
                                warped = FALSE, plot = FALSE, ...)

  e <- est(correct = TRUE)
  expect_identical(e$xi, 0)
  expect_equal(e$dens, est()$dens)
  expect_null(est()$xi)
})
