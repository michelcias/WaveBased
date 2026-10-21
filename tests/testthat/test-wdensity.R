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
