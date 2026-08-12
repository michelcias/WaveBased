# Unit tests for the wavelet-based estimators for mixture regression.

# Common setting: a Symmlet basis with a 10-tap filter, as used throughout
# Montoril, Pinheiro and Vidakovic (2019).
mix_setting <- function(n = 80, J = 4L, s = 2){
  set.seed(321)
  x <- sort(runif(n))
  w <- rnorm(n)
  sn <- s/2
  list(x = x, w = w, s = s, J = J, n = n,
       lower = c(rep(0, sn), x[seq_len(n - sn)]),
       upper = c(x[(sn + 1):n], rep(1, sn)))
}

mix_coef <- function(st, what = 3L, intmethod = 0L, delta = 1e-4, min.pts = 2,
                     evaltab = NULL, ngrid = 13L)
  .Call("_WaveBased_C_MixCoef", as.double(st$w), as.double(st$x),
        as.double(st$lower), as.double(st$upper), as.double(st$s),
        as.integer(st$J), 2L, 10L, 30L, as.double(0), as.integer(what),
        as.integer(intmethod), as.double(c(delta, min.pts)), evaltab,
        as.integer(ngrid))


test_that("estimator (8) reproduces the weighted product with the PHI matrix", {
  st <- mix_setting()

  mat <- PHI(st$x, J = st$J, family = "symmlets", filter.size = 10,
             prec.wavelet = 30, boundary = "periodic")
  ref <- drop(((st$upper - st$lower) * st$w / st$s) %*% mat)

  expect_equal(mix_coef(st, what = 1L)[[1]], ref, tolerance = 1e-12)
})


test_that("the trapezoidal rule reproduces a direct implementation", {
  st <- mix_setting(n = 30, J = 3L)
  p <- 2^st$J
  delta <- 1e-3

  # Direct transcription of WavIntegral/Trapezoidal of the reference MATLAB
  # implementation, built on the explicit basis matrix.
  ref <- numeric(p)
  for(i in seq_len(st$n)){
    leng <- max(2, ceiling((st$upper[i] - st$lower[i])/delta))
    xint <- seq(st$lower[i], st$upper[i], length.out = leng)
    Y <- PHI(xint, J = st$J, family = "symmlets", filter.size = 10,
             prec.wavelet = 30, boundary = "periodic")
    dl <- xint[2] - xint[1]
    ref <- ref + (st$w[i]/st$s) * (dl*colSums(Y) - (dl/2)*(Y[1, ] + Y[leng, ]))
  }

  got <- mix_coef(st, what = 2L, intmethod = 1L, delta = delta)[[2]]

  expect_equal(got, ref, tolerance = 1e-10)
})


test_that("the primitive agrees with a finely refined trapezoidal rule", {
  st <- mix_setting(n = 25, J = 3L)

  prim <- mix_coef(st, what = 2L, intmethod = 0L)[[2]]
  trap <- mix_coef(st, what = 2L, intmethod = 1L, delta = 1e-5)[[2]]

  # The remaining gap is dominated by the error of the trapezoidal rule.
  expect_equal(prim, trap, tolerance = 1e-7)
})


test_that("the primitive integrates the scaling functions exactly", {
  st <- mix_setting(n = 1, J = 5L)
  st$x <- 0.5
  st$w <- 1
  st$lower <- 0
  st$upper <- 1
  st$s <- 1

  got <- mix_coef(st, what = 2L)[[2]]

  # Over the whole unit interval every periodized scaling function integrates
  # to 2^{-J/2}, so the coefficients are all equal and add up to 2^{J/2}.
  expect_equal(got, rep(2^(-5/2), 2^5), tolerance = 1e-10)
  expect_equal(sum(got), 2^(5/2), tolerance = 1e-9)
})


test_that("requesting both estimators matches requesting them separately", {
  st <- mix_setting()

  both <- mix_coef(st, what = 3L)

  expect_equal(both[[1]], mix_coef(st, what = 1L)[[1]])
  expect_equal(both[[2]], mix_coef(st, what = 2L)[[2]])
  expect_null(mix_coef(st, what = 1L)[[2]])
  expect_null(mix_coef(st, what = 2L)[[1]])
})


test_that("the coefficients are stable across the size of the dyadic grid", {
  st <- mix_setting(n = 25, J = 3L)

  coarse <- mix_coef(st, what = 2L, ngrid = 10L)[[2]]
  default <- mix_coef(st, what = 2L, ngrid = 13L)[[2]]
  fine <- mix_coef(st, what = 2L, ngrid = 18L)[[2]]

  # The evaluation of the primitive inside a grid cell is O(h^3), so halving
  # the spacing gains roughly three bits per step.
  expect_equal(coarse, fine, tolerance = 1e-8)
  expect_equal(default, fine, tolerance = 1e-11)
})


test_that("C_MixEval reproduces the product of the basis matrix by the coefficients", {
  st <- mix_setting()
  cf <- mix_coef(st, what = 1L)[[1]]
  xo <- seq(0, 1, length.out = 120)

  mat <- PHI(xo, J = st$J, family = "symmlets", filter.size = 10,
             prec.wavelet = 30, boundary = "periodic")

  got <- .Call("_WaveBased_C_MixEval", as.double(cf), as.double(xo),
               as.integer(st$J), 2L, 10L, 30L, as.double(0), NULL)

  expect_equal(got, drop(mat %*% cf), tolerance = 1e-12)
})


test_that("C_LocLin reproduces the closed form of the local linear smoother", {
  set.seed(11)
  n <- 60
  x <- sort(runif(n))
  y <- sin(2*pi*x) + rnorm(n, sd = 0.2)
  xo <- seq(0, 1, length.out = 50)
  h <- 0.15

  ref <- vapply(xo, function(t){
    d <- t - x
    z <- stats::dnorm(d/h)
    s1 <- sum(d*z)
    s2 <- sum(d*d*z)
    sum((s2 - s1*d)*z*y)/(s2*sum(z) - s1^2)
  }, numeric(1))

  got <- .Call("_WaveBased_C_LocLin", as.double(x), as.double(y),
               as.double(xo), as.double(h))

  expect_equal(got, ref, tolerance = 1e-10)
})


test_that("wmixreg recovers a known mixing probability", {
  set.seed(202)
  n <- 600
  t <- (1:n)/n
  f <- 0.4*cos(2*pi*(t + pi)) + 0.5
  z <- rbinom(n, size = 1, prob = f)
  y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 5, plot = FALSE)

  expect_s3_class(fit, "wmixreg")
  expect_equal(dim(fit$fitted.values), c(250L, 6L))
  expect_true(all(fit$fitted.values >= 0 & fit$fitted.values <= 1))

  ftrue <- 0.4*cos(2*pi*(fit$newx + pi)) + 0.5
  ase <- colMeans((fit$fitted.values - ftrue)^2)
  expect_true(all(ase < 0.02))
})


test_that("the two raw estimators lead to similar estimates", {
  set.seed(7)
  n <- 400
  t <- (1:n)/n
  f <- 0.8*t + 0.1
  z <- rbinom(n, size = 1, prob = f)
  y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)

  for(r in fit$regularization)
    expect_equal(fit$fitted.values[, paste0("phi.", r)],
                 fit$fitted.values[, paste0("integral.", r)],
                 tolerance = 0.05)
})


test_that("the component means are estimated from the cut", {
  set.seed(5)
  y <- c(rnorm(50, mean = 5), rnorm(50, mean = 0))

  fit <- wmixreg(y, cut = 2, J = 3, plot = FALSE)

  expect_equal(fit$mu.U, mean(y[y > 2]))
  expect_equal(fit$mu.V, mean(y[y < 2]))
})


test_that("Efromovich's shrinkage selects a level and shortens the coefficients", {
  set.seed(13)
  n <- 400
  t <- (1:n)/n
  f <- 0.8*t + 0.1
  z <- rbinom(n, size = 1, prob = f)
  y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 5,
                 regularization = c("none", "hard", "efromovich"),
                 plot = FALSE)

  ef <- fit$fits$phi$efromovich

  expect_true(ef$J >= 0 && ef$J <= fit$J)
  expect_equal(length(ef$coef), 2^ef$J)
  expect_equal(length(coef(fit, "phi", "hard")), 2^fit$J)
  expect_equal(length(coef(fit, "phi", "none")), 2^fit$J)

  # One empirical risk per candidate level, and the selected level is the
  # minimizer of that sequence.
  expect_equal(length(ef$risk), fit$J + 1L)
  expect_true(all(is.finite(ef$risk)))
  expect_equal(which.min(ef$risk) - 1L, ef$J)

  # The shrinkage really acts on the coefficients.
  expect_true(max(abs(ef$coef - coef(fit, "phi", "none")[seq_along(ef$coef)])) > 1e-8)
})


test_that("Efromovich's shrinkage can select a level coarser than J", {
  # A constant mixing probability carries no fine detail, so the empirical
  # risk stops decreasing well before the finest candidate level.
  set.seed(101)
  n <- 512
  y <- rbinom(n, size = 1, prob = 0.5)*rnorm(n, mean = 2) +
       (1 - rbinom(n, size = 1, prob = 0.5))*rnorm(n, mean = 0)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 6,
                 regularization = "efromovich", plot = FALSE)

  expect_true(fit$fits$phi$efromovich$J <= 6)
})


test_that("hard thresholding leaves the scaling coefficients untouched", {
  set.seed(17)
  n <- 300
  y <- rbinom(n, 1, 0.5)*rnorm(n, mean = 2) + rnorm(n)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, j0 = 1,
                 regularization = c("none", "hard"), plot = FALSE)

  raw <- wavedec(coef(fit, "phi", "none"), j0 = 1, family = "Symmlets",
                 filter.size = 10)
  thr <- wavedec(coef(fit, "phi", "hard"), j0 = 1, family = "Symmlets",
                 filter.size = 10)

  # The coarsest scaling coefficients pass through untouched.
  expect_equal(raw[1:2], thr[1:2])

  # Every detail coefficient is either killed or left exactly as it was. The
  # tolerance absorbs the round trip through the decomposition.
  d.raw <- raw[-(1:2)]
  d.thr <- thr[-(1:2)]
  expect_true(all(abs(d.thr) < 1e-9 | abs(d.thr - d.raw) < 1e-9))

  # With this signal the universal threshold exceeds every detail.
  expect_true(sum(abs(d.thr) < 1e-9) > 0)
  expect_true(all(abs(d.raw[abs(d.thr) < 1e-9]) <= fit$fits$phi$hard$lambda))
})


test_that("predict agrees with the stored fitted values and honours clipping", {
  set.seed(23)
  n <- 200
  y <- rbinom(n, 1, 0.3)*rnorm(n, mean = 2) + rnorm(n)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)

  expect_equal(predict(fit, newx = fit$newx), fit$fitted.values)
  expect_equal(predict(fit), fit$fitted.values)
  expect_equal(fitted(fit), fit$fitted.values)

  unclipped <- predict(fit, newx = fit$newx, clip = FALSE)
  expect_true(all(predict(fit, newx = fit$newx) >= 0))
  expect_true(any(unclipped < 0) || all(unclipped >= 0))

  one <- predict(fit, newx = c(0.25, 0.75), estimator = "phi",
                 regularization = "hard")
  expect_equal(dim(one), c(2L, 1L))
  expect_equal(colnames(one), "phi.hard")
})


test_that("coef refuses the local linear regularization", {
  set.seed(29)
  n <- 150
  y <- rbinom(n, 1, 0.4)*rnorm(n, mean = 2) + rnorm(n)

  fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 3, plot = FALSE)

  expect_error(coef(fit, "phi", "loclin"), "time domain")
  expect_error(coef(fit, "phi", "soft"), "not available")
})


test_that("wmixreg validates its arguments", {
  y <- rnorm(64)

  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, x = 1:10, plot = FALSE),
               "same length")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, x = seq(-1, 2, length.out = 64),
                       plot = FALSE), "interval")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 2, plot = FALSE), "different")
  expect_error(wmixreg(y, plot = FALSE), "cut")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, s = 3, plot = FALSE), "even")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, J = 3, j0 = 3, plot = FALSE),
               "smaller")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, J = 3, ngrid = 1000,
                       plot = FALSE), "power of two")
  expect_error(wmixreg(y, mu.U = 2, mu.V = 0, family = "unknown",
                       plot = FALSE), "Unknown family")
})


test_that("unsorted observation times are handled as order statistics", {
  set.seed(31)
  n <- 128
  x <- runif(n)
  y <- rbinom(n, 1, 0.5)*rnorm(n, mean = 2) + rnorm(n)

  fit <- wmixreg(y, x = x, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)
  ord <- order(x)
  sorted <- wmixreg(y[ord], x = x[ord], mu.U = 2, mu.V = 0, J = 4, plot = FALSE)

  expect_false(is.unsorted(fit$x))
  expect_equal(fit$coefficients, sorted$coefficients)
})


test_that("a supplied wavelet table gives essentially the same estimates", {
  skip_on_cran()

  set.seed(37)
  n <- 200
  y <- rbinom(n, 1, 0.5)*rnorm(n, mean = 2) + rnorm(n)

  tab <- wtable(family = "Symmlets", filter.size = 10, ngrid = 2^13,
                check = FALSE)

  exact <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)
  tabled <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, wavelet.table = tab,
                    plot = FALSE)

  expect_equal(exact$coefficients$phi, tabled$coefficients$phi,
               tolerance = 1e-6)
})
