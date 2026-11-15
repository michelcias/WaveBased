# Tests for the checks of a user-supplied filter (family = "Own")

# Filter of length 2M from M angles of the lattice of Vaidyanathan and Hoang
# (1988): it sums to sqrt(2) and is orthogonal to its even shifts whenever
# sum(theta) = pi/4, but it may fail Lawton's condition.
lattice_filter <- function(theta) {
  R <- function(t) matrix(c(cos(t), -sin(t), sin(t), cos(t)), 2)
  E <- lapply(1:2, function(i) lapply(1:2, function(j) R(theta[1])[i, j]))
  for (t in theta[-1]) {
    E[[1]] <- lapply(E[[1]], function(p) c(p, 0))
    E[[2]] <- lapply(E[[2]], function(p) c(0, p))
    Rk <- R(t)
    E <- lapply(1:2, function(i) lapply(1:2, function(j)
      Rk[i, 1] * E[[1]][[j]] + Rk[i, 2] * E[[2]][[j]]))
  }
  h <- numeric(2 * length(E[[1]][[1]]))
  h[c(TRUE, FALSE)] <- E[[1]][[1]]
  h[c(FALSE, TRUE)] <- E[[1]][[2]]
  h
}

db4 <- WaveBased:::.wb_filter(1L, 8L)

test_that("the filters of the built-in families pass the checks silently", {
  check <- WaveBased:::.wb_check_filter
  for (L in seq(2, 76, 2)) expect_silent(check(WaveBased:::.wb_filter(1L, L)))
  for (L in seq(4, 40, 2)) expect_silent(check(WaveBased:::.wb_filter(2L, L)))
  for (L in seq(6, 102, 6)) expect_silent(check(WaveBased:::.wb_filter(3L, L)))
  expect_silent(check(rev(db4)))
  expect_silent(check(lattice_filter(c(1, -0.5, pi/4 - 0.5))))
})

test_that("filters that are not orthonormal wavelet filters are refused", {
  set.seed(1)
  x <- sort(runif(20))
  phi <- function(h) PHI(x, J = 3, family = "Own", wavelet.filter = h)

  expect_error(phi(rep(0.1, 5)), "even")
  expect_error(phi(c(db4, NA)), "even")
  expect_error(phi(as.character(db4)), "even")
  # the high-pass filter, the filter summing to 1 and the one summing to
  # -sqrt(2) come with a hint
  expect_error(phi((-1)^(0:7) * rev(db4)), "high-pass")
  expect_error(phi(db4 / sqrt(2)), "multiply it by sqrt(2)", fixed = TRUE)
  expect_error(phi(-db4), "multiply it by -1")
  expect_error(phi(db4 + c(1e-4, rep(0, 7))), "orthonormal wavelet basis")
  # sum and orthogonality hold, but the translations of phi = 1[0,3)/3 are
  # not orthonormal
  expect_error(phi(c(1, 0, 0, 1) / sqrt(2)), "Lawton")
  expect_error(phi(lattice_filter(c(0, pi/4))), "Lawton")
})

test_that("an approximate filter and a nearly degenerate one are accepted with a warning", {
  set.seed(1)
  x <- sort(runif(20))
  phi <- function(h) PHI(x, J = 3, family = "Own", wavelet.filter = h)

  noisy <- db4 + rnorm(8) * 1e-7
  expect_warning(B <- phi(noisy), "up to")
  expect_equal(B, PHI(x, J = 3, family = "Daublets", filter.size = 8),
               tolerance = 1e-4)
  # the second eigenvalue of Lawton's matrix is at about theta^2 from 1
  expect_warning(phi(lattice_filter(c(0.01, pi/4 - 0.01))), "nearly fails")
})

test_that("every entry point that evaluates phi checks Lawton's condition", {
  set.seed(2)
  x <- sort(runif(32))
  y <- rnorm(32)
  bad <- c(1, 0, 0, 1) / sqrt(2)

  expect_error(PSI(x, J = 3, family = "Own", wavelet.filter = bad), "Lawton")
  expect_error(wbasis(x, j0 = 1, J = 3, family = "Own", wavelet.filter = bad),
               "Lawton")
  expect_error(wtable(family = "Own", wavelet.filter = bad), "Lawton")
  yb <- rep(0:1, 16)
  expect_error(wall(cbind(x), yb, J = 3, family = "Own", wavelet.filter = bad),
               "Lawton")
  expect_error(cv.wall(cbind(x), yb, J = 2:3, family = "Own",
                       wavelet.filter = bad), "Lawton")
  expect_error(wdensity(x, wf = function(u) 0.1 + 0.9 * u, J1 = 3,
                        family = "Own", wavelet.filter = bad, plot = FALSE),
               "Lawton")
  expect_error(wmixreg(y, mu.U = 1, mu.V = 0, J = 3, family = "Own",
                       wavelet.filter = bad, plot = FALSE), "Lawton")
  # a filter given alongside a built-in family is checked as well
  expect_error(suppressWarnings(PHI(x, J = 3, family = "d",
                                    wavelet.filter = bad)), "Lawton")
})

test_that("the functions that only transform do not check Lawton's condition", {
  set.seed(2)
  y <- rnorm(64)
  bad <- c(1, 0, 0, 1) / sqrt(2)

  # the other conditions make the transform orthogonal
  d <- expect_silent(wavedec(y, j0 = 2, family = "Own", wavelet.filter = bad))
  expect_equal(sum(d^2), sum(y^2))
  expect_equal(waverec(d, j0 = 2, family = "Own", wavelet.filter = bad), y)
  expect_silent(bayesthresh(y, family = "Own", wavelet.filter = bad))
  expect_silent(bwregime(y, family = "Own", wavelet.filter = bad,
                         nchain = 20, burn = 10, lag = 1))
  expect_silent(bwmixreg(y, family = "Own", wavelet.filter = bad,
                         nchain = 20, burn = 10, lag = 1))
  # the other checks remain
  expect_error(wavedec(y, family = "Own", wavelet.filter = bad / sqrt(2)),
               "multiply it by sqrt(2)", fixed = TRUE)
  expect_error(bayesthresh(y, family = "Own", wavelet.filter = -bad),
               "multiply it by -1")
})

test_that("check.lawton = FALSE accepts a degenerate filter", {
  set.seed(2)
  x <- sort(runif(32))
  y <- rnorm(32)
  bad <- c(1, 0, 0, 1) / sqrt(2)

  # phi = 1[0,3)/3, evaluated exactly by the Daubechies-Lagarias algorithm
  u <- c(0.1, 0.7, 1.4, 2.2, 2.9)
  B <- PHI(u, J = 0, family = "Own", wavelet.filter = bad, boundary = "none",
           check.lawton = FALSE)
  expect_equal(B, outer(u, 0:4, function(t, k) ((t - k + 2) >= 0 &
                                                  (t - k + 2) < 3) / 3),
               ignore_attr = TRUE)

  expect_silent(PSI(x, J = 3, family = "Own", wavelet.filter = bad,
                    check.lawton = FALSE))
  expect_silent(wbasis(x, j0 = 1, J = 3, family = "Own", wavelet.filter = bad,
                       check.lawton = FALSE))
  expect_s3_class(suppressWarnings(wtable(family = "Own", wavelet.filter = bad,
                                          check.lawton = FALSE)),
                  "wavelet_table")
  yb <- rep(0:1, 16)
  # the table warns of the short filter, not of Lawton's condition
  expect_warning(fit <- wall(cbind(x), yb, J = 3, family = "Own",
                             wavelet.filter = bad, use.table = "always",
                             check.lawton = FALSE), "very short filters")
  expect_silent(predict(fit, cbind(x), type = "response"))
  expect_silent(wall(cbind(x), yb, J = 3, family = "Own", wavelet.filter = bad,
                     use.table = "never", sparse = "never",
                     check.lawton = FALSE))
  expect_silent(cv.wall(cbind(x), yb, J = 2:3, family = "Own",
                        wavelet.filter = bad, nfolds = 3,
                        check.lawton = FALSE))
  expect_silent(wdensity(x, wf = function(u) 0.1 + 0.9 * u, J1 = 3,
                         family = "Own", wavelet.filter = bad, warped = FALSE,
                         plot = FALSE, check.lawton = FALSE))
  expect_silent(wmixreg(y, mu.U = 1, mu.V = 0, J = 3, family = "Own",
                        wavelet.filter = bad, plot = FALSE,
                        check.lawton = FALSE))

  # nor the warning of a nearly degenerate filter
  expect_silent(PHI(x, J = 3, family = "Own",
                    wavelet.filter = lattice_filter(c(0.01, pi/4 - 0.01)),
                    check.lawton = FALSE))
  # the other checks remain, and the argument is validated
  expect_error(PHI(x, J = 3, family = "Own", wavelet.filter = -bad,
                   check.lawton = FALSE), "multiply it by -1")
  expect_error(PHI(x, J = 3, family = "Own", wavelet.filter = bad,
                   check.lawton = NA), "TRUE or FALSE")
})

test_that("functions that call others warn once about an approximate filter", {
  set.seed(3)
  x <- matrix(runif(400), 200, 2)
  y <- rbinom(200, 1, plogis(3 * sin(2 * pi * x[, 1]) + 8 * (x[, 2] - 0.5)))
  noisy <- db4 + rnorm(8) * 1e-7
  count <- function(expr) {
    n <- 0L
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("up to", conditionMessage(w))) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    n
  }

  fit <- NULL
  expect_equal(count(fit <- wall(x, y, J = 3, family = "Own",
                                 wavelet.filter = noisy, use.table = "always")),
               1L)
  expect_equal(count(predict(fit, x, type = "response")), 0L)
  expect_equal(count(wall(x, y, J = 3, family = "Own", wavelet.filter = noisy,
                          use.table = "never", sparse = "never")), 1L)
  expect_equal(count(cv.wall(x, y, J = 2:3, family = "Own",
                             wavelet.filter = noisy, nfolds = 3)), 1L)
  t <- (1:600) / 600
  z <- rbinom(600, size = 1, prob = 0.4 * cos(2 * pi * (t + pi)) + 0.5)
  ym <- z * rnorm(600, mean = 2) + (1 - z) * rnorm(600)
  expect_equal(count(wmixreg(ym, mu.U = 2, mu.V = 0, J = 5, family = "Own",
                             wavelet.filter = noisy, plot = FALSE)), 1L)
  expect_equal(count(wdensity(x[, 1], wf = function(u) 0.1 + 0.9 * u, J1 = 4,
                              family = "Own", wavelet.filter = noisy,
                              warped = FALSE, plot = FALSE)), 1L)
})
