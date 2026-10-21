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

test_that("every entry point checks the filter", {
  set.seed(2)
  x <- sort(runif(32))
  y <- rnorm(32)
  bad <- c(1, 0, 0, 1) / sqrt(2)

  expect_error(PSI(x, J = 3, family = "Own", wavelet.filter = bad), "Lawton")
  expect_error(wbasis(x, j0 = 1, J = 3, family = "Own", wavelet.filter = bad),
               "Lawton")
  expect_error(wavedec(y, family = "Own", wavelet.filter = bad), "Lawton")
  expect_error(waverec(y, family = "Own", wavelet.filter = bad), "Lawton")
  expect_error(wtable(family = "Own", wavelet.filter = bad), "Lawton")
  expect_error(bayesthresh(y, family = "Own", wavelet.filter = bad), "Lawton")
  yb <- rep(0:1, 16)
  expect_error(wall(cbind(x), yb, J = 3, family = "Own", wavelet.filter = bad),
               "Lawton")
  expect_error(cv.wall(cbind(x), yb, J = 2:3, family = "Own",
                       wavelet.filter = bad), "Lawton")
  expect_error(wdensity(x, wf = function(u) 0.1 + 0.9 * u, J1 = 3,
                        family = "Own", wavelet.filter = bad, plot = FALSE),
               "Lawton")
  # a filter given alongside a built-in family is checked as well
  expect_error(suppressWarnings(PHI(x, J = 3, family = "d",
                                    wavelet.filter = bad)), "Lawton")
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
