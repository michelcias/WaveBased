# Unit tests for the wavelet basis constructors PHI, PSI and wbasis.

test_that("PHI returns a periodic scaling-function basis with the right shape", {
  set.seed(123)
  x <- sort(runif(50))
  J <- 3

  mat <- PHI(x, J = J, family = "daublets", filter.size = 18,
             prec.wavelet = 30, boundary = "periodic")

  expect_true(is.matrix(mat))
  # With boundary = "periodic" the columns correspond to k = 0, ..., 2^J - 1.
  expect_equal(nrow(mat), length(x))
  expect_equal(ncol(mat), 2^J)
  expect_true(all(is.finite(mat)))
})

test_that("PSI returns a periodic wavelet basis with the right shape", {
  set.seed(123)
  x <- sort(runif(40))
  J <- 3

  mat <- PSI(x, J = J, family = "daublets", filter.size = 18,
             prec.wavelet = 30, boundary = "periodic")

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
                prec.wavelet = 30, boundary = "periodic")

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), length(x))
  # 2^j0 scaling columns + sum_{j=j0}^{J-1} 2^j wavelet columns = 2^J.
  expect_equal(ncol(mat), 2^J)
  expect_true(all(is.finite(mat)))
})

test_that("the Daubechies-Lagarias algorithm gives the known values of phi", {
  # Daubechies' scaling function of four taps at the half-integers of its
  # support [0, 3]: (2 + sqrt(3))/4, 0 and (2 - sqrt(3))/4, in the order of
  # the orientation of the filter
  v <- PHI(0.5, J = 0, family = "Daublets", filter.size = 4,
           boundary = "none", prec.wavelet = 50)
  expect_equal(sort(drop(v)), sort(c(0, (2 - sqrt(3))/4, (2 + sqrt(3))/4)),
               tolerance = 1e-8)

  # the integer translates add up to one at every prec, so the periodic
  # basis at level J adds up to 2^(J/2)
  set.seed(2)
  x <- runif(200)
  for (k in list(c("d", 2), c("d", 10), c("s", 20), c("c", 12))) {
    P <- PHI(x, J = 4, family = k[1], filter.size = as.integer(k[2]))
    expect_equal(rowSums(P), rep(4, length(x)), tolerance = 1e-12)
  }
})

test_that("the products with PHI computed in C agree with its matrix", {
  set.seed(5)
  x <- c(runif(300), 0, 1)
  w <- rexp(length(x))
  h8 <- WaveBased:::.wb_filter(1L, 8L)
  tab <- wtable(family = "Daublets", filter.size = 10)
  cfgs <- list(
    list(J = 5, family = "Daublets", filter.size = 10),
    list(J = 3, family = "Symmlets", filter.size = 20),  # translates share columns
    list(J = 4, family = "Daublets", filter.size = 2),   # Haar
    list(J = 6, family = "Coiflets", filter.size = 12),
    list(J = 5, family = "Own", filter.size = 8, wavelet.filter = h8),
    list(J = 5, family = "Daublets", filter.size = 10, wavelet.table = tab),
    list(J = 6, family = "Symmlets", filter.size = 8, boundary = "interval"),
    list(J = 7, family = "Daublets", filter.size = 6, boundary = "interval"))

  for (k in cfgs) {
    if (is.null(k$boundary))
      k$boundary <- "periodic"
    P <- do.call(PHI, c(list(x = x), k))
    a <- WaveBased:::.phi_args(x, k$J, k$family, k$filter.size,
                               k$wavelet.filter,
                               has.filter = !is.null(k$wavelet.filter),
                               k$boundary, k$wavelet.table)
    cf <- rnorm(2^k$J)
    expect_equal(WaveBased:::.wb_phi_means(x, w, k$J, k$filter.size, 30, a),
                 colMeans(P * w), tolerance = 1e-12)
    expect_equal(WaveBased:::.wb_phi_prod(x, cf, k$J, k$filter.size, 30, a),
                 drop(P %*% cf), tolerance = 1e-12)
  }

  # a point that is not finite gives NA, as its row of PHI does; a weight
  # that is not finite gives NaN means
  a <- WaveBased:::.phi_args(x, 5, "Daublets", 10, NULL, has.filter = FALSE,
                             "periodic")
  expect_warning(m <- WaveBased:::.wb_phi_means(c(x[1:5], NA), rep(1, 6), 5,
                                                10, 30, a), "not finite")
  expect_true(all(is.na(m)))
  expect_warning(v <- WaveBased:::.wb_phi_prod(c(0.2, NA), rnorm(32), 5, 10,
                                               30, a), "not finite")
  expect_true(is.finite(v[1]) && is.na(v[2]))
  expect_true(all(is.nan(WaveBased:::.wb_phi_means(x[1:5], c(1, 1, Inf, 1, 1),
                                                   5, 10, 30, a))))
})
