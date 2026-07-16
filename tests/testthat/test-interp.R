# Tests for the interpolated evaluation of PHI, PSI, wbasis and wdensity
# via the wavelet.table argument.

test_that("interpolated PHI/PSI match the exact evaluation across families and levels", {
  set.seed(42)
  x <- sort(runif(200))

  cases <- list(list(family = "daublets", filter.size = 8),
                list(family = "daublets", filter.size = 20),
                list(family = "symmlets", filter.size = 20),
                list(family = "coiflets", filter.size = 12))

  for(cs in cases){
    tab <- wtable(family = cs$family, filter.size = cs$filter.size,
                  ngrid = 2^13, check = FALSE)

    # A single table (independent of J) must serve every level. Matrix
    # entries carry a sqrt(2^J) factor, which scales the tolerance.
    for(J in c(1, 3, 5)){
      tol <- sqrt(2^J)*1e-5

      exact <- PHI(x, J = J, family = cs$family, filter.size = cs$filter.size)
      apx <- PHI(x, J = J, family = cs$family, filter.size = cs$filter.size,
                 wavelet.table = tab)
      expect_lt(max(abs(exact - apx)), tol)

      exact <- PSI(x, J = J, family = cs$family, filter.size = cs$filter.size)
      apx <- PSI(x, J = J, family = cs$family, filter.size = cs$filter.size,
                 wavelet.table = tab)
      expect_lt(max(abs(exact - apx)), tol)
    }
  }
})

test_that("interpolated wbasis matches the exact evaluation", {
  set.seed(7)
  x <- sort(runif(150))
  tab <- wtable(family = "daublets", filter.size = 12, ngrid = 2^13,
                check = FALSE)

  for(j0 in c(0, 2)){
    for(J in c(3, 5)){
      exact <- wbasis(x, j0 = j0, J = J, family = "daublets",
                      filter.size = 12)
      apx <- wbasis(x, j0 = j0, J = J, family = "daublets",
                    filter.size = 12, wavelet.table = tab)
      expect_lt(max(abs(exact - apx)), sqrt(2^J)*1e-5)
    }
  }

  # j0 == J: scaling functions only
  exact <- wbasis(x, j0 = 3, J = 3, family = "daublets", filter.size = 12)
  apx <- wbasis(x, j0 = 3, J = 3, family = "daublets", filter.size = 12,
                wavelet.table = tab)
  expect_lt(max(abs(exact - apx)), sqrt(2^3)*1e-5)

  # Non-periodic case
  exact <- wbasis(x, j0 = 0, J = 4, family = "daublets", filter.size = 12,
                  boundary = "none")
  apx <- wbasis(x, j0 = 0, J = 4, family = "daublets", filter.size = 12,
                boundary = "none", wavelet.table = tab)
  expect_lt(max(abs(exact - apx)), sqrt(2^4)*1e-5)
})

test_that("a mismatched wavelet table is rejected", {
  set.seed(1)
  x <- sort(runif(20))
  tab <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
                check = FALSE)

  expect_error(PHI(x, J = 3, family = "symmlets", filter.size = 8,
                   wavelet.table = tab), "family")
  expect_error(PHI(x, J = 3, family = "daublets", filter.size = 20,
                   wavelet.table = tab), "filter.size")
  expect_error(PHI(x, J = 3, family = "daublets", filter.size = 8,
                   wavelet.table = list(a = 1)), "wtable")
  expect_error(wbasis(x, j0 = 0, J = 3, family = "symmlets", filter.size = 8,
                      wavelet.table = tab), "family")
})

test_that("wdensity accepts a wavelet table and matches the exact estimate", {
  set.seed(11)
  d <- rbeta(400, 2, 5)
  tab <- wtable(family = "symmlets", filter.size = 20, ngrid = 2^12,
                check = FALSE)

  exact <- wdensity(d, wf = function(x) rep_len(1, length(x)),
                    power.dens = 0.5, J1 = 4, family = "symmlets",
                    filter.size = 20, warped = FALSE, rescale = TRUE,
                    thresh = "linear", plot = FALSE)
  apx <- wdensity(d, wf = function(x) rep_len(1, length(x)),
                  power.dens = 0.5, J1 = 4, family = "symmlets",
                  filter.size = 20, warped = FALSE, rescale = TRUE,
                  thresh = "linear", plot = FALSE, wavelet.table = tab)

  expect_equal(exact$y, apx$y)
  expect_lt(max(abs(exact$dens - apx$dens)), 1e-3)

  # Thresholded variant runs end to end with the table
  expect_silent(wdensity(d, wf = function(x) rep_len(1, length(x)),
                         power.dens = 0.5, J1 = 4, family = "symmlets",
                         filter.size = 20, warped = FALSE, rescale = TRUE,
                         thresh = "hard", plot = FALSE, wavelet.table = tab))
})
