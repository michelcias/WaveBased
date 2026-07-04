# Unit tests for the wavelet interpolation table generator wtable().

test_that("wtable returns a well-formed wavelet_table object", {
  tab <- wtable(family = "daublets", filter.size = 8, ngrid = 512,
                check = FALSE)

  expect_s3_class(tab, "wavelet_table")
  expect_equal(dim(tab$phi), c(7, 513))
  expect_equal(dim(tab$psi), c(7, 513))
  expect_true(all(is.finite(tab$phi)))
  expect_true(all(is.finite(tab$psi)))
  expect_equal(tab$family, "Daublets")
  expect_equal(tab$filter.size, 8)
  expect_true(is.na(tab$max.error))
})

test_that("tabulated phi and psi satisfy the basic wavelet identities", {
  tab <- wtable(family = "symmlets", filter.size = 8, ngrid = 512,
                check = FALSE)
  G <- tab$ngrid

  # Partition of unity: for each grid point w the column holds phi(w + c)
  # over all integer offsets c covering the support, so sum_k phi(x - k) = 1
  # implies unit column sums.
  expect_lt(max(abs(colSums(tab$phi) - 1)), 1e-6)

  # The table entries cover each support point exactly once on a grid of
  # step 1/G, so Riemann sums approximate the integrals: int(phi) = 1,
  # int(psi) = 0 and unit L2 norms.
  expect_lt(abs(sum(tab$phi[, seq_len(G)])/G - 1), 1e-6)
  expect_lt(abs(sum(tab$psi[, seq_len(G)])/G), 1e-6)
  expect_lt(abs(sum(tab$phi[, seq_len(G)]^2)/G - 1), 1e-6)
  expect_lt(abs(sum(tab$psi[, seq_len(G)]^2)/G - 1), 1e-6)
})

test_that("wtable measures the interpolation error when check = TRUE", {
  tab <- wtable(family = "symmlets", filter.size = 20, ngrid = 2^12,
                check = TRUE, check.points = 200)

  expect_false(is.na(tab$max.error))
  expect_lt(tab$max.error, 1e-5)
})

test_that("wtable validates its inputs", {
  expect_error(wtable(family = "not-a-family"), "Unknown family")
  expect_error(wtable(family = "Own"), "Provide your own filter")
  expect_error(wtable(family = "daublets", filter.size = 8, ngrid = 1),
               "ngrid")
  expect_warning(wtable(family = "daublets", filter.size = 4, ngrid = 256,
                        check = FALSE),
                 "not recommended")
})

test_that("print method reports the table summary", {
  tab <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
                check = FALSE)
  expect_output(print(tab), "Wavelet interpolation table")
  expect_output(print(tab), "Daublets")
  expect_output(print(tab), "not measured")
})
