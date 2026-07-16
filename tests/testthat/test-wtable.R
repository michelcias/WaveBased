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

# The 8-tap Daublet filter, used to exercise the user-supplied filter path.
d8 <- c( 0.230377813308896500863291183044070850,
         0.714846570552915647089921955273992604,
         0.630880767929858907881716338300615220,
        -0.027983769416859854211413747180075385,
        -0.187034811719093084079570672789081420,
         0.030841381835560763627219362534959050,
         0.032883011666885199735407513549244389,
        -0.010597401785069032104883208524027229)

test_that("wtable supports user-supplied filters (family = 'Own')", {
  tab.own <- wtable(family = "Own", wavelet.filter = d8, ngrid = 256,
                    check = FALSE)
  tab.d8 <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
                   check = FALSE)

  expect_equal(tab.own$family, "Own")
  expect_equal(tab.own$filter.size, 8)
  # Tabulated values must coincide with the built-in Daublet-8 tables
  expect_equal(tab.own$phi, tab.d8$phi, tolerance = 1e-12)
  expect_equal(tab.own$psi, tab.d8$psi, tolerance = 1e-12)

  # End to end through PHI/wbasis with the custom-filter table
  set.seed(3)
  x <- sort(runif(50))
  expect_equal(PHI(x, J = 3, family = "Own", wavelet.filter = d8,
                   wavelet.table = tab.own),
               PHI(x, J = 3, family = "daublets", filter.size = 8,
                   wavelet.table = tab.d8),
               tolerance = 1e-12)
  expect_equal(wbasis(x, j0 = 0, J = 3, family = "Own", wavelet.filter = d8,
                      wavelet.table = tab.own),
               wbasis(x, j0 = 0, J = 3, family = "daublets", filter.size = 8,
                      wavelet.table = tab.d8),
               tolerance = 1e-12)

  # A table built with a different custom filter is rejected
  expect_error(PHI(x, J = 3, family = "Own", wavelet.filter = rev(d8),
                   wavelet.table = tab.own), "different")

  # Odd-length filters are rejected
  expect_error(wtable(family = "Own", wavelet.filter = rep(0.1, 5)), "even")
})

test_that("wtable caches tables on disk and reuses them", {
  old <- options(WaveBased.cache.dir = file.path(tempdir(), "wb-cache-test"))
  on.exit({
    suppressMessages(wtable_cache(clear = TRUE))
    options(old)
  })

  # First call generates and saves; second call must return the same table
  t1 <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
               check = FALSE, cache = TRUE)
  listing <- wtable_cache()
  expect_equal(nrow(listing), 1)
  expect_match(listing$file, "^wtable_daublets_fs8_prec30_g256\\.rds$")

  t2 <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
               check = FALSE, cache = TRUE)
  expect_identical(t1$phi, t2$phi)
  expect_identical(t1$psi, t2$psi)

  # check = TRUE on a cached table measures the error and persists it; the
  # coarse test grid makes the self-check warn, which is asserted too
  expect_warning(
    t3 <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
                 check = TRUE, cache = TRUE),
    "exceeds")
  expect_false(is.na(t3$max.error))
  reread <- readRDS(file.path(getOption("WaveBased.cache.dir"),
                              listing$file[1]))
  expect_false(is.na(reread$max.error))

  # Different configurations get different files; Own filters are hashed
  wtable(family = "symmlets", filter.size = 8, ngrid = 256,
         check = FALSE, cache = TRUE)
  wtable(family = "Own", wavelet.filter = d8, ngrid = 256,
         check = FALSE, cache = TRUE)
  expect_equal(nrow(wtable_cache()), 3)

  # A corrupted cache file is regenerated instead of failing
  bad <- file.path(getOption("WaveBased.cache.dir"),
                   "wtable_daublets_fs8_prec30_g256.rds")
  saveRDS(list(rubbish = TRUE), bad)
  t4 <- wtable(family = "daublets", filter.size = 8, ngrid = 256,
               check = FALSE, cache = TRUE)
  expect_s3_class(t4, "wavelet_table")
  expect_identical(t4$phi, t1$phi)

  # Clearing empties the cache
  expect_message(wtable_cache(clear = TRUE), "Removed 3")
  expect_equal(nrow(suppressMessages(wtable_cache())), 0)
})
