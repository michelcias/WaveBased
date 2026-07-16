# Unit tests for the boundary-corrected (Cohen-Daubechies-Vial) basis,
# boundary = "interval".

# The construction is validated through mathematical invariants: the derived
# blocks must assemble into an orthogonal one-step transform, the C kernels
# must match the pure R reference exactly, and the basis must reproduce
# polynomials up to degree filter.size/2 - 1 on the whole interval.

test_that("CDV blocks assemble into an orthogonal one-step transform", {
  for (case in list(c(1, 8), c(2, 12))) {
    h <- WaveBased:::.wb_filter(case[1], case[2])
    blk <- WaveBased:::.cdv_blocks(h)
    Nv <- case[2] / 2

    # two-scale orthonormality of the edge blocks
    expect_lt(max(abs(blk$BL %*% t(blk$BL) + blk$bL %*% t(blk$bL) - diag(Nv))),
              1e-9)
    expect_lt(max(abs(blk$UL %*% t(blk$UL) + blk$uL %*% t(blk$uL) - diag(Nv))),
              1e-9)

    # full one-step analysis matrix at the smallest admissible size
    jmin <- WaveBased:::.cdv_min_level(case[2], blk$uwidth, "decompose")
    n_in <- 2^(jmin + 1)
    O <- WaveBased:::.cdv_step_matrix(h, blk, n_in)
    expect_lt(max(abs(O %*% t(O) - diag(n_in))), 1e-9)
  }
})

test_that("wbasis with boundary = 'interval' matches the pure R reference", {
  set.seed(7)
  x <- c(0, 1, sort(runif(80)))
  h <- WaveBased:::.wb_filter(1, 8)
  blk <- WaveBased:::.cdv_blocks(h)
  jmin <- WaveBased:::.cdv_min_level(8, blk$uwidth, "decompose")
  Jf <- jmin + 1

  V <- wbasis(x, j0 = Jf, J = Jf, family = "d", filter.size = 8,
              boundary = "interval")
  expect_equal(ncol(V), 2^Jf)
  expect_true(all(is.finite(V)))

  # the C cascade must agree with the R step matrix applied to the finest
  # scaling basis (this validates WaveDec1CDV and the full column layout)
  W <- wbasis(x, j0 = jmin, J = Jf, family = "d", filter.size = 8,
              boundary = "interval")
  O <- WaveBased:::.cdv_step_matrix(h, blk, 2^Jf)
  expect_lt(max(abs(W - V %*% t(O))), 1e-12)

  # PHI is the j0 == J case
  expect_equal(PHI(x, J = Jf, family = "d", filter.size = 8,
                   boundary = "interval"), V)

  # PSI must match the detail block, up to the Daubechies-Lagarias precision
  P <- PSI(x, J = jmin, family = "d", filter.size = 8, prec.wavelet = 50,
           boundary = "interval")
  V50 <- wbasis(x, j0 = Jf, J = Jf, family = "d", filter.size = 8,
                prec.wavelet = 50, boundary = "interval")
  expect_lt(max(abs(P - (V50 %*% t(O))[, 2^jmin + seq_len(2^jmin)])), 1e-10)
})

test_that("the interval basis reproduces polynomials exactly", {
  set.seed(11)
  x <- c(0, 1, sort(runif(150)))
  for (case in list(list(fam = "d", fs = 6), list(fam = "s", fs = 10))) {
    fam <- match(case$fam, c("d", "s"))
    h <- WaveBased:::.wb_filter(fam, case$fs)
    blk <- WaveBased:::.cdv_blocks(h)
    jmin <- WaveBased:::.cdv_min_level(case$fs, blk$uwidth, "decompose")

    W <- wbasis(x, j0 = jmin, J = jmin + 1, family = case$fam,
                filter.size = case$fs, boundary = "interval")
    for (a in 0:(case$fs / 2 - 1)) {
      fit <- lm.fit(W, x^a)
      expect_lt(max(abs(fit$residuals)), 1e-6)
    }
  }
})

test_that("the interval basis is orthonormal (quadrature check)", {
  ng <- 2^14
  xg <- (seq_len(ng) - 0.5) / ng
  W <- wbasis(xg, j0 = 4, J = 5, family = "d", filter.size = 8,
              boundary = "interval")
  G <- crossprod(W) / ng
  expect_lt(max(abs(G - diag(ncol(W)))), 1e-3)
})

test_that("interior columns agree with the periodic basis away from edges", {
  set.seed(3)
  J <- 5
  L <- 8
  # points where no boundary or wrapped function is active
  x <- runif(40, min = (L - 1) / 2^J + 0.01, max = 1 - (L - 1) / 2^J - 0.01)

  Wi <- PHI(x, J = J, family = "d", filter.size = L, boundary = "interval")
  Wp <- PHI(x, J = J, family = "d", filter.size = L, boundary = "periodic")

  # interval interior column m corresponds to periodic column k = m
  for (m in seq_len(2^J - L)) {
    expect_lt(max(abs(Wi[, L / 2 + m - 1 + 1] - Wp[, m + 1])), 1e-12)
  }
  # edge columns must vanish at interior points
  edge <- c(seq_len(L / 2), 2^J - seq_len(L / 2) + 1)
  expect_lt(max(abs(Wi[, edge])), 1e-12)
})

test_that("boundary argument validation and backward compatibility", {
  set.seed(5)
  x <- sort(runif(30))

  # backward compatible defaults
  expect_equal(wbasis(x, j0 = 3, J = 4, family = "d", filter.size = 8),
               wbasis(x, j0 = 3, J = 4, family = "d", filter.size = 8,
                      boundary = "periodic"))
  expect_equal(wbasis(x, j0 = 3, J = 4, family = "d", filter.size = 8,
                      periodic = FALSE),
               wbasis(x, j0 = 3, J = 4, family = "d", filter.size = 8,
                      boundary = "none"))
  expect_warning(wbasis(x, j0 = 3, J = 4, family = "d", filter.size = 8,
                        periodic = FALSE, boundary = "periodic"),
                 "overrides")

  # interval-specific validation
  expect_error(wbasis(x, j0 = 2, J = 5, family = "d", filter.size = 8,
                      boundary = "interval"), "j0 >=")
  expect_error(PSI(x, J = 2, family = "d", filter.size = 8,
                   boundary = "interval"), "J >=")
  expect_error(wbasis(x - 2, j0 = 4, J = 5, family = "d", filter.size = 8,
                      boundary = "interval"), "\\[0, 1\\]")
  expect_error(wbasis(x, j0 = 6, J = 7, family = "c", filter.size = 12,
                      boundary = "interval"), "Coiflets")
  expect_error(wbasis(x, j0 = 7, J = 8, family = "d", filter.size = 20,
                      boundary = "interval"), "double precision")
})

test_that("interval basis handles boundary points and NA observations", {
  x <- c(0, 1e-14, 0.5, 1 - 1e-14, 1)
  W <- wbasis(x, j0 = 4, J = 5, family = "d", filter.size = 8,
              boundary = "interval")
  expect_true(all(is.finite(W)))
  # the first left-edge function must be active at x = 0, none of the
  # right-edge functions may be
  expect_gt(sum(abs(W[1, 1:4])), 0)
  expect_equal(sum(abs(W[1, 29:32])), 0)
  expect_gt(sum(abs(W[5, 29:32])), 0)

  expect_warning(Wna <- wbasis(c(0.3, NA, 0.7), j0 = 4, J = 4, family = "d",
                               filter.size = 8, boundary = "interval"),
                 "not finite")
  expect_true(all(is.na(Wna[2, ])))
  expect_true(all(is.finite(Wna[c(1, 3), ])))
})
