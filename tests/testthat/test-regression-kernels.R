# Numerical regression tests for the computational kernels.
#
# These tests lock the EXACT numerical output (bit-level, via serialized
# snapshots) of the periodic DWT kernels (WaveDec1/WaveRec1), the basis
# construction (wav_basis.c, exact and table-lookup paths) and the wall
# design-matrix pipeline. Their purpose is to detect any unintended numeric
# change introduced by performance work on the C code or the R helpers.
#
# First run: the snapshots are recorded under tests/testthat/_snaps/ (with a
# warning); run the suite a second time to confirm, and COMMIT the _snaps/
# files. After an INTENTIONAL numerical change (e.g., replacing the periodic
# cascade by direct per-level evaluation, whose results differ by the
# truncation error of the Daubechies-Lagarias iteration, ~1e-8 at the
# default prec.wavelet = 30 -- an error both methods share), inspect the
# reported differences and refresh the snapshots with
# testthat::snapshot_accept("regression-kernels").

test_that("periodic wavedec/waverec are numerically unchanged", {
  set.seed(2026)
  n <- 64
  x <- cumsum(rnorm(n))

  dec <- wavedec(x, j0 = 0, family = "Daublets", filter.size = 20)
  rec <- waverec(dec, j0 = 0, family = "Daublets", filter.size = 20)
  expect_snapshot_value(dec, style = "serialize")
  expect_snapshot_value(rec, style = "serialize")

  # Odd filter alignment / small signal exercises the index wrap-around.
  dec8 <- wavedec(x[1:8], j0 = 0, family = "Symmlets", filter.size = 16)
  expect_snapshot_value(dec8, style = "serialize")
})

test_that("periodic decomposed basis (exact evaluation) is numerically unchanged", {
  set.seed(11)
  u <- runif(30)

  B <- wbasis(u, j0 = 0, J = 4, family = "Daublets", filter.size = 8,
              prec.wavelet = 30)
  expect_snapshot_value(B, style = "serialize")

  # j0 > 0 takes the multi-level cascade branch in C_WavBasis.
  B2 <- wbasis(u, j0 = 1, J = 5, family = "Symmlets", filter.size = 8,
               prec.wavelet = 30)
  expect_snapshot_value(B2, style = "serialize")
})

test_that("periodic decomposed basis (table lookup) is numerically unchanged", {
  set.seed(12)
  u <- runif(30)

  tab <- wtable(family = "Daublets", filter.size = 8, prec.wavelet = 30,
                check = FALSE)
  Bt <- wbasis(u, j0 = 0, J = 4, family = "Daublets", filter.size = 8,
               wavelet.table = tab)
  expect_snapshot_value(Bt, style = "serialize")
})

test_that("wall design-matrix pipeline is numerically unchanged", {
  set.seed(42)
  n <- 40
  x <- matrix(runif(2*n), n, 2, dimnames = list(NULL, c("X1", "X2")))
  J <- rep_len(4L, 2L)
  j0 <- 0L

  eps <- WaveBased:::.wall_eps(NULL, J, j0, "periodic", TRUE)
  rs <- WaveBased:::.wall_rescale_pars(x, eps, TRUE, "periodic")
  spec <- list(J = J, j0 = j0, boundary = "periodic", family = "Daublets",
               filter.size = 8, prec.wavelet = 30, wavelet.filter = NULL,
               wavelet.table = NULL, rescale = TRUE,
               location = rs$location, scale = rs$scale, eps = eps,
               drop.phi = TRUE, sparse = FALSE, xnames = colnames(x))

  design <- WaveBased:::.wall_design(x, spec, clip = FALSE)
  expect_snapshot_value(rs, style = "serialize")
  expect_snapshot_value(design, style = "serialize")

  # Prediction path (clip = TRUE) with points beyond the training range.
  xnew <- x
  xnew[1L, ] <- c(-0.2, 1.3)
  design.new <- WaveBased:::.wall_design(xnew, spec, clip = TRUE)
  expect_snapshot_value(design.new, style = "serialize")
})

test_that("direct sparse design equals the dense construction bit for bit", {
  set.seed(77)
  n <- 60
  x <- matrix(runif(2*n), n, 2, dimnames = list(NULL, c("X1", "X2")))

  build_spec <- function(J, j0, sparse, wtab = NULL){
    J <- rep_len(as.integer(J), 2L)
    eps <- WaveBased:::.wall_eps(NULL, J, j0, "periodic", TRUE)
    rs <- WaveBased:::.wall_rescale_pars(x, eps, TRUE, "periodic")
    list(J = J, j0 = j0, boundary = "periodic", family = "Daublets",
         filter.size = 8, prec.wavelet = 30, wavelet.filter = NULL,
         wavelet.table = wtab, rescale = TRUE,
         location = rs$location, scale = rs$scale, eps = eps,
         drop.phi = j0 == 0L, sparse = sparse, xnames = colnames(x))
  }

  # The old sparse path: dense evaluation converted by Matrix(, sparse).
  old_sparse <- function(spec){
    blocks <- lapply(1:2, function(l){
      u <- (x[, l] - spec$location[l])/spec$scale[l]
      B <- WaveBased:::.wall_wbasis(u, spec, spec$J[l])
      if(spec$drop.phi)
        B <- B[, -1L, drop = FALSE]
      colnames(B) <- WaveBased:::.wall_colnames(spec$xnames[l], spec$j0,
                                                spec$J[l], spec$drop.phi)
      Matrix::Matrix(B, sparse = TRUE)
    })
    do.call(cbind, blocks)
  }

  for(j0 in c(0L, 1L)){
    spec <- build_spec(4L, j0, sparse = TRUE)
    new <- WaveBased:::.wall_design(x, spec, clip = FALSE)
    old <- old_sparse(spec)
    expect_s4_class(new, "dgCMatrix")
    # Same values, same structure (bit-level agreement of the slots).
    expect_identical(new@i, old@i)
    expect_identical(new@p, old@p)
    expect_identical(new@x, old@x)
    expect_identical(dimnames(new), dimnames(old))
    # And bitwise equal to the dense design.
    dense <- WaveBased:::.wall_design(x, build_spec(4L, j0, sparse = FALSE),
                                      clip = FALSE)
    expect_identical(as.matrix(new), dense)
  }

  # Table-lookup path.
  tab <- wtable(family = "Daublets", filter.size = 8, prec.wavelet = 30,
                check = FALSE)
  spec <- build_spec(5L, 0L, sparse = TRUE, wtab = tab)
  expect_identical(as.matrix(WaveBased:::.wall_design(x, spec, clip = FALSE)),
                   WaveBased:::.wall_design(x, build_spec(5L, 0L, FALSE, tab),
                                            clip = FALSE))

  # Prediction path (clip = TRUE) with points beyond the training range.
  xnew <- x
  xnew[1L, ] <- c(-0.2, 1.3)
  spec <- build_spec(4L, 0L, sparse = TRUE)
  expect_identical(as.matrix(WaveBased:::.wall_design(xnew, spec, clip = TRUE)),
                   WaveBased:::.wall_design(xnew, build_spec(4L, 0L, FALSE),
                                            clip = TRUE))
})

test_that("share.design reuses the max(J) design without changing results", {
  set.seed(99)
  n <- 200
  x <- matrix(runif(2*n), n, 2)
  h <- 3*sin(2*pi*x[, 1]) + 4*(x[, 2] - 0.5)
  y <- rbinom(n, 1, 1/(1 + exp(-h)))
  foldid <- sample(rep_len(1:5, n))

  # Same fixed eps in both: the column subsetting must be exact, so the
  # whole cross-validation and the selected fit must coincide.
  epsfix <- 1.9^(-4)
  fit1 <- cv.wall(x, y, J = 2:4, filter.size = 8, foldid = foldid,
                  eps = epsfix, use.table = "never", sparse = "never")
  fit2 <- cv.wall(x, y, J = 2:4, filter.size = 8, foldid = foldid,
                  eps = epsfix, use.table = "never", sparse = "never",
                  share.design = TRUE)

  expect_equal(fit1$cvtab, fit2$cvtab)
  expect_identical(fit1$J.min, fit2$J.min)
  expect_equal(fit1$lambda.min, fit2$lambda.min)
  expect_equal(fit1$cvm.min, fit2$cvm.min)
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(predict(fit1, x[1:7, ], type = "response"),
               predict(fit2, x[1:7, ], type = "response"))
})

test_that("share.design defaults eps to 1.9^(-max(J))", {
  set.seed(100)
  n <- 150
  x <- matrix(runif(2*n), n, 2)
  y <- rbinom(n, 1, 1/(1 + exp(-5*(x[, 1] - 0.5))))
  foldid <- sample(rep_len(1:5, n))

  fit <- cv.wall(x, y, J = 2:4, filter.size = 8, foldid = foldid,
                 use.table = "never", share.design = TRUE)
  expect_equal(fit$wall.fit$eps, rep(1.9^(-4), 2))

  # And the stored fit predicts with the fixed eps, like a direct wall()
  # fit at (J.min, eps fixed).
  ref <- wall(x, y, J = fit$J.min, filter.size = 8, eps = 1.9^(-4),
              use.table = "never",
              lambda = fit$wall.fit$lambda)
  expect_equal(predict(fit, x[1:5, ], s = fit$lambda.min),
               predict(ref, x[1:5, ], s = fit$lambda.min))
})

test_that("parallel.J gives the same results as the sequential grid", {
  skip_if_not_installed("foreach")
  set.seed(101)
  n <- 150
  x <- matrix(runif(2*n), n, 2)
  y <- rbinom(n, 1, 1/(1 + exp(-4*(x[, 1] - 0.5))))
  foldid <- sample(rep_len(1:5, n))

  # registerDoSEQ() exercises the foreach path (collection and selection
  # after the loop) without requiring a real parallel backend.
  foreach::registerDoSEQ()
  fit1 <- cv.wall(x, y, J = 2:4, filter.size = 8, foldid = foldid,
                  use.table = "never", sparse = "never")
  fit2 <- cv.wall(x, y, J = 2:4, filter.size = 8, foldid = foldid,
                  use.table = "never", sparse = "never", parallel.J = TRUE)

  expect_equal(fit1$cvtab, fit2$cvtab)
  expect_identical(fit1$J.min, fit2$J.min)
  expect_equal(fit1$lambda.min, fit2$lambda.min)
  expect_equal(coef(fit1), coef(fit2))
})

test_that("precomputed ranges do not change the rescaling parameters", {
  set.seed(7)
  x <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, paste0("X", 1:3)))
  eps <- WaveBased:::.wall_eps(NULL, rep_len(3L, 3L), 0L, "periodic", TRUE)

  rs1 <- WaveBased:::.wall_rescale_pars(x, eps, TRUE, "periodic")
  rs2 <- WaveBased:::.wall_rescale_pars(x, eps, TRUE, "periodic",
                                        ranges = apply(x, 2L, range))
  expect_identical(rs1, rs2)
})
