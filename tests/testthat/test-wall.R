# Unit tests for the wavelet-based additive logistic LASSO classifier
# (wall and cv.wall).

# Small simulated additive logistic model, shared by the tests below.
.wall_sim <- function(n = 200, seed = 123){
  set.seed(seed)
  x <- matrix(runif(2*n), n, 2)
  h <- 3*sin(2*pi*x[, 1]) + 8*(x[, 2] - 0.5)
  y <- rbinom(n, 1, 1/(1 + exp(-h)))
  list(x = x, y = y)
}

test_that("wall fits and returns a well-formed object", {
  dat <- .wall_sim()

  fit <- wall(dat$x, dat$y, J = 3, family = "Daublets", filter.size = 8)

  expect_s3_class(fit, "wall")
  expect_s3_class(fit$glmnet.fit, "glmnet")
  # With j0 = 0 and the periodized basis, the constant scaling column of
  # each covariate is dropped: p = d*(2^J - 1).
  expect_equal(fit$nvars, 2*(2^3 - 1))
  expect_equal(fit$J, c(3L, 3L))
  expect_true(fit$drop.phi)
  expect_equal(fit$eps, rep(1.9^(-3), 2))
  expect_equal(length(fit$penalty.factor), fit$nvars)
  expect_true(all(fit$penalty.factor == 1))  # all psi with j0 = 0
  expect_output(print(fit), "wall")
})

test_that("predict.wall returns consistent link, response and class", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)

  eta <- predict(fit, dat$x, s = 0.01)
  prob <- predict(fit, dat$x, s = 0.01, type = "response")
  cls <- predict(fit, dat$x, s = 0.01, type = "class")

  expect_equal(dim(eta), c(nrow(dat$x), 1L))
  expect_equal(drop(prob), drop(1/(1 + exp(-eta))), tolerance = 1e-10)
  expect_equal(drop(cls), as.character(as.integer(drop(eta) > 0)))
  # The classifier must beat random guessing by far on the training data.
  expect_lt(mean(drop(cls) != dat$y), 0.3)
})

test_that("factor responses keep their labels in class predictions", {
  dat <- .wall_sim()
  yf <- factor(ifelse(dat$y == 1, "case", "control"),
               levels = c("control", "case"))

  fit <- wall(dat$x, yf, J = 3, filter.size = 8)
  cls <- predict(fit, dat$x, s = 0.01, type = "class")

  expect_equal(fit$classnames, c("control", "case"))
  expect_true(all(cls %in% c("control", "case")))
  # Same fit as with the 0/1 coding of the second level.
  fit01 <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  expect_equal(coef(fit, s = 0.01), coef(fit01, s = 0.01))
})

test_that("coef.wall names the basis functions", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 2, filter.size = 8)

  beta <- coef(fit, s = 0.01)
  expect_equal(nrow(beta), fit$nvars + 1L)  # + intercept
  expect_true(all(grepl("^X[12]\\.psi[01]\\.", rownames(beta)[-1L])))
})

test_that("per-variable resolution levels are accepted", {
  dat <- .wall_sim()

  fit <- wall(dat$x, dat$y, J = c(3, 2), filter.size = 8)

  expect_equal(fit$J, c(3L, 2L))
  expect_equal(fit$nvars, (2^3 - 1) + (2^2 - 1))
  expect_equal(fit$eps, 1.9^(-c(3, 2)))
  pred <- predict(fit, dat$x[1:5, ], s = 0.01)
  expect_true(all(is.finite(pred)))
})

test_that("sparse and dense designs give the same fit", {
  dat <- .wall_sim()

  fit.d <- wall(dat$x, dat$y, J = 4, filter.size = 8, sparse = "never")
  fit.s <- wall(dat$x, dat$y, J = 4, filter.size = 8, sparse = "always")

  expect_false(fit.d$sparse)
  expect_true(fit.s$sparse)
  expect_equal(predict(fit.d, dat$x, s = 0.01),
               predict(fit.s, dat$x, s = 0.01), tolerance = 1e-6)
})

test_that("table lookup and exact evaluations agree", {
  dat <- .wall_sim()
  tab <- wtable(family = "Daublets", filter.size = 8, check = FALSE)

  fit.e <- wall(dat$x, dat$y, J = 3, filter.size = 8, use.table = "never")
  fit.t <- wall(dat$x, dat$y, J = 3, filter.size = 8, wavelet.table = tab)

  expect_null(fit.e$wavelet.table)
  expect_s3_class(fit.t$wavelet.table, "wavelet_table")
  expect_equal(predict(fit.e, dat$x, s = 0.01),
               predict(fit.t, dat$x, s = 0.01), tolerance = 1e-3)
})

test_that("new observations beyond the training range are truncated", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)

  outside <- rbind(c(-1, 2), c(2, -1))
  pred <- predict(fit, outside, s = 0.01)
  expect_true(all(is.finite(pred)))
  # Truncation maps them to the same predictions of the extreme points.
  extremes <- rbind(c(min(dat$x[, 1]), max(dat$x[, 2])),
                    c(max(dat$x[, 1]), min(dat$x[, 2])))
  expect_equal(pred, predict(fit, extremes, s = 0.01), tolerance = 1e-10)
})

test_that("the interval boundary basis is supported", {
  dat <- .wall_sim()

  fit <- wall(dat$x, dat$y, J = 4, j0 = 3, filter.size = 4,
              boundary = "interval")

  expect_false(fit$drop.phi)
  expect_equal(fit$eps, c(0, 0))            # no boundary strip needed
  expect_equal(fit$nvars, 2*2^4)            # phi blocks retained
  # Scaling functions unpenalized, wavelets penalized.
  expect_equal(sum(fit$penalty.factor == 0), 2*2^3)
  pred <- predict(fit, dat$x[1:5, ], s = 0.01)
  expect_true(all(is.finite(pred)))
})

test_that("wall validates its inputs", {
  dat <- .wall_sim()

  expect_error(wall(dat$x, dat$y), "'J' must be provided")
  expect_error(wall(dat$x, dat$y, J = 0), "larger than 'j0'")
  expect_error(wall(dat$x, dat$y, J = c(2, 3, 4)), "length 1 or one entry")
  expect_error(wall(dat$x, rep(1, nrow(dat$x)), J = 3), "two classes")
  expect_error(wall(dat$x, c(dat$y[-1], 2), J = 3), "two classes")
  expect_error(wall(cbind(dat$x, 1), dat$y, J = 3), "constant")
  xna <- dat$x; xna[1, 1] <- NA
  expect_error(wall(xna, dat$y, J = 3), "finite")
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  expect_error(predict(fit, dat$x[, 1]), "column")
})

test_that("cv.wall selects (J, lambda) and predicts", {
  dat <- .wall_sim(n = 300)

  set.seed(1)
  cvfit <- cv.wall(dat$x, dat$y, J = 1:4, filter.size = 8, nfolds = 5)

  expect_s3_class(cvfit, "cv.wall")
  expect_true(cvfit$J.min %in% 1:4)
  expect_s3_class(cvfit$wall.fit, "wall")
  expect_equal(cvfit$wall.fit$J, rep(cvfit$J.min, 2))
  expect_equal(nrow(cvtab <- cvfit$cvtab), 4L)
  expect_true(all(c("J", "lambda.min", "class", "sd", "nzero") %in%
                    names(cvtab)))
  # The stored full-data fit predicts without refitting.
  cls <- predict(cvfit, dat$x, type = "class")
  expect_lt(mean(drop(cls) != dat$y), 0.3)
  expect_equal(predict(cvfit, dat$x, s = "lambda.min"),
               predict(cvfit$wall.fit, dat$x, s = cvfit$lambda.min))
  beta <- coef(cvfit)
  expect_equal(nrow(beta), cvfit$wall.fit$nvars + 1L)
  expect_output(print(cvfit), "Selected: J =")
})

test_that("cv.wall shares the folds across the grid of J", {
  dat <- .wall_sim(n = 300)

  foldid <- rep_len(1:5, nrow(dat$x))
  cv1 <- cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8, foldid = foldid)
  cv2 <- cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8, foldid = foldid)

  expect_equal(cv1$foldid, foldid)
  expect_equal(cv1$cvtab, cv2$cvtab)  # deterministic given the folds
})

test_that("cv.wall supports the deviance and auc measures", {
  dat <- .wall_sim(n = 300)

  set.seed(2)
  cvd <- cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8, nfolds = 5,
                 type.measure = "deviance")
  expect_equal(cvd$type.measure, "deviance")
  expect_true(cvd$cvm.min > 0)

  set.seed(3)
  cva <- suppressWarnings(cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8,
                                  nfolds = 5, type.measure = "auc"))
  expect_equal(cva$type.measure, "auc")
  # AUC is maximized: the selected pair attains the largest cvm.min.
  expect_equal(cva$cvm.min, max(cva$cvtab$auc))
  expect_gt(cva$cvm.min, 0.5)
})

test_that("cv.wall builds a default grid of J and validates inputs", {
  dat <- .wall_sim()

  set.seed(4)
  cvfit <- cv.wall(dat$x, dat$y, filter.size = 8, nfolds = 5)
  expect_equal(cvfit$J, seq_len(max(2, ceiling(log2(nrow(dat$x))/2))))

  expect_error(cv.wall(dat$x, dat$y, J = 0:2, filter.size = 8),
               "larger than 'j0'")
  expect_error(cv.wall(dat$x, dat$y, J = 2:3, eps = c(0.1, 0.2)),
               "single value")
  expect_error(cv.wall(dat$x, dat$y, J = 2:3, foldid = c(1, 2)),
               "one entry per observation")
})

test_that("plot.cv.wall draws without errors", {
  dat <- .wall_sim()

  set.seed(5)
  cvfit <- cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8, nfolds = 5)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(cvfit))
  expect_invisible(plot(cvfit, se = FALSE, legend.pos = NULL))
})
