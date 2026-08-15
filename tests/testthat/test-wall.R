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
  expect_equal(plot(cvfit), cvfit)
})

test_that("plot.cv.wall delegates the path and the components to plot.wall", {
  dat <- .wall_sim()

  set.seed(5)
  cvfit <- cv.wall(dat$x, dat$y, J = 2:3, filter.size = 8, nfolds = 5)

  pdf(NULL)
  on.exit(dev.off())

  # The selected lambda is passed to plot.wall, on the fit of the selected J.
  expect_equal(plot(cvfit, type = "path"),
               plot(cvfit$wall.fit, s = cvfit$lambda.min))
  expect_equal(plot(cvfit, type = "components")$s, cvfit$lambda.min)
  expect_equal(plot(cvfit, type = "components", s = "lambda.1se")$s,
               cvfit$lambda.1se)
  expect_equal(plot(cvfit, type = "components", s = 0.05)$s, 0.05)
  expect_equal(plot(cvfit, type = "network")$s, cvfit$lambda.min)

  # Further arguments of plot.wall are honored.
  cp <- plot(cvfit, type = "components", which = "X1")
  expect_equal(cp$which, 1L)
  expect_equal(nrow(cp$components), 256L)
  expect_length(plot(cvfit, type = "path", max.vars = 1)$which, 1L)
})

test_that("plot.wall summarizes the coefficient path of each covariate", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  beta <- as.matrix(fit$glmnet.fit$beta)

  pdf(NULL)
  on.exit(dev.off())

  p <- plot(fit)
  expect_invisible(plot(fit))
  expect_equal(p$type, "path")
  expect_equal(dim(p$norms), c(2L, length(fit$lambda)))
  expect_equal(rownames(p$norms), c("X1", "X2"))
  # Each row is the norm of the block of coefficients of one covariate.
  expect_equal(p$norms["X1", ], sqrt(colSums(beta[1:7, ]^2)),
               ignore_attr = TRUE)
  expect_equal(p$norms["X2", ], sqrt(colSums(beta[8:14, ]^2)),
               ignore_attr = TRUE)

  p1 <- plot(fit, norm = "l1")
  expect_equal(p1$norms["X2", ], colSums(abs(beta[8:14, ])),
               ignore_attr = TRUE)

  # The covariates are ranked by importance, and 'which' overrides it.
  expect_equal(p$which, order(rowMeans(p$norms), decreasing = TRUE))
  expect_equal(plot(fit, which = "X2", legend.pos = NULL)$which, 2L)
  expect_equal(plot(fit, which = c(FALSE, TRUE))$which, 2L)
  expect_equal(plot(fit, max.vars = 1)$which, p$which[1L])

  expect_error(plot(fit, which = "X9"), "Unknown covariate")
  expect_error(plot(fit, which = 5), "must index")
  expect_error(plot(fit, which = c(TRUE, TRUE, FALSE)), "one entry per")
  expect_error(plot(fit, max.vars = 0), "single positive value")
  expect_error(plot(fit, s = c(0.1, 0.2)), "single non-negative value")
})

test_that("plot.wall draws the fitted components of the log-odds", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  s <- 0.02

  pdf(NULL)
  on.exit(dev.off())

  cp <- plot(fit, type = "components", s = s, n.grid = 2048, center = FALSE)
  expect_equal(cp$type, "components")
  expect_equal(cp$s, s)
  expect_equal(dim(cp$components), c(2048L, 2L))
  expect_equal(colnames(cp$components), c("X1", "X2")[cp$which])
  # The grid covers the range of the covariates observed in the fit.
  expect_equal(apply(cp$x, 2L, range), apply(dat$x[, cp$which], 2L, range),
               tolerance = 1e-3, ignore_attr = TRUE)

  # The components add up to the fitted log-odds, up to the intercept: this
  # checks the alignment between the blocks of coefficients, the rescaling
  # of each covariate and the basis evaluated on the grid.
  b0 <- as.numeric(coef(fit, s = s)[1L, 1L])
  eta <- as.numeric(predict(fit, dat$x[1:20, ], s = s))
  fl <- vapply(seq_along(cp$which),
               function(i) approx(cp$x[, i], cp$components[, i],
                                  xout = dat$x[1:20, cp$which[i]])$y,
               numeric(20L))
  expect_equal(b0 + rowSums(fl), eta, tolerance = 1e-4)

  # Centering only shifts each component by its mean over the grid.
  cpc <- plot(fit, type = "components", s = s, n.grid = 2048)
  expect_equal(colMeans(cpc$components), c(X1 = 0, X2 = 0)[cp$which],
               ignore_attr = TRUE)
  expect_equal(sweep(cp$components, 2L, colMeans(cp$components)),
               cpc$components)

  # Without 's', the least penalized fit of the path is displayed.
  expect_equal(plot(fit, type = "components")$s, min(fit$lambda))
  expect_error(plot(fit, type = "components", n.grid = 1),
               "larger than one")
})

test_that("plot.wall draws the classifier as a layered network", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  s <- 0.02

  pdf(NULL)
  on.exit(dev.off())

  np <- plot(fit, type = "network", s = s)
  expect_invisible(plot(fit, type = "network", s = s))
  expect_equal(np$type, "network")
  expect_equal(np$s, s)
  expect_equal(plot(fit, type = "network")$s, min(fit$lambda))

  # One node per resolution level of each covariate, holding the block of
  # coefficients of that level: with j0 = 0 and the periodized basis, the
  # blocks of a covariate are W0 (1 coefficient), W1 (2) and W2 (4).
  b <- as.numeric(coef(fit, s = s)[, 1L])
  idx <- list(W0 = 1L, W1 = 2:3, W2 = 4:7)
  for(l in seq_len(2L)){
    bl <- b[1L + (l - 1L)*7L + seq_len(7L)]
    nz <- vapply(idx, function(i) sum(bl[i] != 0), 0L)
    nd <- np$nodes[np$nodes$var == c("X1", "X2")[l], ]
    expect_equal(nd$label, names(idx)[nz > 0])
    expect_equal(nd$level, (0:2)[nz > 0])
    expect_equal(nd$nnz, unname(nz[nz > 0]))
    expect_equal(nd$norm,
                 unname(vapply(idx, function(i) sqrt(sum(bl[i]^2)), 0)[nz > 0]))
  }

  # Without an observation, the display is structural: the components are
  # summarized by the norm of their coefficients and there is nothing to
  # propagate through the last two layers.
  expect_equal(np$intercept, b[1L])
  expect_true(all(is.na(np$nodes$value)))
  expect_true(is.na(np$h) && is.na(np$prob))
  expect_equal(np$f, vapply(split(np$nodes$norm, np$nodes$var),
                            function(z) sqrt(sum(z^2)), 0)[names(np$f)])

  # Every resolution level is kept with nonzero = FALSE, even the ones the
  # LASSO emptied.
  expect_equal(nrow(plot(fit, type = "network", s = s, nonzero = FALSE)$nodes),
               2L*3L)
})

test_that("plot.wall propagates one observation through the network", {
  dat <- .wall_sim()
  fit <- wall(dat$x, dat$y, J = 3, filter.size = 8)
  s <- 0.02
  x0 <- dat$x[1L, ]

  pdf(NULL)
  on.exit(dev.off())

  q <- plot(fit, type = "network", s = s, newx = x0)

  # The nodes hold the forward pass of the observation: the values of the
  # levels of a covariate add up to its component, the components and the
  # intercept add up to the log-odds, and the logistic link closes the
  # network with the predicted probability.
  expect_equal(vapply(split(q$nodes$value, q$nodes$var), sum, 0)[names(q$f)],
               q$f)
  expect_equal(sum(q$f) + q$intercept, q$h)
  expect_equal(q$h, as.numeric(predict(fit, dat$x[1L, , drop = FALSE], s = s)))
  expect_equal(q$prob, as.numeric(predict(fit, dat$x[1L, , drop = FALSE],
                                          s = s, type = "response")))

  # A vector, a one-row matrix and a one-row data frame are the same
  # observation. The structure of the graph does not depend on it.
  expect_equal(q, plot(fit, type = "network", s = s,
                       newx = dat$x[1L, , drop = FALSE]))
  expect_equal(q, plot(fit, type = "network", s = s,
                       newx = as.data.frame(dat$x)[1L, ]))
  expect_equal(q$nodes[setdiff(names(q$nodes), "value")],
               plot(fit, type = "network", s = s)$nodes[
                 setdiff(names(q$nodes), "value")])

  expect_error(plot(fit, type = "network", s = s, newx = dat$x[1:2, ]),
               "single observation")
  expect_error(plot(fit, type = "network", s = s, newx = dat$x[1L, 1L]),
               "must have 2 column")
})

test_that("plot.wall keeps the display readable with many covariates", {
  set.seed(6)
  n <- 200
  x <- matrix(runif(20*n), n, 20)
  y <- rbinom(n, 1, 1/(1 + exp(-(8*(x[, 3] - 0.5)))))
  fit <- wall(x, y, J = 2, filter.size = 8)

  pdf(NULL)
  on.exit(dev.off())

  # The path shows every covariate, but only 'max.vars' of them are named.
  p <- plot(fit, max.vars = 4)
  expect_equal(nrow(p$norms), 20L)
  expect_length(p$which, 4L)

  # The components are restricted to the ones that are worth a panel, with a
  # message reporting the omitted covariates.
  expect_message(cp <- plot(fit, type = "components", s = 0.05, max.vars = 2),
                 "Showing 2 of 20 covariates")
  expect_length(cp$which, 2L)
  expect_true(all(apply(cp$components, 2L, function(f) diff(range(f))) > 0))

  # Covariates with a zero component are dropped, unless asked otherwise.
  s0 <- max(fit$lambda[colSums(as.matrix(fit$glmnet.fit$beta) != 0) > 0])
  expect_message(cp1 <- plot(fit, type = "components", s = s0),
                 "1 with a nonzero component")
  expect_length(cp1$which, 1L)
  expect_length(plot(fit, type = "components", s = s0, nonzero = FALSE,
                     max.vars = Inf)$which, 20L)
  expect_error(plot(fit, type = "components", s = max(fit$lambda)),
               "identically zero")

  # An explicit selection is honored as given, without a message.
  expect_silent(cp2 <- plot(fit, type = "components", s = 0.05,
                            which = c("X5", "X3")))
  expect_equal(cp2$which, c(5L, 3L))
  expect_equal(colnames(cp2$components), c("X5", "X3"))

  # The network shares that selection, and displays no node for the
  # covariates whose coefficients were all set to zero.
  expect_message(np <- plot(fit, type = "network", s = 0.05, max.vars = 2),
                 "Showing 2 of 20 covariates")
  expect_equal(np$which, cp$which)
  expect_silent(np2 <- plot(fit, type = "network", s = 0.2, which = c(5L, 3L),
                            legend.pos = NULL))
  expect_equal(np2$which, c(5L, 3L))
  expect_equal(sort(unique(np2$nodes$var)), sort(names(np2$f)[np2$f > 0]))

  # By default, the components fill a 3 x 3 grid of panels, laid out in
  # three columns when there are fewer of them. The layout in force is read
  # from inside the first panel, before it is restored on exit.
  e <- new.env()
  grab <- function(...)
    plot(fit, type = "components", s = 0.05, nonzero = FALSE, ...,
         panel.first = assign("layout", par("mfrow"), envir = e))

  expect_length(suppressMessages(grab())$which, 9L)
  expect_equal(e$layout, c(3, 3))
  suppressMessages(grab(max.vars = 4))
  expect_equal(e$layout, c(2, 3))
  suppressMessages(grab(max.vars = 2))
  expect_equal(e$layout, c(1, 2))
  suppressMessages(grab(max.vars = 4, mfrow = c(4, 1)))
  expect_equal(e$layout, c(4, 1))
  expect_error(grab(mfrow = 3), "'mfrow' must be")
})
