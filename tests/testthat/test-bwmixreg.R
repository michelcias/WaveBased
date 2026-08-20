# Unit tests for the Bayesian wavelet-based mixture regression of Motta and
# Montoril (2026).

# Transcription of Algorithm 1 of the paper, drawing the same variates in the
# same order as the compiled sampler, so that the two chains must coincide.
ref_gibbs <- function(y, padidx, nchain, burn, lag, prior, init, thr,
                      family, filter.size, cprior = "gamma"){

  n <- length(y)
  mu <- init$mu
  tau2 <- init$tau2

  weights <- function(mu){
    w <- (y[padidx] - mu[1L])/(mu[2L] - mu[1L])
    f <- bayesthresh(w, j0 = thr$j0, alpha = thr$alpha, beta = thr$beta,
                     C1 = thr$C1, C2 = thr$C2, C1.start = thr$C1.start,
                     family = family, filter.size = filter.size)$fit
    pmin(1, pmax(0, f[seq_len(n)]))
  }

  allocate <- function(alpha, mu, tau2){
    u <- runif(n)
    s <- 1/sqrt(tau2)
    l0 <- log1p(-alpha) + dnorm(y, mu[1L], s[1L], log = TRUE)
    l1 <- log(alpha) + dnorm(y, mu[2L], s[2L], log = TRUE)
    z <- as.integer(u < 1/(1 + exp(l0 - l1)))
    z[alpha <= 0] <- 0L
    z[alpha >= 1] <- 1L
    z
  }

  # The auxiliary variable of the scale mixture that represents the half-t
  # prior, drawn from its own full conditional given the precision of its slot.
  halft_aux <- function(tau2k, k)
    rgamma(1, shape = (prior$df[k] + 1)/2,
           rate = prior$df[k]*tau2k + 1/prior$scale[k]^2)

  aux <- c(0, 0)

  if(cprior != "gamma")
    aux <- vapply(1:2, function(k) halft_aux(init$tau2[k], k), 0)

  alpha <- weights(mu)
  z <- allocate(alpha, mu, tau2)

  total <- burn + lag*nchain
  out <- list(mu = matrix(NA_real_, nchain, 2),
              tau2 = matrix(NA_real_, nchain, 2),
              alpha = matrix(NA_real_, nchain, n))
  keep <- 0L

  for(i in seq_len(total)){

    for(k in 1:2){
      g <- k - 1L
      nk <- sum(z == g)
      Bk <- 1/(nk*tau2[k] + 1/prior$var[k])
      bk <- Bk*(sum(y[z == g])*tau2[k] + prior$mean[k]/prior$var[k])
      mu[k] <- rnorm(1, bk, sqrt(Bk))

      # Under the half-t prior the gamma full conditional is the same one, with
      # the shape df/2 and the rate df*aux of the scale mixture.
      shape <- if(cprior == "gamma") prior$shape[k] else prior$df[k]/2
      rate <- if(cprior == "gamma") prior$rate[k] else prior$df[k]*aux[k]

      tau2[k] <- rgamma(1, shape = shape + nk/2,
                        rate = rate + sum((y[z == g] - mu[k])^2)/2)
    }

    if(mu[2L] < mu[1L]){
      mu <- mu[2:1]
      tau2 <- tau2[2:1]
    }

    # The auxiliary variables are refreshed after the labels are permuted, from
    # the precision that ended up in each slot.
    if(cprior != "gamma")
      aux <- vapply(1:2, function(k) halft_aux(tau2[k], k), 0)

    alpha <- weights(mu)
    z <- allocate(alpha, mu, tau2)

    if(i > burn && (i - burn) %% lag == 0){
      keep <- keep + 1L
      out$mu[keep, ] <- mu
      out$tau2[keep, ] <- tau2
      out$alpha[keep, ] <- alpha
    }
  }

  out
}


# A small mixture, shared by the tests below.
mix_data <- function(n = 128, seed = 99){
  set.seed(seed)
  a <- 0.2 + 0.6*((1:n)/n > 0.5)
  z <- rbinom(n, size = 1, prob = a)
  list(y = z*rnorm(n, mean = 3, sd = 0.7) + (1 - z)*rnorm(n, sd = 0.7),
       alpha = a)
}


test_that("the compiled sampler reproduces a transcription of Algorithm 1", {

  d <- mix_data()
  y <- d$y
  n <- length(y)

  # The same defaults that bwmixreg() builds.
  prior <- list(mean = unname(quantile(y, c(0.25, 0.75))),
                var = rep(var(y), 2), shape = rep(0.01, 2),
                rate = rep(0.01, 2), df = rep(1, 2), scale = rep(sd(y), 2))
  init <- list(mu = prior$mean, tau2 = rep(1/var(y), 2))
  thr <- list(j0 = 3, alpha = 0.5, beta = 1, C1 = NA, C2 = NA, C1.start = 100)

  set.seed(2024)
  got <- bwmixreg(y, nchain = 20, burn = 10, lag = 2, plot = FALSE)

  set.seed(2024)
  ref <- ref_gibbs(y, padidx = seq_len(n), nchain = 20, burn = 10, lag = 2,
                   prior = prior, init = init, thr = thr,
                   family = "Coiflets", filter.size = 30)

  expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
  expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
  expect_equal(got$draws$alpha, ref$alpha, tolerance = 1e-10)
})


test_that("the half-t prior of the component scales follows the scale mixture", {

  y <- mix_data()$y
  n <- length(y)

  init <- list(mu = unname(quantile(y, c(0.25, 0.75))), tau2 = rep(1/var(y), 2))
  thr <- list(j0 = 3, alpha = 0.5, beta = 1, C1 = NA, C2 = NA, C1.start = 100)

  # The half-Cauchy prior is the half-t prior with one degree of freedom, and
  # both are drawn through the auxiliary variable of Wand et al. (2011).
  for(cprior in c("halfcauchy", "halft")){

    df <- if(cprior == "halfcauchy") c(1, 1) else c(3, 7)
    prior <- list(mean = init$mu, var = rep(var(y), 2), shape = rep(0.01, 2),
                  rate = rep(0.01, 2), df = df, scale = rep(sd(y), 2))

    set.seed(2024)
    got <- bwmixreg(y, cprior = cprior,
                    prior = if(cprior == "halft") list(df = df) else list(),
                    nchain = 20, burn = 10, lag = 2, plot = FALSE)

    set.seed(2024)
    ref <- ref_gibbs(y, padidx = seq_len(n), nchain = 20, burn = 10, lag = 2,
                     prior = prior, init = init, thr = thr, family = "Coiflets",
                     filter.size = 30, cprior = cprior)

    expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
    expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
    expect_equal(got$draws$alpha, ref$alpha, tolerance = 1e-10)

    expect_equal(got$cprior, cprior)
    expect_equal(got$prior$df, df)
    expect_equal(got$prior$scale, rep(sd(y), 2))
    expect_true(all(got$draws$tau2 > 0))
  }

  # The default is the prior of the paper, which the option leaves untouched.
  set.seed(2024)
  gam <- bwmixreg(y, nchain = 20, burn = 10, lag = 2, plot = FALSE)

  expect_equal(gam$cprior, "gamma")
  expect_false(isTRUE(all.equal(gam$draws$tau2, got$draws$tau2)))
})


test_that("the half-t prior recovers the mixture as the gamma prior does", {

  # The weights are estimated from the component means alone, so the prior of
  # the scales reaches them only through the allocation variables, and the two
  # priors are expected to agree closely on a well separated mixture.
  set.seed(20)
  n <- 512
  a <- 0.1 + 0.8*((1:n)/n > 0.4)
  z <- rbinom(n, size = 1, prob = a)
  y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)

  set.seed(11)
  fit_hc <- bwmixreg(y, cprior = "halfcauchy", nchain = 300, burn = 300,
                     plot = FALSE)
  set.seed(11)
  fit_ht <- bwmixreg(y, cprior = "halft", prior = list(df = c(4, 4)),
                     nchain = 300, burn = 300, plot = FALSE)

  for(fit in list(fit_hc, fit_ht)){
    expect_equal(unname(fit$mu[1]), 0, tolerance = 0.15)
    expect_equal(unname(fit$mu[2]), 2, tolerance = 0.15)
    expect_equal(unname(1/sqrt(fit$tau2)), c(0.5, 0.5), tolerance = 0.15)

    expect_true(mean(fit$alpha[1:200]) < 0.3)
    expect_true(mean(fit$alpha[300:512]) > 0.7)
    expect_true(mean((fit$alpha - a)^2) < 0.06)
  }

  expect_equal(fit_ht$prior$df, c(4, 4))
})


test_that("the prior of the component scales validates its settings", {

  y <- mix_data(n = 64)$y

  expect_error(bwmixreg(y, cprior = "inverse-gamma", plot = FALSE),
               "should be one of")
  expect_error(bwmixreg(y, prior = list(df = c(0, 1)), plot = FALSE),
               "must be positive")
  expect_error(bwmixreg(y, prior = list(scale = c(1, -1)), plot = FALSE),
               "must be positive")
  expect_error(bwmixreg(y, prior = list(scale = 1), plot = FALSE),
               "two finite values")

  # The half-Cauchy prior has a single degree of freedom by definition, so
  # asking for another number of them is a contradiction and not a preference.
  expect_error(bwmixreg(y, cprior = "halfcauchy", prior = list(df = c(3, 3)),
                        plot = FALSE),
               "fixed at one")
  expect_silent(bwmixreg(y, cprior = "halfcauchy", prior = list(df = c(1, 1)),
                         nchain = 5, burn = 2, lag = 1, plot = FALSE))

  # The printout names the prior in use.
  set.seed(3)
  fit <- bwmixreg(y, cprior = "halft", prior = list(df = c(4, 4)), nchain = 5,
                  burn = 2, lag = 1, plot = FALSE)

  expect_output(print(fit), "half-t prior of the standard deviations \\(df = 4, 4\\)")
  expect_output(print(bwmixreg(y, nchain = 5, burn = 2, lag = 1, plot = FALSE)),
                "gamma prior of the precisions")
})


test_that("the reflected padding is the one of the reference implementation", {

  d <- mix_data(n = 100)
  y <- d$y
  n <- length(y)

  # The tail reversal c(y, rev(y)[1:m]) of the reference code.
  padidx <- c(seq_len(n), rev(seq_len(n))[seq_len(128 - n)])

  prior <- list(mean = unname(quantile(y, c(0.25, 0.75))),
                var = rep(var(y), 2), shape = rep(0.01, 2),
                rate = rep(0.01, 2), df = rep(1, 2), scale = rep(sd(y), 2))
  init <- list(mu = prior$mean, tau2 = rep(1/var(y), 2))
  thr <- list(j0 = 3, alpha = 0.5, beta = 1, C1 = NA, C2 = NA, C1.start = 100)

  set.seed(31)
  got <- bwmixreg(y, nchain = 10, burn = 5, lag = 1, plot = FALSE)

  set.seed(31)
  ref <- ref_gibbs(y, padidx = padidx, nchain = 10, burn = 5, lag = 1,
                   prior = prior, init = init, thr = thr,
                   family = "Coiflets", filter.size = 30)

  expect_equal(got$npad, 128)
  expect_equal(got$draws$alpha, ref$alpha, tolerance = 1e-10)
})


test_that("the chains are reproducible and respect the model constraints", {

  y <- mix_data()$y

  set.seed(1)
  a <- bwmixreg(y, nchain = 30, burn = 20, plot = FALSE)
  set.seed(1)
  b <- bwmixreg(y, nchain = 30, burn = 20, plot = FALSE)

  expect_equal(a$draws, b$draws)

  # The identifying restriction is enforced at every sweep, the weights are
  # probabilities, and the precisions are positive.
  expect_true(all(a$draws$mu[, 1] <= a$draws$mu[, 2]))
  expect_true(all(a$draws$alpha >= 0 & a$draws$alpha <= 1))
  expect_true(all(a$draws$tau2 > 0))
  expect_true(all(is.finite(a$draws$sigma)))
})


test_that("the summaries have the announced shapes and meanings", {

  y <- mix_data()$y
  n <- length(y)

  set.seed(4)
  fit <- bwmixreg(y, nchain = 40, burn = 20, level = 0.9, plot = FALSE)

  expect_s3_class(fit, "bwmixreg")
  expect_equal(length(fit$alpha), n)
  expect_equal(dim(fit$alpha.hpd), c(n, 2L))
  expect_equal(dim(fit$draws$alpha), c(40L, n))
  expect_equal(rownames(fit$estimates), c("mu1", "tau2.1", "mu2", "tau2.2"))
  expect_equal(colnames(fit$estimates), c("median", "mean", "lower", "upper"))

  expect_equal(fit$alpha, apply(fit$draws$alpha, 2, median))
  expect_equal(fit$alpha.mean, colMeans(fit$draws$alpha))
  expect_equal(unname(fit$alpha.hpd), unname(hpdi(fit$draws$alpha, 0.9)))
  expect_equal(unname(fit$estimates[, "median"]),
               c(median(fit$draws$mu[, 1]), median(fit$draws$tau2[, 1]),
                 median(fit$draws$mu[, 2]), median(fit$draws$tau2[, 2])))

  # The point estimates lie inside their intervals.
  expect_true(all(fit$estimates[, "lower"] <= fit$estimates[, "median"]))
  expect_true(all(fit$estimates[, "median"] <= fit$estimates[, "upper"]))
  expect_true(all(fit$alpha.hpd[, "lower"] <= fit$alpha.hpd[, "upper"]))

  # The allocation probabilities are means of Bernoulli draws.
  expect_true(all(fit$z >= 0 & fit$z <= 1))
  expect_true(all(fit$z*40 == round(fit$z*40)))
})


test_that("the methods extract and interpolate the estimates", {

  y <- mix_data()$y
  n <- length(y)

  set.seed(6)
  fit <- bwmixreg(y, nchain = 30, burn = 20, plot = FALSE)

  expect_equal(fitted(fit), fit$alpha)
  expect_equal(fitted(fit, estimate = "mean"), fit$alpha.mean)
  expect_equal(unname(coef(fit)), unname(fit$estimates[, "median"]))
  expect_equal(predict(fit), fit$alpha)

  expect_equal(dim(predict(fit, what = "all")), c(n, 3L))
  expect_equal(colnames(predict(fit, what = "all")),
               c("estimate", "lower", "upper"))

  # The observation times are returned untouched, and the intermediate ones
  # are interpolated between their neighbours.
  expect_equal(predict(fit, newx = fit$x[c(3, 10)]), fit$alpha[c(3, 10)])

  mid <- (fit$x[3] + fit$x[4])/2
  expect_equal(predict(fit, newx = mid), mean(fit$alpha[3:4]))

  # Points outside the range of the data are not extrapolated.
  expect_true(is.na(predict(fit, newx = 2)))

  expect_equal(predict(fit, newx = fit$x[5], what = "lower"),
               unname(fit$alpha.hpd[5, "lower"]))

  # The plot methods run without error, on a null device.
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fit))
  expect_silent(plot(fit, data = TRUE, band = FALSE, estimate = "mean"))
  expect_output(print(fit), "Bayesian wavelet-based mixture regression")
})


test_that("the sampler recovers a known constant weight", {

  # A well separated mixture with a constant weight, which the sampler should
  # estimate near its true value, along with the component parameters.
  set.seed(20)
  n <- 512
  z <- rbinom(n, size = 1, prob = 0.75)
  y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)

  fit <- bwmixreg(y, nchain = 500, burn = 500, plot = FALSE)

  expect_equal(unname(fit$mu[1]), 0, tolerance = 0.15)
  expect_equal(unname(fit$mu[2]), 2, tolerance = 0.15)
  expect_equal(unname(fit$tau2[1]), 4, tolerance = 1.5)
  expect_equal(unname(fit$tau2[2]), 4, tolerance = 1.5)
  expect_equal(mean(fit$alpha), 0.75, tolerance = 0.1)

  # The credible intervals cover the truth.
  expect_true(fit$estimates["mu1", "lower"] <= 0 &&
              fit$estimates["mu1", "upper"] >= 0)
  expect_true(fit$estimates["mu2", "lower"] <= 2 &&
              fit$estimates["mu2", "upper"] >= 2)
})


test_that("the padding schemes behave as documented", {

  y <- mix_data(n = 100)$y

  set.seed(8)
  ref <- bwmixreg(y, nchain = 10, burn = 5, plot = FALSE)
  set.seed(8)
  per <- bwmixreg(y, nchain = 10, burn = 5, padding = "periodic", plot = FALSE)

  expect_equal(ref$padding, "reflect")
  expect_equal(per$npad, 128)
  expect_false(isTRUE(all.equal(ref$alpha, per$alpha)))

  expect_error(bwmixreg(y, padding = "none", plot = FALSE), "power of two")

  y2 <- mix_data(n = 128)$y
  set.seed(8)
  none <- bwmixreg(y2, nchain = 10, burn = 5, padding = "none", plot = FALSE)
  set.seed(8)
  full <- bwmixreg(y2, nchain = 10, burn = 5, plot = FALSE)

  # With a dyadic sample size there is nothing to pad, so the schemes agree.
  expect_equal(none$draws, full$draws)
})


test_that("bwmixreg validates its arguments", {

  y <- mix_data(n = 64)$y

  expect_error(bwmixreg(y, x = 1:10, plot = FALSE), "same length")
  expect_error(bwmixreg(c(y[-1], NA), plot = FALSE), "finite")
  expect_error(bwmixreg(y, x = rev(seq_along(y)), plot = FALSE), "sorted")
  expect_error(bwmixreg(y, nchain = 1, plot = FALSE), "at least two")
  expect_error(bwmixreg(y, burn = -1, plot = FALSE), "non-negative")
  expect_error(bwmixreg(y, lag = 0, plot = FALSE), "positive")
  expect_error(bwmixreg(y, level = 1, plot = FALSE), "between 0 and 1")
  expect_error(bwmixreg(y, thresholding = list(j0 = 6), plot = FALSE),
               "smaller than 6")
  expect_error(bwmixreg(y, prior = list(nonsense = 1), plot = FALSE),
               "Unknown component")
  expect_error(bwmixreg(y, prior = list(var = c(-1, 1)), plot = FALSE),
               "must be positive")
  expect_error(bwmixreg(y, init = list(mu = c(3, 1)), plot = FALSE),
               "mu\\[1\\] < mu\\[2\\]")
  expect_error(bwmixreg(y, verbose = -1, plot = FALSE), "non-negative integer")
})


test_that("the pyramid routines shared with the sampler are unchanged", {

  # The Gibbs sampler calls the same periodic transform as wavedec/waverec,
  # which was factored out of them. A round trip must be exact.
  set.seed(2)
  for(n in c(64, 256)){
    x <- rnorm(n)
    for(j0 in c(0, 2, log2(n) - 1))
      expect_equal(waverec(wavedec(x, j0 = j0, family = "Coiflets",
                                   filter.size = 30),
                           j0 = j0, family = "Coiflets", filter.size = 30),
                   x, tolerance = 1e-10)
  }
})
