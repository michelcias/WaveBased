# Unit tests for the identification of regime switches of Motta and Montoril
# (2026b).

# Transcription of Algorithm 1 of the paper, drawing the same variates in the
# same order as the compiled sampler, so that the two chains must coincide.
ref_regime <- function(y, padidx, nchain, burn, lag, slab, prior, shr, init,
                       family = "Daublets", filter.size = 20){

  n <- length(y)
  npad <- length(padidx)
  J <- round(log2(npad))
  mu <- init$mu
  tau2 <- init$tau2
  sigma <- 1

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

  # The inversion of the truncated distribution function, on the log scale.
  tnorm <- function(m, pos, u){
    lp <- log1p(-u)
    ifelse(pos,
           m - sigma*qnorm(lp + pnorm(m/sigma, log.p = TRUE), log.p = TRUE),
           m + sigma*qnorm(lp + pnorm(-m/sigma, log.p = TRUE), log.p = TRUE))
  }

  logbf <- function(d, par){
    if(slab == "gaussian"){
      r <- par/sigma^2
      0.5*(d*d*r/(sigma^2*(1 + r)) - log1p(r))
    }
    else{
      zm <- d/sigma - par*sigma
      zp <- d/sigma + par*sigma
      a1 <- pnorm(zm, log.p = TRUE) - dnorm(zm, log = TRUE)
      a2 <- pnorm(zp, log.p = TRUE, lower.tail = FALSE) - dnorm(zp, log = TRUE)
      log(0.5*par*sigma) + pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))
    }
  }

  theta <- numeric(npad)
  fit <- numeric(npad)
  gam <- integer(npad)
  alpha <- rep(pnorm(0), n)

  out <- list(mu = matrix(NA_real_, nchain, 2),
              tau2 = matrix(NA_real_, nchain, 2),
              alpha = matrix(NA_real_, nchain, n),
              pi = matrix(NA_real_, nchain, J),
              slab = matrix(NA_real_, nchain, J))
  keep <- 0L

  z <- allocate(alpha, mu, tau2)

  for(i in seq_len(burn + lag*nchain)){

    for(k in 1:2){
      g <- k - 1L
      nk <- sum(z == g)
      Bk <- 1/(nk*tau2[k] + 1/prior$var[k])
      bk <- Bk*(sum(y[z == g])*tau2[k] + prior$mean[k]/prior$var[k])
      mu[k] <- rnorm(1, bk, sqrt(Bk))
      tau2[k] <- rgamma(1, shape = prior$shape[k] + nk/2,
                        rate = prior$rate[k] + sum((y[z == g] - mu[k])^2)/2)
    }

    if(mu[2L] < mu[1L]){
      mu <- mu[2:1]
      tau2 <- tau2[2:1]
    }

    z <- allocate(alpha, mu, tau2)

    lat <- tnorm(fit, z[padidx] == 1L, runif(npad))
    dstar <- wavedec(lat, j0 = 0, family = family, filter.size = filter.size)

    theta[1L] <- dstar[1L] + sigma*rnorm(1)

    pil <- numeric(J)
    parl <- numeric(J)

    for(j in seq_len(J) - 1L){

      m <- 2^j
      idx <- (m + 1L):(2L*m)
      n1 <- sum(gam[idx])

      pil[j + 1L] <- rbeta(1, shr$zeta + n1, shr$rho + m - n1)
      parl[j + 1L] <- if(slab == "gaussian")
          1/rgamma(1, shape = shr$kappa + m/2,
                   rate = shr$xi + sum(theta[idx]^2)/2)
        else
          rgamma(1, shape = shr$kappa + m, rate = shr$xi + sum(abs(theta[idx])))

      lodds <- log(pil[j + 1L]) - log1p(-pil[j + 1L])
      par <- parl[j + 1L]

      for(k in idx){
        pp <- 1/(1 + exp(-lodds - logbf(dstar[k], par)))
        if(pp < shr$cut)
          pp <- 0
        gam[k] <- as.integer(runif(1) < pp)

        if(slab == "gaussian"){
          lam <- par/(par + sigma^2)
          x <- lam*dstar[k] + sigma*sqrt(lam)*rnorm(1)
        }
        else{
          shift <- par*sigma^2
          lup <- -par*dstar[k] + pnorm((dstar[k] - shift)/sigma, log.p = TRUE)
          ldn <- par*dstar[k] + pnorm((dstar[k] + shift)/sigma, log.p = TRUE,
                                      lower.tail = FALSE)
          pos <- runif(1) < 1/(1 + exp(ldn - lup))
          x <- tnorm(if(pos) dstar[k] - shift else dstar[k] + shift, pos,
                     runif(1))
        }

        theta[k] <- if(gam[k] == 1L) x else 0
      }
    }

    fit <- waverec(theta, j0 = 0, family = family, filter.size = filter.size)
    alpha <- pnorm(fit[seq_len(n)])

    if(i > burn && (i - burn) %% lag == 0){
      keep <- keep + 1L
      out$mu[keep, ] <- mu
      out$tau2[keep, ] <- tau2
      out$alpha[keep, ] <- alpha
      out$pi[keep, ] <- pil
      out$slab[keep, ] <- parl
    }
  }

  out
}


# The defaults that bwregime() builds, and a small mixture to run it on.
ref_settings <- function(y)
  list(prior = list(mean = unname(quantile(y, c(0.25, 0.75))),
                    var = rep(var(y), 2), shape = rep(0.01, 2),
                    rate = rep(0.01, 2)),
       init = list(mu = unname(quantile(y, c(0.25, 0.75))),
                   tau2 = rep(1/var(y), 2)),
       shr = list(zeta = 1, rho = 1, kappa = 1, xi = 100, cut = 0))

mix_data <- function(n = 128, seed = 99){
  set.seed(seed)
  a <- 0.2 + 0.6*((1:n)/n > 0.5)
  z <- rbinom(n, size = 1, prob = a)
  list(y = z*rnorm(n, mean = 3, sd = 0.7) + (1 - z)*rnorm(n, sd = 0.7),
       alpha = a)
}


test_that("the compiled sampler reproduces a transcription of Algorithm 1", {

  y <- mix_data()$y
  n <- length(y)
  s <- ref_settings(y)

  for(slab in c("laplace", "gaussian")){

    set.seed(2024)
    got <- bwregime(y, slab = slab, nchain = 20, burn = 10, lag = 2,
                    plot = FALSE)

    set.seed(2024)
    ref <- ref_regime(y, padidx = seq_len(n), nchain = 20, burn = 10, lag = 2,
                      slab = slab, prior = s$prior, shr = s$shr, init = s$init)

    expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
    expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
    expect_equal(unname(got$draws$alpha), ref$alpha, tolerance = 1e-10)
    expect_equal(unname(got$draws$pi), ref$pi, tolerance = 1e-10)
    expect_equal(unname(got$draws$slab), ref$slab, tolerance = 1e-10)
  }
})


test_that("the reflected padding is the one of the reference implementation", {

  y <- mix_data(n = 100)$y
  n <- length(y)
  s <- ref_settings(y)

  # The tail reversal c(y, rev(y)[1:m]) of the Remark of Section 5.2.
  padidx <- c(seq_len(n), rev(seq_len(n))[seq_len(128 - n)])

  set.seed(31)
  got <- bwregime(y, nchain = 10, burn = 5, lag = 1, plot = FALSE)

  set.seed(31)
  ref <- ref_regime(y, padidx = padidx, nchain = 10, burn = 5, lag = 1,
                    slab = "laplace", prior = s$prior, shr = s$shr,
                    init = s$init)

  expect_equal(got$npad, 128)
  expect_equal(unname(got$draws$alpha), ref$alpha, tolerance = 1e-10)
})


test_that("the chains are reproducible and respect the model constraints", {

  y <- mix_data()$y

  set.seed(1)
  a <- bwregime(y, nchain = 30, burn = 20, lag = 2, plot = FALSE)
  set.seed(1)
  b <- bwregime(y, nchain = 30, burn = 20, lag = 2, plot = FALSE)

  expect_equal(a$draws, b$draws)

  # The identifying restriction is enforced at every sweep, the weights are
  # probabilities, and the parameters of the priors are positive.
  expect_true(all(a$draws$mu[, 1] <= a$draws$mu[, 2]))
  expect_true(all(a$draws$alpha >= 0 & a$draws$alpha <= 1))
  expect_true(all(a$draws$tau2 > 0))
  expect_true(all(a$draws$pi > 0 & a$draws$pi < 1))
  expect_true(all(a$draws$slab > 0))

  # The scaling coefficient is never shrunk, and the detail coefficients are.
  expect_equal(unname(a$inclusion[1]), 1)
  expect_true(all(a$inclusion >= 0 & a$inclusion <= 1))
  expect_true(mean(a$inclusion[-1]) < 0.5)
})


test_that("the summaries have the announced shapes and meanings", {

  y <- mix_data()$y
  n <- length(y)

  set.seed(4)
  fit <- bwregime(y, nchain = 40, burn = 20, lag = 2, level = 0.9,
                  plot = FALSE)

  expect_s3_class(fit, "bwregime")
  expect_equal(length(fit$alpha), n)
  expect_equal(dim(fit$alpha.hpd), c(n, 2L))
  expect_equal(dim(fit$draws$alpha), c(40L, n))
  expect_equal(dim(fit$draws$pi), c(40L, log2(n)))
  expect_equal(colnames(fit$draws$pi), paste0("j", 0:(log2(n) - 1)))
  expect_equal(rownames(fit$estimates), c("mu1", "tau2.1", "mu2", "tau2.2"))
  expect_equal(colnames(fit$estimates), c("median", "mean", "lower", "upper"))
  expect_equal(fit$slab, "laplace")
  expect_equal(fit$link, "probit")

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


test_that("the cut of the inclusion probabilities behaves as documented", {

  y <- mix_data(n = 256)$y

  # The default keeps every coefficient the model selects, as in the paper.
  set.seed(15)
  none <- bwregime(y, nchain = 50, burn = 200, lag = 2, plot = FALSE)
  set.seed(15)
  cut <- bwregime(y, nchain = 50, burn = 200, lag = 2,
                  shrinkage = list(cut = 0.05), plot = FALSE)

  expect_equal(none$shrinkage$cut, 0)

  # Without the cut the sparsity parameters drift upwards, so that many more
  # coefficients survive and the weights are pushed towards zero and one.
  expect_true(mean(none$inclusion[-1]) > mean(cut$inclusion[-1]))
  expect_true(mean(none$alpha > 0.99 | none$alpha < 0.01) >
              mean(cut$alpha > 0.99 | cut$alpha < 0.01))

  # A cut of one half is stronger still.
  set.seed(15)
  strong <- bwregime(y, nchain = 50, burn = 200, lag = 2,
                     shrinkage = list(cut = 0.5), plot = FALSE)
  expect_true(mean(strong$inclusion[-1]) <= mean(cut$inclusion[-1]))

  # The wider slab and the sparsity favouring prior of the Details reach the
  # same end without excluding anything by hand.
  set.seed(15)
  wide <- bwregime(y, nchain = 50, burn = 200, lag = 2,
                   shrinkage = list(xi = 1000), plot = FALSE)
  set.seed(15)
  sparse <- bwregime(y, nchain = 50, burn = 200, lag = 2,
                     shrinkage = list(rho = 50), plot = FALSE)

  expect_true(mean(wide$inclusion[-1]) < mean(none$inclusion[-1]))
  expect_true(mean(sparse$inclusion[-1]) < mean(none$inclusion[-1]))
})


test_that("the methods extract and interpolate the estimates", {

  y <- mix_data()$y
  n <- length(y)

  set.seed(6)
  fit <- bwregime(y, nchain = 30, burn = 20, lag = 2, plot = FALSE)

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
  expect_output(print(fit), "Regime switches through Bayesian wavelet")
})


test_that("the sampler recovers a two-regime weight function", {

  # A well separated mixture whose weights switch once, which the sampler
  # should recover along with the component parameters.
  set.seed(20)
  n <- 512
  a <- 0.1 + 0.8*((1:n)/n > 0.4)
  z <- rbinom(n, size = 1, prob = a)
  y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)

  for(slab in c("laplace", "gaussian")){

    fit <- bwregime(y, slab = slab, nchain = 500, burn = 500, lag = 2,
                    plot = FALSE)

    expect_equal(unname(fit$mu[1]), 0, tolerance = 0.15)
    expect_equal(unname(fit$mu[2]), 2, tolerance = 0.15)
    expect_equal(unname(fit$tau2[1]), 4, tolerance = 1.5)
    expect_equal(unname(fit$tau2[2]), 4, tolerance = 1.5)

    # The two regimes are told apart, and the estimates follow the truth. The
    # default excludes nothing, so the weights sit closer to zero and one than
    # the curve does, which is what keeps the tolerance below loose.
    expect_true(mean(fit$alpha[1:200]) < 0.3)
    expect_true(mean(fit$alpha[300:512]) > 0.7)
    expect_true(mean((fit$alpha - a)^2) < 0.06)

    # The credible intervals cover the truth.
    expect_true(fit$estimates["mu1", "lower"] <= 0 &&
                fit$estimates["mu1", "upper"] >= 0)
    expect_true(fit$estimates["mu2", "lower"] <= 2 &&
                fit$estimates["mu2", "upper"] >= 2)
  }
})


test_that("the padding schemes behave as documented", {

  y <- mix_data(n = 100)$y

  set.seed(8)
  ref <- bwregime(y, nchain = 10, burn = 5, lag = 1, plot = FALSE)
  set.seed(8)
  per <- bwregime(y, nchain = 10, burn = 5, lag = 1, padding = "periodic",
                  plot = FALSE)

  expect_equal(ref$padding, "reflect")
  expect_equal(per$npad, 128)
  expect_false(isTRUE(all.equal(ref$alpha, per$alpha)))

  expect_error(bwregime(y, padding = "none", plot = FALSE), "power of two")

  y2 <- mix_data(n = 128)$y
  set.seed(8)
  none <- bwregime(y2, nchain = 10, burn = 5, lag = 1, padding = "none",
                   plot = FALSE)
  set.seed(8)
  full <- bwregime(y2, nchain = 10, burn = 5, lag = 1, plot = FALSE)

  # With a dyadic sample size there is nothing to pad, so the schemes agree.
  expect_equal(none$draws, full$draws)
})


test_that("bwregime validates its arguments", {

  y <- mix_data(n = 64)$y

  expect_error(bwregime(y, x = 1:10, plot = FALSE), "same length")
  expect_error(bwregime(c(y[-1], NA), plot = FALSE), "finite")
  expect_error(bwregime(y, x = rev(seq_along(y)), plot = FALSE), "sorted")
  expect_error(bwregime(y, nchain = 1, plot = FALSE), "at least two")
  expect_error(bwregime(y, burn = -1, plot = FALSE), "non-negative")
  expect_error(bwregime(y, lag = 0, plot = FALSE), "positive")
  expect_error(bwregime(y, level = 1, plot = FALSE), "between 0 and 1")
  expect_error(bwregime(y, slab = "cauchy", plot = FALSE), "should be one of")
  expect_error(bwregime(y, link = "logit", plot = FALSE), "should be .*probit")
  expect_error(bwregime(y, prior = list(nonsense = 1), plot = FALSE),
               "Unknown component")
  expect_error(bwregime(y, prior = list(var = c(-1, 1)), plot = FALSE),
               "must be positive")
  expect_error(bwregime(y, init = list(mu = c(3, 1)), plot = FALSE),
               "mu\\[1\\] < mu\\[2\\]")
  expect_error(bwregime(y, shrinkage = list(xi = 0), plot = FALSE),
               "single positive value")
  expect_error(bwregime(y, shrinkage = list(cut = 1), plot = FALSE),
               "in \\[0, 1\\)")
  expect_error(bwregime(y, verbose = -1, plot = FALSE), "non-negative integer")
})
