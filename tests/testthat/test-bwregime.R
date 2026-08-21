# Unit tests for the identification of regime switches of Motta and Montoril
# (2026b).

# Transcription of Devroye's sampler for PG(1, z), consuming the same variates
# in the same order as the compiled one, so that a logit chain can be followed
# draw by draw. It is the sampler of src/wav_polyagamma.c, line by line.
ref_pg1 <- function(z){

  trunc <- 0.64
  half <- 0.5*abs(z)
  fz <- 0.125*pi^2 + 0.5*half^2

  coef <- function(n, x){
    k <- (n + 0.5)*pi
    if(x > trunc) return(k*exp(-0.5*k*k*x))
    if(x > 0) return(exp(-1.5*(log(0.5*pi) + log(x)) + log(k)
                         - 2*(n + 0.5)^2/x))
    0
  }

  # Probability that the proposal is the truncated exponential one.
  b <- sqrt(1/trunc)*(trunc*half - 1)
  a <- -sqrt(1/trunc)*(trunc*half + 1)
  x0 <- log(fz) + fz*trunc
  p <- 1/(1 + 4/pi*(exp(x0 - half + pnorm(b, log.p = TRUE))
                    + exp(x0 + half + pnorm(a, log.p = TRUE))))

  # Inverse Gaussian truncated to (0, trunc).
  rtigauss <- function(){
    x <- trunc + 1
    if(half < 1/trunc){
      alpha <- 0
      repeat{
        if(runif(1) <= alpha) break
        e1 <- rexp(1); e2 <- rexp(1)
        while(e1*e1 > 2*e2/trunc){ e1 <- rexp(1); e2 <- rexp(1) }
        x <- 1 + e1*trunc
        x <- trunc/(x*x)
        alpha <- exp(-0.5*half*half*x)
      }
    }
    else{
      mu <- 1/half
      while(x > trunc){
        y <- rnorm(1)^2
        x <- mu + 0.5*mu*mu*y - 0.5*mu*sqrt(4*mu*y + (mu*y)^2)
        if(runif(1) > mu/(mu + x)) x <- mu*mu/x
      }
    }
    x
  }

  repeat{
    x <- if(runif(1) < p) trunc + rexp(1)/fz else rtigauss()

    s <- coef(0, x)
    u <- runif(1)*s
    n <- 0

    repeat{
      n <- n + 1
      if(n %% 2 == 1){
        s <- s - coef(n, x)
        if(u <= s) return(0.25*x)
      }
      else{
        s <- s + coef(n, x)
        if(u > s) break
      }
    }
  }
}

# The working response of one sweep, which is where the two links differ. The
# logit branch is the Polya-Gamma augmentation completed to a common scale by
# the orthogonal data augmentation, and it returns that scale.
ref_working <- function(link, fit, zt, tnorm){

  npad <- length(fit)

  if(link == "probit")
    return(list(v = tnorm(fit, zt == 1L, runif(npad)), sigma = 1))

  om <- vapply(fit, ref_pg1, 0)
  cc <- max(om)
  d <- cc - om

  list(v = (zt - 0.5 + d*fit + sqrt(d)*rnorm(npad))/cc, sigma = 1/sqrt(cc))
}


# Transcription of Algorithm 1 of the paper, drawing the same variates in the
# same order as the compiled sampler, so that the two chains must coincide.
ref_regime <- function(y, padidx, nchain, burn, lag, slab, prior, shr, init,
                       cprior = "gamma", wprior = "spikeslab",
                       link = "probit", family = "Daublets", filter.size = 20){

  n <- length(y)
  npad <- length(padidx)
  J <- round(log2(npad))
  mu <- init$mu
  tau2 <- init$tau2

  # The scale of the working response, which the link supplies at every sweep.
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

  # The auxiliary variable of the scale mixture that represents the half-t
  # prior, drawn from its own full conditional given the precision of its slot.
  halft_aux <- function(tau2k, k)
    rgamma(1, shape = (prior$df[k] + 1)/2,
           rate = prior$df[k]*tau2k + 1/prior$scale[k]^2)

  theta <- numeric(npad)
  fit <- numeric(npad)
  gam <- integer(npad)
  alpha <- rep(pnorm(0), n)
  aux <- c(0, 0)

  # The state the horseshoe prior carries between sweeps: the local scales of
  # the coefficients and the global scale of each level, both starting at one.
  lam2 <- rep(1, npad)
  tau2l <- rep(1, J)

  if(cprior != "gamma")
    aux <- vapply(1:2, function(k) halft_aux(init$tau2[k], k), 0)

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

    z <- allocate(alpha, mu, tau2)

    wrk <- ref_working(link, fit, z[padidx], tnorm)
    sigma <- wrk$sigma
    dstar <- wavedec(wrk$v, j0 = 0, family = family, filter.size = filter.size)

    theta[1L] <- dstar[1L] + sigma*rnorm(1)

    pil <- numeric(J)
    parl <- numeric(J)

    for(j in seq_len(J) - 1L){

      m <- 2^j
      idx <- (m + 1L):(2L*m)

      if(wprior == "horseshoe"){

        # The inverse gamma scale mixture of Makalic and Schmidt (2016): the
        # auxiliary of a half-Cauchy scale is drawn from the scale itself, so
        # only the scales are carried from the previous sweep.
        ax <- 1/rgamma(1, shape = 1, rate = 1/shr$scale^2 + 1/tau2l[j + 1L])
        tau2l[j + 1L] <- 1/rgamma(1, shape = (m + 1)/2,
                                  rate = 1/ax + sum(theta[idx]^2/lam2[idx])/
                                                (2*sigma^2))

        kap <- numeric(m)

        for(u in seq_along(idx)){
          k <- idx[u]
          ax <- 1/rgamma(1, shape = 1, rate = 1 + 1/lam2[k])
          lam2[k] <- 1/rgamma(1, shape = 1,
                              rate = 1/ax +
                                     theta[k]^2/(2*tau2l[j + 1L]*sigma^2))

          kap[u] <- lam2[k]*tau2l[j + 1L]/(1 + lam2[k]*tau2l[j + 1L])
          theta[k] <- kap[u]*dstar[k] + sigma*sqrt(kap[u])*rnorm(1)
          gam[k] <- 1L
        }

        pil[j + 1L] <- mean(kap)
        parl[j + 1L] <- tau2l[j + 1L]

        next
      }

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
    alpha <- if(link == "probit") pnorm(fit[seq_len(n)])
             else plogis(fit[seq_len(n)])

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
                    rate = rep(0.01, 2), df = rep(1, 2),
                    scale = rep(sd(y), 2)),
       init = list(mu = unname(quantile(y, c(0.25, 0.75))),
                   tau2 = rep(1/var(y), 2)),
       shr = list(zeta = 1, rho = 1, kappa = 1, xi = 100, cut = 0, scale = 1))

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


test_that("the logit link reproduces its transcription", {

  y <- mix_data()$y
  n <- length(y)
  s <- ref_settings(y)

  for(wp in c("spikeslab", "horseshoe")){

    set.seed(31)
    got <- bwregime(y, link = "logit", wavelet.prior = wp, nchain = 15,
                    burn = 5, lag = 2, plot = FALSE)

    set.seed(31)
    ref <- ref_regime(y, padidx = seq_len(n), nchain = 15, burn = 5, lag = 2,
                      slab = "laplace", prior = s$prior, shr = s$shr,
                      init = s$init, wprior = wp, link = "logit")

    expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
    expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
    expect_equal(unname(got$draws$alpha), ref$alpha, tolerance = 1e-10)
    expect_equal(unname(got$draws[[4L]]), ref$pi, tolerance = 1e-10)
    expect_equal(unname(got$draws[[5L]]), ref$slab, tolerance = 1e-10)
  }
})


test_that("the orthogonal augmentation targets the exact posterior", {

  # The algebra of the logit branch of ref_working(), with the weights held
  # fixed so that the posterior of the coefficients is available in closed
  # form. Together with the transcription test above, which ties the compiled
  # sampler to that reference variate by variate, this is what says that the
  # compiled logit chain has the right stationary distribution.
  set.seed(41)
  m <- 16

  X <- vapply(seq_len(m), function(k){
    e <- numeric(m); e[k] <- 1
    waverec(e, j0 = 0, family = "Daublets", filter.size = 8)
  }, numeric(m))

  expect_equal(crossprod(X), diag(m), tolerance = 1e-10)

  # Weights deliberately spread over two orders of magnitude, which is the case
  # the augmentation is least comfortable with.
  omega <- runif(m, 0.02, 1.2)
  yw <- rnorm(m)
  dvar <- runif(m, 0.1, 4)

  prec <- crossprod(X, omega*X) + diag(1/dvar)
  covar <- solve(prec)
  mean0 <- covar %*% crossprod(X, omega*yw)

  cc <- max(omega)
  kap <- dvar/(dvar + 1/cc)
  theta <- numeric(m)
  keep <- matrix(NA_real_, 20000, m)

  for(i in seq_len(20000)){
    eta <- as.vector(X %*% theta)
    d <- cc - omega
    v <- (omega*yw + d*eta + sqrt(d)*rnorm(m))/cc
    theta <- kap*wavedec(v, j0 = 0, family = "Daublets", filter.size = 8) +
             sqrt(kap/cc)*rnorm(m)
    keep[i, ] <- theta
  }

  # Within the Monte Carlo error of a chain of this length.
  expect_lt(max(abs(colMeans(keep) - mean0)/apply(keep, 2, sd)), 0.1)
  expect_lt(max(abs(apply(keep, 2, sd)/sqrt(diag(covar)) - 1)), 0.05)
})


test_that("the logit link estimates the same regimes as the probit one", {

  d <- mix_data(n = 256, seed = 11)

  set.seed(5)
  pr <- bwregime(d$y, nchain = 300, burn = 500, lag = 2, plot = FALSE)
  set.seed(5)
  lg <- bwregime(d$y, link = "logit", nchain = 300, burn = 500, lag = 2,
                 plot = FALSE)

  expect_equal(lg$link, "logit")
  expect_equal(pr$link, "probit")

  # The link is a modelling choice, not a different estimand: what the two must
  # agree on is the answer the method gives, which is where the weights cross
  # one half, and the components behind it. The weights themselves differ by
  # more than the chain error at this length, since the two links map the same
  # expansion through different functions.
  expect_lt(mean((lg$alpha - d$alpha)^2), 0.15)
  expect_gt(mean((lg$alpha > 0.5) == (pr$alpha > 0.5)), 0.9)
  expect_equal(lg$mu, pr$mu, tolerance = 0.1)

  expect_output(print(lg), "logit link")
})


test_that("the horseshoe prior reproduces its transcription", {

  y <- mix_data()$y
  n <- length(y)
  s <- ref_settings(y)

  set.seed(2024)
  got <- bwregime(y, wavelet.prior = "horseshoe", nchain = 20, burn = 10,
                  lag = 2, plot = FALSE)

  set.seed(2024)
  ref <- ref_regime(y, padidx = seq_len(n), nchain = 20, burn = 10, lag = 2,
                    slab = "laplace", prior = s$prior, shr = s$shr,
                    init = s$init, wprior = "horseshoe")

  expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
  expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
  expect_equal(unname(got$draws$alpha), ref$alpha, tolerance = 1e-10)
  expect_equal(unname(got$draws$kappa), ref$pi, tolerance = 1e-10)
  expect_equal(unname(got$draws$scale), ref$slab, tolerance = 1e-10)

  # The slab is not consulted, so choosing another one changes nothing.
  set.seed(2024)
  other <- bwregime(y, wavelet.prior = "horseshoe", slab = "gaussian",
                    nchain = 20, burn = 10, lag = 2, plot = FALSE)

  expect_equal(other$draws$alpha, got$draws$alpha)
})


test_that("the horseshoe prior shrinks without selecting", {

  y <- mix_data(n = 256)$y

  set.seed(7)
  hs <- bwregime(y, wavelet.prior = "horseshoe", nchain = 40, burn = 100,
                 lag = 2, plot = FALSE)
  set.seed(7)
  ss <- bwregime(y, nchain = 40, burn = 100, lag = 2, plot = FALSE)

  expect_equal(hs$wavelet.prior, "horseshoe")
  expect_equal(ss$wavelet.prior, "spikeslab")

  # Nothing is set to zero: the weights are shrinkage factors in the open unit
  # interval, and not the frequencies of an indicator.
  expect_true(all(hs$inclusion > 0 & hs$inclusion <= 1))
  expect_false(any(hs$inclusion[-1L]*40 == round(hs$inclusion[-1L]*40)))
  expect_true(all(ss$inclusion*40 == round(ss$inclusion*40)))
  expect_true(any(ss$inclusion < 1))

  # The scaling coefficient carries a diffuse prior under either prior, and is
  # never shrunk.
  expect_equal(unname(hs$inclusion[1L]), 1)
  expect_equal(unname(ss$inclusion[1L]), 1)

  # The level summaries are named after what they hold, and the mean shrinkage
  # weight of a level is the average of the weights of its coefficients.
  expect_named(hs$draws, c("mu", "tau2", "alpha", "kappa", "scale"))
  expect_named(ss$draws, c("mu", "tau2", "alpha", "pi", "slab"))
  expect_true(all(hs$draws$kappa > 0 & hs$draws$kappa < 1))
  expect_true(all(hs$draws$scale > 0))
  expect_equal(colnames(hs$draws$kappa), paste0("j", 0:(log2(256) - 1)))

  # The two summaries are averages of the same weights, over the coefficients
  # of a level and over the sweeps, so they must agree level by level.
  for(j in 0:(log2(256) - 1)){
    idx <- (2^j + 1L):(2^(j + 1L))
    expect_equal(mean(hs$inclusion[idx]), mean(hs$draws$kappa[, j + 1L]))
  }

  # The finest levels hold noise, and are shrunk harder than the coarsest ones.
  kap <- colMeans(hs$draws$kappa)
  expect_lt(kap[length(kap)], kap[1L])

  # The printout names the prior that was used.
  expect_output(print(hs), "horseshoe prior")
  expect_output(print(hs), "left unshrunk")
  expect_output(print(ss), "laplace slab")
})


test_that("the horseshoe prior recovers a two-regime weight function", {

  d <- mix_data(n = 256, seed = 11)

  set.seed(3)
  fit <- bwregime(d$y, wavelet.prior = "horseshoe", nchain = 200, burn = 500,
                  lag = 2, plot = FALSE)

  expect_lt(mean((fit$alpha - d$alpha)^2), 0.05)
  expect_lt(fit$mu[1L], fit$mu[2L])
  expect_true(mean(fit$alpha[1:128]) < mean(fit$alpha[129:256]))
})


test_that("the half-t prior of the component scales follows the scale mixture", {

  y <- mix_data()$y
  n <- length(y)
  s <- ref_settings(y)

  # The half-Cauchy prior is the half-t prior with one degree of freedom, and
  # both are drawn through the auxiliary variable of Wand et al. (2011).
  for(cprior in c("halfcauchy", "halft")){

    df <- if(cprior == "halfcauchy") c(1, 1) else c(3, 7)
    prior <- s$prior
    prior$df <- df

    set.seed(2024)
    got <- bwregime(y, scale.prior = cprior,
                    components = if(cprior == "halft") list(df = df)
                                 else list(),
                    nchain = 20, burn = 10, lag = 2, plot = FALSE)

    set.seed(2024)
    ref <- ref_regime(y, padidx = seq_len(n), nchain = 20, burn = 10, lag = 2,
                      slab = "laplace", cprior = cprior, prior = prior,
                      shr = s$shr, init = s$init)

    expect_equal(unname(got$draws$mu), ref$mu, tolerance = 1e-10)
    expect_equal(unname(got$draws$tau2), ref$tau2, tolerance = 1e-10)
    expect_equal(unname(got$draws$alpha), ref$alpha, tolerance = 1e-10)

    expect_equal(got$scale.prior, cprior)
    expect_equal(got$components$df, df)
    expect_equal(got$components$scale, rep(sd(y), 2))
    expect_true(all(got$draws$tau2 > 0))
  }

  # The default is the prior of the paper, which the option leaves untouched.
  set.seed(2024)
  gam <- bwregime(y, nchain = 20, burn = 10, lag = 2, plot = FALSE)

  expect_equal(gam$scale.prior, "gamma")
  expect_false(isTRUE(all.equal(gam$draws$tau2, got$draws$tau2)))
})


test_that("the half-t prior leaves the component variances freer", {

  # A short second regime, so that the upper component holds few observations
  # and the prior of its scale is felt. The gamma prior of the precision keeps
  # it away from the origin, whereas the half-t prior of the standard deviation
  # is flat there.
  set.seed(20)
  n <- 512
  a <- 0.02 + 0.9*((1:n)/n > 0.85)
  z <- rbinom(n, size = 1, prob = a)
  y <- z*rnorm(n, mean = 3, sd = 0.4) + (1 - z)*rnorm(n, mean = 0, sd = 0.4)

  set.seed(11)
  fit <- bwregime(y, scale.prior = "halfcauchy", nchain = 500, burn = 500,
                  lag = 2,
                  plot = FALSE)

  expect_equal(unname(fit$mu[1]), 0, tolerance = 0.15)
  expect_equal(unname(fit$mu[2]), 3, tolerance = 0.25)
  expect_equal(unname(1/sqrt(fit$tau2)), c(0.4, 0.4), tolerance = 0.2)

  # The two regimes are told apart, as they are under the prior of the paper.
  expect_true(mean(fit$alpha[1:400]) < 0.3)
  expect_true(mean(fit$alpha[450:512]) > 0.7)

  # The degrees of freedom of the half-t prior are the ones asked for, and a
  # large number of them approaches a half-normal prior, which is lighter
  # tailed and therefore shrinks the standard deviations more.
  set.seed(11)
  heavy <- bwregime(y, scale.prior = "halft", components = list(df = c(50, 50)),
                    nchain = 200, burn = 200, lag = 2, plot = FALSE)

  expect_equal(heavy$components$df, c(50, 50))
  expect_true(all(heavy$draws$tau2 > 0))
})


test_that("the prior of the component scales validates its settings", {

  y <- mix_data(n = 64)$y

  expect_error(bwregime(y, scale.prior = "inverse-gamma", plot = FALSE),
               "should be one of")
  expect_error(bwregime(y, components = list(df = c(0, 1)), plot = FALSE),
               "must be positive")
  expect_error(bwregime(y, components = list(scale = c(1, -1)), plot = FALSE),
               "must be positive")
  expect_error(bwregime(y, components = list(scale = 1), plot = FALSE),
               "two finite values")

  # The half-Cauchy prior has a single degree of freedom by definition, so
  # asking for another number of them is a contradiction and not a preference.
  expect_error(bwregime(y, scale.prior = "halfcauchy",
                        components = list(df = c(3, 3)), plot = FALSE),
               "fixed at one")
  expect_silent(bwregime(y, scale.prior = "halfcauchy",
                         components = list(df = c(1, 1)), nchain = 5, burn = 2,
                         lag = 1, plot = FALSE))

  # The printout names the prior in use.
  set.seed(3)
  fit <- bwregime(y, scale.prior = "halft", components = list(df = c(4, 4)),
                  nchain = 5,
                  burn = 2, lag = 1, plot = FALSE)

  expect_output(print(fit), "half-t prior of the standard deviations \\(df = 4, 4\\)")
  expect_output(print(bwregime(y, nchain = 5, burn = 2, lag = 1, plot = FALSE)),
                "gamma prior of the precisions")
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
  expect_error(bwregime(y, link = "cloglog", plot = FALSE), "should be one of")
  expect_error(bwregime(y, components = list(nonsense = 1), plot = FALSE),
               "Unknown component")
  expect_error(bwregime(y, components = list(var = c(-1, 1)), plot = FALSE),
               "must be positive")
  expect_error(bwregime(y, init = list(mu = c(3, 1)), plot = FALSE),
               "mu\\[1\\] < mu\\[2\\]")
  expect_error(bwregime(y, shrinkage = list(xi = 0), plot = FALSE),
               "single positive value")
  expect_error(bwregime(y, shrinkage = list(cut = 1), plot = FALSE),
               "in \\[0, 1\\)")
  expect_error(bwregime(y, verbose = -1, plot = FALSE), "non-negative integer")

  # An argument that the sampler does not have reaches the plot method through
  # '...', where the graphics engine would drop it without a word, and the fit
  # that came back would be the one of the defaults.
  expect_error(bwregime(y, nonsense = 1, plot = FALSE), "Unknown argument")
  expect_error(bwregime(y, cprior = "halft", plot = FALSE), "renamed to")
  expect_error(bwregime(y, prior = list(), plot = FALSE), "renamed to")
  expect_silent(bwregime(y, col = 2, lwd = 3, nchain = 5, burn = 0, lag = 1,
                   plot = FALSE))
})
