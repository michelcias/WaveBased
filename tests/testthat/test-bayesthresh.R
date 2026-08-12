# Unit tests for the Bayesian wavelet thresholding of Abramovich, Sapatinas and
# Silverman (1998), as implemented by bayesthresh().

# Transcription of the "BayesThresh" policy of wavethresh::threshold.wd, acting
# on a WaveBased decomposition (j0 = 0) of length n. It is the specification
# that the compiled routine must reproduce.
ref_bayesthresh <- function(wc, jlo = 3, alpha = 0.5, beta = 1, C1 = NA,
                            C2 = NA, C1.start = 100){

  n <- length(wc)
  J <- log2(n)
  levs <- jlo:(J - 1)
  nthresh <- length(levs)
  idx <- function(j) (2^j + 1):(2^(j + 1))

  noise.level <- mad(wc[idx(J - 1)])
  v <- 2^(-alpha*levs)
  nsignal <- rep(0, nthresh)
  sum2 <- rep(0, nthresh)

  d <- wc[(2^jlo + 1):n]
  thr <- sqrt(2*log(length(d)))*mad(d)

  for(i in seq_len(nthresh)){
    dj <- wc[idx(levs[i])]
    dun <- dj*(abs(dj) > thr)
    nsignal[i] <- sum(abs(dun) > 1e-10)
    if(nsignal[i] > 0)
      sum2[i] <- sum(dun[abs(dun) > 0]^2)
  }

  if(is.na(C1)){
    if(sum(nsignal) == 0){
      wc[(2^jlo + 1):n] <- 0
      return(list(coef = wc, sigma = noise.level, C1 = 0, C2 = 0))
    }
    fn <- function(C)
      sum(nsignal*(log(noise.level^2 + C^2*v) -
          2*pnorm(-noise.level*sqrt(2*log(n))/sqrt(noise.level^2 + C^2*v),
                  log.p = TRUE)) + sum2/(noise.level^2 + C^2*v))
    C1 <- optimize(fn, interval = c(0, 50*sqrt(C1.start)))$minimum^2
  }

  tau2 <- C1*v

  if(is.na(C2)){
    p <- 2*pnorm(-noise.level*sqrt(2*log(n))/sqrt(noise.level^2 + tau2))
    C2 <- if(beta == 1) sum(nsignal/p)/J
          else (1 - 2^(1 - beta))/(1 - 2^((1 - beta)*J))*sum(nsignal/p)
  }

  if(C1 == 0 || C2 == 0)
    wc[(2^jlo + 1):n] <- 0
  else{
    pr <- pmin(1, C2*2^(-beta*levs))
    rat <- tau2/(noise.level^2 + tau2)
    for(i in seq_len(nthresh)){
      dj <- wc[idx(levs[i])]
      w <- (1 - pr[i])/pr[i]/sqrt(noise.level^2*rat[i]/tau2[i])*
             exp(-rat[i]*dj^2/2/noise.level^2)
      z <- 0.5*(1 + pmin(w, 1))
      wc[idx(levs[i])] <- sign(dj)*pmax(0, rat[i]*abs(dj) -
                                          noise.level*sqrt(rat[i])*qnorm(z))
    }
  }

  list(coef = wc, sigma = noise.level, C1 = C1, C2 = C2)
}


test_that("BayesThresh reproduces the policy of wavethresh::threshold.wd", {

  set.seed(7)

  settings <- list(list(n = 256, fam = "Coiflets", fs = 30, j0 = 3),
                   list(n = 512, fam = "Symmlets", fs = 10, j0 = 4),
                   list(n = 128, fam = "Daublets", fs = 8,  j0 = 2),
                   list(n = 256, fam = "Coiflets", fs = 30, j0 = 0))

  for(st in settings){

    tt <- (1:st$n)/st$n
    y <- sin(4*pi*tt) + rnorm(st$n, sd = 0.3)

    got <- bayesthresh(y, j0 = st$j0, family = st$fam, filter.size = st$fs)
    ref <- ref_bayesthresh(wavedec(y, j0 = 0, family = st$fam,
                                   filter.size = st$fs), jlo = st$j0)

    # The estimation of C1 goes through a minimizer written to reproduce
    # optimize(), so the hyperparameters must agree as well.
    expect_equal(got$sigma, ref$sigma, tolerance = 1e-12)
    expect_equal(got$C1, ref$C1, tolerance = 1e-8)
    expect_equal(got$C2, ref$C2, tolerance = 1e-8)
    expect_equal(got$coefficients, ref$coef, tolerance = 1e-10)
    expect_equal(got$fit, waverec(ref$coef, j0 = 0, family = st$fam,
                                  filter.size = st$fs), tolerance = 1e-10)
  }
})


test_that("BayesThresh honours the prior parameters and fixed hyperparameters", {

  set.seed(13)
  n <- 256
  y <- sin(4*pi*(1:n)/n) + rnorm(n, sd = 0.3)
  wc <- wavedec(y, j0 = 0, family = "Coiflets", filter.size = 30)

  for(ab in list(c(1, 1), c(0.5, 2), c(2, 0.5), c(0, 0))){
    got <- bayesthresh(y, alpha = ab[1], beta = ab[2])
    ref <- ref_bayesthresh(wc, alpha = ab[1], beta = ab[2])
    expect_equal(got$coefficients, ref$coef, tolerance = 1e-10)
    expect_equal(got$C2, ref$C2, tolerance = 1e-8)
  }

  got <- bayesthresh(y, C1 = 2.5, C2 = 0.3)
  ref <- ref_bayesthresh(wc, C1 = 2.5, C2 = 0.3)
  expect_equal(got$coefficients, ref$coef, tolerance = 1e-12)

  # A null hyperparameter empties the thresholded levels.
  got <- bayesthresh(y, j0 = 3, C1 = 0, C2 = 0.3)
  expect_true(all(got$coefficients[(2^3 + 1):n] == 0))
  expect_equal(got$coefficients[1:(2^3)], wc[1:(2^3)], tolerance = 1e-12)
})


test_that("BayesThresh leaves the coarse levels untouched", {

  set.seed(17)
  n <- 256
  y <- rnorm(n)
  wc <- wavedec(y, j0 = 0, family = "Coiflets", filter.size = 30)

  for(j0 in c(0, 2, 3, 5)){
    got <- bayesthresh(y, j0 = j0)
    expect_equal(got$coefficients[seq_len(2^j0)], wc[seq_len(2^j0)],
                 tolerance = 1e-12)
  }
})


test_that("BayesThresh matches wavethresh itself, on the same coefficients", {

  skip_if_not_installed("wavethresh")

  set.seed(11)
  n <- 256
  y <- sin(4*pi*(1:n)/n) + rnorm(n, sd = 0.3)

  w <- wavethresh::wd(y, filter.number = 5, family = "Coiflets")

  # The coefficients of the wd object, laid out in the ordering of WaveBased.
  as_wb <- function(obj){
    J <- wavethresh::nlevelsWT(obj)
    out <- numeric(2^J)
    out[1] <- wavethresh::accessC(obj, level = 0)
    for(j in 0:(J - 1))
      out[(2^j + 1):(2^(j + 1))] <- wavethresh::accessD(obj, level = j)
    out
  }

  # A signal carrying exactly those coefficients, so that both implementations
  # threshold the same numbers despite the different filter conventions.
  wc <- as_wb(w)
  x <- waverec(wc, j0 = 0, family = "Coiflets", filter.size = 30)
  expect_equal(wavedec(x, j0 = 0, family = "Coiflets", filter.size = 30), wc,
               tolerance = 1e-10)

  got <- bayesthresh(x, j0 = 3, family = "Coiflets", filter.size = 30)
  ref <- wavethresh::threshold(w, levels = 3:(wavethresh::nlevelsWT(w) - 1),
                               type = "soft", policy = "BayesThresh",
                               by.level = FALSE, dev = wavethresh::madmad,
                               boundary = FALSE, alpha = 0.5, beta = 1,
                               C1 = NA, C2 = NA, C1.start = 100)

  expect_equal(got$coefficients, as_wb(ref), tolerance = 1e-10)
})


test_that("bayesthresh validates its arguments", {

  y <- rnorm(256)

  expect_error(bayesthresh(rnorm(100)), "power of two")
  expect_error(bayesthresh(y, j0 = 8), "smaller than log2")
  expect_error(bayesthresh(y, j0 = -1), "non-negative")
  expect_error(bayesthresh(y, alpha = -1), "non-negative")
  expect_error(bayesthresh(y, C1 = -1), "non-negative")
  expect_error(bayesthresh(c(y[-1], NA)), "finite")
  expect_error(bayesthresh(y, family = "Nonsense"), "Unknown family")
})
