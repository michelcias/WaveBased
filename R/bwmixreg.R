#' @title Bayesian wavelet-based mixture regression
#'
#' @description Estimates jointly the dynamic mixture weights and the component
#' parameters of a two-component Gaussian mixture, by the Gibbs sampling
#' algorithm of Motta and Montoril (2026). It is the Bayesian counterpart of
#' \command{\link{wmixreg}}, which takes the component means as known.
#'
#' @param y Vector with the observed mixture. At each index the process behaves
#'   as \eqn{N(\mu_2, \tau_2^{-2})} with probability \eqn{\alpha_t}, and as
#'   \eqn{N(\mu_1, \tau_1^{-2})} otherwise, with \eqn{\mu_1 < \mu_2}.
#' @param x Vector of observation times, of the same length as \code{y}, used
#'   only to place the estimates on an axis and to interpolate them in
#'   \command{\link{predict.bwmixreg}}. It defaults to the equally spaced grid
#'   \eqn{i/n, i = 1, \ldots, n}. The discrete wavelet transform requires
#'   equally spaced samples, so \code{x} is assumed to be sorted and equally
#'   spaced. See Details.
#' @param nchain Number of draws to retain.
#' @param burn Number of sweeps discarded as burn-in.
#' @param lag Every \code{lag}-th sweep after the burn-in is retained, which
#'   reduces the autocorrelation of the chains. The sampler therefore performs
#'   \code{burn + lag*nchain} sweeps.
#' @param prior A list with the hyperparameters of the priors of the component
#'   parameters, any of which may be omitted. The entries are \code{mean} and
#'   \code{var}, the means and the variances of the normal priors of
#'   \eqn{\mu_1} and \eqn{\mu_2}, and \code{shape} and \code{rate}, the
#'   parameters of the gamma priors of \eqn{\tau_1^2} and \eqn{\tau_2^2}, each
#'   of them of length two. The defaults are the choices of Section 5 of the
#'   paper, described in Details.
#' @param init A list with the initial values \code{mu} and \code{tau2} of the
#'   chains, each of length two. The defaults are the means of the priors of
#'   \eqn{\mu_k}, and the reciprocal of the sample variance for both
#'   \eqn{\tau_k^2}.
#' @param thresholding A list with the settings of the Bayesian wavelet
#'   thresholding applied at each sweep, any of which may be omitted. The
#'   entries are \code{j0}, \code{alpha}, \code{beta}, \code{C1}, \code{C2} and
#'   \code{C1.start}, exactly as in \command{\link{bayesthresh}}, whose defaults
#'   they share. The names \code{alpha} and \code{beta} refer to the prior of
#'   the wavelet coefficients, and not to the mixture weights \eqn{\alpha_t}.
#' @param family The family of wavelets to use, as in \command{\link{wavedec}}.
#'   The default follows the paper, which uses a Coiflet basis.
#' @param filter.size The size of the wavelet filter. The default is the 30-tap
#'   Coiflet filter of the reference implementation of the paper, which
#'   corresponds to \code{filter.number = 5} in the notation of \pkg{wavethresh}
#'   (the filters of that family have length six times the filter number).
#' @param wavelet.filter Use this to provide your own filter, together with
#'   \code{family = "Own"}.
#' @param padding How the sample is extended to a dyadic length, which the
#'   discrete wavelet transform requires. Either \code{"reflect"} (the default,
#'   a mirror reflection of the tail of the sample, as in the reference
#'   implementation), \code{"periodic"} (a periodic extension) or \code{"none"}
#'   (which requires the sample size to be a power of two).
#' @param level The probability mass of the highest posterior density intervals
#'   reported.
#' @param verbose Logical, or a positive integer. Should the progress of the
#'   sampler be reported? An integer sets the reporting interval, in sweeps,
#'   whereas \code{TRUE} reports every 1000 sweeps.
#' @param plot Logical. Should the estimated mixture weights be plotted?
#' @param ... Further arguments passed to \command{\link{plot.bwmixreg}}.
#'
#' @details
#' The model of expression (2) of the paper is
#' \deqn{y_t = (1 - z_t)x_{1t} + z_t x_{2t}, \qquad
#'       x_{kt} \mid \mu_k, \tau_k^2 \sim N(\mu_k, \tau_k^{-2}), \qquad
#'       z_t \mid \alpha_t \sim \mathrm{Bern}(\alpha_t),}
#' where the \eqn{z_t} are allocation variables and the mixture weights
#' \eqn{\alpha_t} vary with the index \eqn{t}. Note that \eqn{\tau_k^2} is a
#' precision, so that the variance of the \eqn{k}-th component is
#' \eqn{\tau_k^{-2}}.
#'
#' Each sweep of the Gibbs sampler (Algorithm 1 of the paper) draws the
#' component parameters from the conjugate full conditionals
#' \deqn{\mu_k \mid \cdots \sim N(B_k b_k, B_k), \qquad
#'       \tau_k^2 \mid \cdots \sim \Gamma(v_k, V_k),}
#' permutes the labels of the pairs \eqn{(\mu_k, \tau_k^2)} whenever the draw
#' violates the identifying restriction \eqn{\mu_1 < \mu_2}, and draws the
#' allocation variables from their Bernoulli full conditionals. The mixture
#' weights are then rebuilt by turning the mixture into a nonparametric
#' regression: the transformed responses
#' \deqn{m_t = (y_t - \mu_1)/(\mu_2 - \mu_1)}
#' satisfy \eqn{E(m_t) = \alpha_t}, so \eqn{\alpha_t} is estimated by the
#' Bayesian wavelet regularization of \command{\link{bayesthresh}} applied to
#' \eqn{m}, truncated to \eqn{[0, 1]} because it is a probability.
#'
#' The default priors are the ones used in Section 5 of the paper, namely
#' \eqn{\mu_1 \sim N(q_1, s^2)} and \eqn{\mu_2 \sim N(q_3, s^2)}, where
#' \eqn{q_1} and \eqn{q_3} are the quartiles and \eqn{s^2} the sample variance
#' of the data, and \eqn{\tau_k^2 \sim \Gamma(0.01, 0.01)}. Following the
#' paper, the point estimates are posterior medians, which are optimal under
#' the absolute loss, and they are reported together with highest posterior
#' density intervals computed by \command{\link{hpdi}}.
#'
#' The whole sampler runs in compiled code, including the wavelet transforms
#' and the thresholding, so that no memory is allocated during the sweeps. It
#' uses the random number generator of R, and the chains are therefore
#' reproducible through \command{\link{set.seed}}.
#'
#' Two remarks on the assumptions. First, the discrete wavelet transform
#' requires equally spaced observations, and the sample size to be a power of
#' two; the latter is dealt with by \code{padding}, and the former is assumed.
#' Workarounds for non-equally spaced designs are discussed in Section 4.2 of
#' the paper and are not implemented here. Second, the estimates of
#' \eqn{\alpha_t} live on the grid of the observations, and not on a wavelet
#' expansion, so \command{\link{predict.bwmixreg}} interpolates them, unlike
#' the corresponding method of \command{\link{wmixreg}}.
#'
#' This function was implemented based on the reference R code of the authors,
#' available at
#' \url{https://github.com/flaviamotta/Bayesian-wavelet-mixture-regression},
#' which uses the packages \pkg{wavethresh} and \pkg{changepoint}. Two
#' deliberate differences are worth recording: the allocation variables are
#' drawn on the logarithmic scale, which keeps the full conditional (4) defined
#' when a weight reaches zero or one or when both normal densities underflow;
#' and the retained sweeps are those of indices
#' \code{burn + lag}, \code{burn + 2*lag}, and so on, rather than being offset
#' by the initial values. Moreover, the wavelet filters of this package and
#' those of \pkg{wavethresh} follow different conventions, so that the chains
#' produced here are statistically equivalent to, but not numerically identical
#' with, the ones of the reference implementation.
#'
#' @return An object of class \code{"bwmixreg"}, that is, a list with components
#' \item{call}{The matched call.}
#' \item{alpha}{The posterior medians of the mixture weights, of length
#'   \eqn{n}.}
#' \item{alpha.mean}{Their posterior means.}
#' \item{alpha.hpd}{Matrix with the highest posterior density intervals of the
#'   mixture weights, with one row per observation.}
#' \item{estimates}{Matrix with the posterior medians, means and highest
#'   posterior density intervals of \eqn{\mu_1}, \eqn{\tau_1^2}, \eqn{\mu_2}
#'   and \eqn{\tau_2^2}.}
#' \item{mu, tau2}{The posterior medians of the component parameters.}
#' \item{z}{The posterior means of the allocation variables, that is, the
#'   probability that each observation belongs to the second component.}
#' \item{draws}{List with the retained draws of \code{mu} and \code{tau2}
#'   (\code{nchain} by 2 matrices), of \code{alpha} (a \code{nchain} by
#'   \eqn{n} matrix, whose columns are ready for \command{\link{hpdi}}), and
#'   with the noise level \code{sigma} and the hyperparameters \code{C1} and
#'   \code{C2} selected by the thresholding at each retained sweep.}
#' \item{x, y}{The observation times and the observed mixture.}
#' \item{nobs, npad, padding}{The sample size, the padded sample size and the
#'   padding scheme.}
#' \item{nchain, burn, lag, level}{The settings of the sampler.}
#' \item{prior, init, thresholding}{The priors, the initial values and the
#'   thresholding settings actually used.}
#' \item{family, filter.size, wavelet.filter}{The wavelet basis used.}
#'
#' @references
#' Abramovich, F., Sapatinas, T. and Silverman, B. W. (1998). Wavelet
#' thresholding via a Bayesian approach. \emph{Journal of the Royal Statistical
#' Society Series B}, 60(4), 725--749,
#' \doi{10.1111/1467-9868.00151}.
#'
#' Lai, W. R., Johnson, M. D., Kucherlapati, R. and Park, P. J. (2005).
#' Comparative analysis of algorithms for identifying amplifications and
#' deletions in array CGH data. \emph{Bioinformatics}, 21(19), 3763--3770,
#' \doi{10.1093/bioinformatics/bti611}.
#'
#' Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2019). Wavelet-based
#' estimators for mixture regression. \emph{Scandinavian Journal of Statistics},
#' 46(1), 215--234, \doi{10.1111/sjos.12344}.
#'
#' Motta, F. C. and Montoril, M. H. (2026). A Bayesian estimation approach for
#' the wavelet-based mixture regression. \emph{Communications in Statistics -
#' Simulation and Computation}, 55(6), 2426--2434,
#' \doi{10.1080/03610918.2025.2470797}.
#'
#' @seealso \command{\link{wmixreg}}, \command{\link{bayesthresh}},
#'   \command{\link{hpdi}}, \command{\link{plot.bwmixreg}},
#'   \command{\link{predict.bwmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # Application of Section 5.2 of Motta and Montoril (2026).
#' #
#' # The data consist of array Comparative Genomic Hybridization values from a
#' # glioblastoma multiforme study (Lai et al., 2005). The peaks of the
#' # estimated curve indicate a higher probability of amplification of the DNA
#' # copy number at that chromosomal position.
#' #
#' if(requireNamespace("changepoint", quietly = TRUE)){
#'
#'   y <- changepoint::Lai2005fig4$GBM29
#'
#'   set.seed(123)
#'   fit <- bwmixreg(y, nchain = 1000, burn = 1000, lag = 5)
#'
#'   fit
#'   plot(fit)
#'
#'   # The three amplifications reported in the paper.
#'   which(fit$alpha > 0.5)
#' }
#'
#' #
#' # A simulated example with a known sinusoidal weight function, which is one
#' # of the scenarios of Section 5.1 of the paper.
#' #
#' set.seed(321)
#' n <- 512
#' t <- (1:n)/n
#' a <- 0.4*cos(2*pi*(t + pi)) + 0.5
#' z <- rbinom(n, size = 1, prob = a)
#' y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)
#'
#' fit <- bwmixreg(y, nchain = 500, burn = 500, plot = FALSE)
#'
#' plot(fit)
#' lines(t, a, col = 2, lwd = 2, lty = 2)
#' legend("topright", legend = c("Estimate", "Real function"), col = 1:2,
#'        lty = c(1, 2), lwd = 2, bty = "n")
#'
#' @keywords smooth
#' @importFrom stats quantile var median
#' @export
bwmixreg <- function(y, x = seq_along(y)/length(y),
                     nchain = 1000, burn = 1000, lag = 5,
                     prior = list(), init = list(), thresholding = list(),
                     family = "Coiflets", filter.size = 30,
                     wavelet.filter = NULL,
                     padding = c("reflect", "periodic", "none"),
                     level = 0.95, verbose = FALSE, plot = TRUE, ...){

  cl <- match.call()

  if(is.complex(y)){
    y <- Re(y)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

  y <- as.double(y)
  x <- as.double(x)
  n <- length(y)

  if(length(x) != n)
    stop("'x' and 'y' must have the same length.")
  if(any(!is.finite(y)) || any(!is.finite(x)))
    stop("'x' and 'y' must be finite numeric vectors.")
  if(n < 4)
    stop("The sample is too small.")
  if(is.unsorted(x))
    stop("'x' must be sorted, and the observations are assumed to be equally spaced.")

  nchain <- as.integer(nchain)
  burn <- as.integer(burn)
  lag <- as.integer(lag)

  if(is.na(nchain) || nchain < 2)
    stop("'nchain' must be at least two, so that the intervals can be computed.")
  if(is.na(burn) || burn < 0)
    stop("'burn' must be a non-negative integer.")
  if(is.na(lag) || lag < 1)
    stop("'lag' must be a positive integer.")
  if(!is.finite(level) || level <= 0 || level >= 1)
    stop("'level' must be a single value strictly between 0 and 1.")

  # --- The dyadic extension of the sample ---------------------------------
  padding <- match.arg(tolower(padding[1L]), c("reflect", "periodic", "none"))
  npad <- 2^ceiling(log2(n))

  if(padding == "none" && npad != n)
    stop("With padding = \"none\", the sample size must be a power of two.")

  # Zero-based positions of y, which C gathers to build the padded sample.
  # The reflection is the whole-point mirror of period 2n, which reproduces
  # the tail reversal of the reference implementation and keeps working when
  # more than n points have to be added.
  idx <- seq_len(npad) - 1L
  padidx <- switch(padding,
    none = idx,
    periodic = idx %% n,
    reflect = ifelse(idx %% (2L*n) < n, idx %% (2L*n),
                     2L*n - 1L - (idx %% (2L*n))))

  # --- Priors, initial values and thresholding settings -------------------
  q <- stats::quantile(y, probs = c(0.25, 0.75), names = FALSE)
  s2 <- stats::var(y)

  prior <- .wb_merge(prior, list(mean = q, var = c(s2, s2),
                                 shape = c(0.01, 0.01), rate = c(0.01, 0.01)),
                     "prior")
  init <- .wb_merge(init, list(mu = prior$mean, tau2 = c(1/s2, 1/s2)), "init")
  thresholding <- .wb_merge(thresholding,
                            list(j0 = 3, alpha = 0.5, beta = 1, C1 = NA,
                                 C2 = NA, C1.start = 100), "thresholding")

  for(nm in c("mean", "var", "shape", "rate"))
    if(length(prior[[nm]]) != 2L || any(!is.finite(prior[[nm]])))
      stop("The component '", nm, "' of 'prior' must hold two finite values.")
  if(any(prior$var <= 0) || any(prior$shape <= 0) || any(prior$rate <= 0))
    stop("The prior variances, shapes and rates must be positive.")

  for(nm in c("mu", "tau2"))
    if(length(init[[nm]]) != 2L || any(!is.finite(init[[nm]])))
      stop("The component '", nm, "' of 'init' must hold two finite values.")
  if(any(init$tau2 <= 0))
    stop("The initial precisions must be positive.")
  if(init$mu[1L] >= init$mu[2L])
    stop("The initial values must satisfy mu[1] < mu[2].")

  J <- as.integer(round(log2(npad)))
  j0 <- as.integer(thresholding$j0)
  if(is.na(j0) || j0 < 0 || j0 >= J)
    stop("The component 'j0' of 'thresholding' must be a non-negative integer smaller than ",
         J, ", the number of resolution levels of the padded sample.")
  if(!is.finite(thresholding$alpha) || thresholding$alpha < 0 ||
     !is.finite(thresholding$beta) || thresholding$beta < 0)
    stop("The components 'alpha' and 'beta' of 'thresholding' must be non-negative.")

  fam <- .wb_family_code(family, wavelet.filter)

  nverb <- if(is.logical(verbose)) (if(isTRUE(verbose)) 1000L else 0L)
           else as.integer(verbose)
  if(is.na(nverb) || nverb < 0)
    stop("'verbose' must be TRUE, FALSE or a non-negative integer.")

  # --- The Gibbs sampler ---------------------------------------------------
  draws <- .Call("_WaveBased_C_BayesMixReg", y, as.integer(padidx),
                 fam$fam, as.integer(filter.size), fam$filter, j0,
                 as.double(c(thresholding$alpha, thresholding$beta,
                             thresholding$C1, thresholding$C2,
                             thresholding$C1.start)),
                 as.double(c(prior$mean[1L], prior$var[1L], prior$shape[1L],
                             prior$rate[1L], prior$mean[2L], prior$var[2L],
                             prior$shape[2L], prior$rate[2L])),
                 as.double(c(init$mu, init$tau2)),
                 c(nchain, burn, lag, nverb))

  dimnames(draws$mu) <- list(NULL, c("mu1", "mu2"))
  dimnames(draws$tau2) <- list(NULL, c("tau2.1", "tau2.2"))

  # --- Posterior summaries -------------------------------------------------
  # The point estimates are posterior medians, which minimize the absolute
  # loss, as in Section 5 of the paper.
  pars <- cbind(draws$mu[, 1L], draws$tau2[, 1L], draws$mu[, 2L],
                draws$tau2[, 2L])
  hpd <- hpdi(pars, prob = level)

  estimates <- cbind(median = apply(pars, 2L, stats::median),
                     mean = colMeans(pars), lower = hpd[, 1L],
                     upper = hpd[, 2L])
  rownames(estimates) <- c("mu1", "tau2.1", "mu2", "tau2.2")

  alpha <- apply(draws$alpha, 2L, stats::median)

  out <- list(call = cl,
              alpha = alpha,
              alpha.mean = colMeans(draws$alpha),
              alpha.hpd = hpdi(draws$alpha, prob = level),
              estimates = estimates,
              mu = estimates[c("mu1", "mu2"), "median"],
              tau2 = estimates[c("tau2.1", "tau2.2"), "median"],
              z = draws$z,
              draws = draws[c("mu", "tau2", "alpha", "sigma", "C1", "C2")],
              x = x, y = y,
              nobs = n, npad = npad, padding = padding,
              nchain = nchain, burn = burn, lag = lag, level = level,
              prior = prior, init = init, thresholding = thresholding,
              family = fam$name, fam = fam$fam, filter.size = filter.size,
              wavelet.filter = if(fam$fam == 4L) wavelet.filter else NULL)

  names(out$mu) <- c("mu1", "mu2")
  names(out$tau2) <- c("tau2.1", "tau2.2")

  class(out) <- "bwmixreg"

  if(plot)
    plot(out, ...)

  out
}


#' @rdname bwmixreg
#' @param digits Significant digits in the printout.
#' @export
print.bwmixreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nBayesian wavelet-based mixture regression (bwmixreg)\n\n")
  cat("Call:", deparse(x$call), "\n\n", sep = " ")
  cat("  Wavelet basis  :",
      if(is.null(x$wavelet.filter))
        paste0(x$family, " (filter size ", x$filter.size, ")")
      else
        paste0("own filter (size ", length(x$wavelet.filter), ")"),
      "\n")
  cat("  Thresholding   : levels", x$thresholding$j0, "to",
      round(log2(x$npad)) - 1L, "| alpha =", x$thresholding$alpha,
      "| beta =", x$thresholding$beta, "\n")
  cat("  Data           : n =", x$nobs,
      if(x$npad != x$nobs)
        paste0("observations, padded to ", x$npad, " (", x$padding, ")")
      else "observations",
      "\n")
  cat("  MCMC           :", x$nchain, "draws kept, burn-in of", x$burn,
      "and lag of", x$lag, "\n\n")

  cat("  Component parameters (posterior median and ",
      format(100*x$level), "% HPD interval):\n", sep = "")

  labs <- c(mu1 = "mu[1]", tau2.1 = "tau^2[1]", mu2 = "mu[2]",
            tau2.2 = "tau^2[2]")

  for(p in rownames(x$estimates))
    cat("    ", format(labs[[p]], width = 10),
        format(signif(x$estimates[p, "median"], digits), width = 10),
        paste0("(", signif(x$estimates[p, "lower"], digits), "; ",
               signif(x$estimates[p, "upper"], digits), ")"), "\n", sep = "")

  cat("\n  Mixture weights: median in [",
      paste(signif(range(x$alpha), digits), collapse = ", "), "]\n\n", sep = "")

  invisible(x)
}


#' @title Plot method for the Bayesian wavelet-based mixture regression
#'
#' @description Plots the estimated dynamic mixture weights of a
#' \code{"bwmixreg"} object, together with their highest posterior density
#' band, as in Figure 3(b) of Motta and Montoril (2026).
#'
#' @param x A fitted object of class \code{"bwmixreg"}.
#' @param band Logical. Should the highest posterior density band be drawn?
#' @param data Logical. Should the observed mixture be added to the plot, on a
#'   secondary axis? It is useful to relate the peaks of the weights to the
#'   observations, as in Figure 3 of the paper.
#' @param estimate Which point estimate to draw, the posterior \code{"median"}
#'   (the default, as in the paper) or the posterior \code{"mean"}.
#' @param main,xlab,ylab,type,col,band.col Graphical parameters with sensible
#'   defaults.
#' @param ... Further graphical parameters passed to \command{\link{plot}}.
#'
#' @details
#' The function draws a single panel, so further curves may be superimposed on
#' it afterwards with \command{\link{lines}} or \command{\link{points}}, which
#' is how the true weight function is added in the examples of
#' \command{\link{bwmixreg}}.
#'
#' @return Nothing. The function is called for its side effect.
#'
#' @seealso \command{\link{bwmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(321)
#' n <- 256
#' t <- (1:n)/n
#' a <- 0.2*(t < 0.3) + 0.8*(t >= 0.3 & t < 0.7) + 0.3*(t >= 0.7)
#' z <- rbinom(n, size = 1, prob = a)
#' y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)
#'
#' fit <- bwmixreg(y, nchain = 300, burn = 300, plot = FALSE)
#'
#' plot(fit)
#' plot(fit, data = TRUE, band = FALSE)
#'
#' @keywords smooth
#' @method plot bwmixreg
#' @importFrom graphics par points polygon lines axis mtext
#' @export
plot.bwmixreg <- function(x, band = TRUE, data = FALSE,
                          estimate = c("median", "mean"), main,
                          xlab = "Index", ylab = "Estimated probability",
                          type = "l", col = 1, band.col = "grey80", ...){

  estimate <- match.arg(tolower(estimate[1L]), c("median", "mean"))
  fit <- if(estimate == "median") x$alpha else x$alpha.mean
  mn <- if(missing(main)) "" else main

  plot(x$x, fit, type = "n", main = mn, xlab = xlab, ylab = ylab,
       ylim = c(0, 1), ...)

  if(band)
    polygon(c(x$x, rev(x$x)),
            c(x$alpha.hpd[, "lower"], rev(x$alpha.hpd[, "upper"])),
            col = band.col, border = NA)

  if(data){
    # The observations live on their own scale, so they are rescaled to the
    # unit interval of the weights and given a secondary axis.
    rg <- range(x$y)
    points(x$x, (x$y - rg[1L])/diff(rg), col = "grey50", cex = 0.5)
    at <- pretty(rg)
    at <- at[at >= rg[1L] & at <= rg[2L]]
    axis(4, at = (at - rg[1L])/diff(rg), labels = at)
    mtext("Observed values", side = 4, line = 2.5, cex = par("cex.lab"))
  }

  lines(x$x, fit, type = type, col = col)

  invisible(NULL)
}


#' @title Predictions and estimates from a Bayesian mixture regression
#'
#' @description Evaluates the estimated mixture weights at arbitrary points, and
#' extracts the estimated component parameters.
#'
#' @param object A fitted object of class \code{"bwmixreg"}.
#' @param newx Vector of points at which the weights are evaluated. It defaults
#'   to the observation times used in the fit.
#' @param what Which quantity to return, the point estimate (\code{"estimate"},
#'   the default) or the bounds of the highest posterior density band
#'   (\code{"lower"}, \code{"upper"}), or all of them (\code{"all"}).
#' @param estimate Which point estimate to return, the posterior
#'   \code{"median"} (the default) or the posterior \code{"mean"}.
#' @param ... Not used.
#'
#' @details
#' Unlike \command{\link{predict.wmixreg}}, whose estimates are a scaling
#' function expansion that can be evaluated anywhere, here the mixture weights
#' are estimated on the grid of the observations, because the discrete wavelet
#' transform delivers them there. Values at other points are therefore obtained
#' by linear interpolation, and points outside the range of the data are not
#' extrapolated.
#'
#' The function \code{coef} returns the posterior medians of the component
#' parameters, and \code{fitted} returns the estimated mixture weights at the
#' observation times.
#'
#' @return A vector of estimated probabilities, or a matrix with the columns
#'   requested by \code{what}. For \code{coef}, a named vector with the
#'   estimates of \eqn{\mu_1}, \eqn{\tau_1^2}, \eqn{\mu_2} and \eqn{\tau_2^2}.
#'
#' @seealso \command{\link{bwmixreg}}, \command{\link{plot.bwmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(321)
#' n <- 256
#' t <- (1:n)/n
#' z <- rbinom(n, size = 1, prob = 0.8*t + 0.1)
#' y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)
#'
#' fit <- bwmixreg(y, nchain = 300, burn = 300, plot = FALSE)
#'
#' predict(fit, newx = c(0.1, 0.5, 0.9))
#' predict(fit, newx = c(0.1, 0.5, 0.9), what = "all")
#' coef(fit)
#'
#' @keywords smooth
#' @method predict bwmixreg
#' @importFrom stats predict approx
#' @export
predict.bwmixreg <- function(object, newx = NULL,
                             what = c("estimate", "lower", "upper", "all"),
                             estimate = c("median", "mean"), ...){

  what <- match.arg(tolower(what[1L]),
                    c("estimate", "lower", "upper", "all"))
  estimate <- match.arg(tolower(estimate[1L]), c("median", "mean"))

  fit <- if(estimate == "median") object$alpha else object$alpha.mean
  cols <- cbind(estimate = fit, lower = object$alpha.hpd[, "lower"],
                upper = object$alpha.hpd[, "upper"])

  keep <- if(what == "all") colnames(cols) else what
  cols <- cols[, keep, drop = FALSE]

  if(is.null(newx))
    return(if(what == "all") cols else as.vector(cols))

  newx <- as.double(newx)
  out <- vapply(keep,
                function(k) stats::approx(object$x, cols[, k], xout = newx,
                                          rule = 1)$y,
                numeric(length(newx)))
  out <- matrix(out, nrow = length(newx),
                dimnames = list(NULL, keep))

  # A single quantity is returned as a plain vector, whatever the number of
  # points requested.
  if(what == "all") out else as.vector(out)
}


#' @rdname predict.bwmixreg
#' @method coef bwmixreg
#' @importFrom stats coef
#' @export
coef.bwmixreg <- function(object, estimate = c("median", "mean"), ...){

  estimate <- match.arg(tolower(estimate[1L]), c("median", "mean"))
  object$estimates[, estimate]
}


#' @rdname predict.bwmixreg
#' @method fitted bwmixreg
#' @importFrom stats fitted
#' @export
fitted.bwmixreg <- function(object, estimate = c("median", "mean"), ...){

  estimate <- match.arg(tolower(estimate[1L]), c("median", "mean"))
  if(estimate == "median") object$alpha else object$alpha.mean
}


#==============================================================================
# Internal helpers
#==============================================================================

# Completes a user-supplied list of settings with the defaults, rejecting
# unknown names so that a misspelled entry is not silently ignored.
#' @keywords internal
#' @importFrom utils modifyList
#' @noRd
.wb_merge <- function(user, default, what){

  if(is.null(user))
    return(default)

  if(!is.list(user) || (length(user) > 0L && is.null(names(user))))
    stop("'", what, "' must be a named list.")

  unknown <- setdiff(names(user), names(default))
  if(length(unknown))
    stop("Unknown component", if(length(unknown) > 1L) "s" else "", " of '",
         what, "': ", paste(unknown, collapse = ", "), ". The components ",
         "available are ", paste(names(default), collapse = ", "), ".")

  utils::modifyList(default, user)
}
