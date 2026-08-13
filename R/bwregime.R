#' @title Identification of regime switches by Bayesian wavelet estimation
#'
#' @description Estimates jointly the dynamic mixture weights and the component
#' parameters of a two-component Gaussian mixture, by the Gibbs sampling
#' algorithm of Motta and Montoril (2026b). The weights are written as a probit
#' transformation of a wavelet expansion, which is estimated by data
#' augmentation under a spike and slab prior, so that regime switches are
#' identified without a functional form being imposed on the weights.
#'
#' @param y Vector with the observed mixture. At each index the process behaves
#'   as \eqn{N(\mu_2, \tau_2^{-2})} with probability \eqn{\alpha_t}, and as
#'   \eqn{N(\mu_1, \tau_1^{-2})} otherwise, with \eqn{\mu_1 < \mu_2}.
#' @param x Vector of observation times, of the same length as \code{y}, used
#'   only to place the estimates on an axis and to interpolate them in
#'   \command{\link{predict.bwregime}}. It defaults to the equally spaced grid
#'   \eqn{i/n, i = 1, \ldots, n}. The discrete wavelet transform requires
#'   equally spaced samples, so \code{x} is assumed to be sorted and equally
#'   spaced.
#' @param nchain Number of draws to retain.
#' @param burn Number of sweeps discarded as burn-in.
#' @param lag Every \code{lag}-th sweep after the burn-in is retained, which
#'   reduces the autocorrelation of the chains. The sampler therefore performs
#'   \code{burn + lag*nchain} sweeps. The default is the thinning used in the
#'   applications of the paper.
#' @param slab The slab component of the prior of the wavelet coefficients,
#'   either \code{"laplace"} (the default), the heavy tailed slab of Johnstone
#'   and Silverman (2005), or \code{"gaussian"}. They are the priors called SSL
#'   and SSG in the paper, which found the former to separate close probability
#'   peaks better in both applications.
#' @param link The link between the wavelet expansion and the mixture weights.
#'   Only the probit link of Albert and Chib (1993), used by the paper, is
#'   currently available.
#' @param prior A list with the hyperparameters of the priors of the component
#'   parameters, any of which may be omitted. The entries are \code{mean} and
#'   \code{var}, the means and the variances of the normal priors of
#'   \eqn{\mu_1} and \eqn{\mu_2}, and \code{shape} and \code{rate}, the
#'   parameters of the gamma priors of \eqn{\tau_1^2} and \eqn{\tau_2^2}, each
#'   of them of length two. The defaults are the choices of Section 4 of the
#'   paper, described in Details.
#' @param init A list with the initial values \code{mu} and \code{tau2} of the
#'   chains, each of length two. The defaults are the means of the priors of
#'   \eqn{\mu_k}, and the reciprocal of the sample variance for both
#'   \eqn{\tau_k^2}.
#' @param shrinkage A list with the settings of the spike and slab prior of the
#'   wavelet coefficients, any of which may be omitted. The entries are
#'   \code{zeta} and \code{rho}, the parameters of the beta prior of the
#'   sparsity parameters \eqn{\pi_j}, \code{kappa} and \code{xi}, the shape and
#'   the rate of the gamma prior of the slab parameter of each level, and
#'   \code{cut}, the smallest posterior probability of being non-null that does
#'   not exclude a coefficient outright. The default \code{cut = 0} is the
#'   method of the paper, and a positive value may improve the estimates
#'   considerably. See Details.
#' @param family The family of wavelets to use, as in \command{\link{wavedec}}.
#'   The default follows the applications of the paper, which use a Daubechies
#'   extremal phase basis.
#' @param filter.size The size of the wavelet filter. The default is the 20-tap
#'   extremal phase filter of the applications of the paper, which corresponds
#'   to \code{filter.number = 10} in the notation of \pkg{wavethresh} (the
#'   filters of that family have length twice the filter number).
#' @param wavelet.filter Use this to provide your own filter, together with
#'   \code{family = "Own"}.
#' @param padding How the sample is extended to a dyadic length, which the
#'   discrete wavelet transform requires. Either \code{"reflect"} (the default,
#'   a mirror reflection of the tail of the sample, as in the Remark of Section
#'   5.2 of the paper), \code{"periodic"} (a periodic extension) or
#'   \code{"none"} (which requires the sample size to be a power of two).
#' @param level The probability mass of the highest posterior density intervals
#'   reported.
#' @param verbose Logical, or a positive integer. Should the progress of the
#'   sampler be reported? An integer sets the reporting interval, in sweeps,
#'   whereas \code{TRUE} reports every 1000 sweeps.
#' @param plot Logical. Should the estimated mixture weights be plotted?
#' @param ... Further arguments passed to \command{\link{plot.bwregime}}.
#'
#' @details
#' The model of expression (6) of the paper is
#' \deqn{y_t = (1 - z_t)x_{1t} + z_t x_{2t}, \qquad
#'       x_{kt} \mid \mu_k, \tau_k^2 \sim N(\mu_k, \tau_k^{-2}), \qquad
#'       z_t \mid \alpha_t \sim \mathrm{Bern}(\alpha_t),}
#' where the \eqn{z_t} are allocation variables and the mixture weights
#' \eqn{\alpha_t} vary with the index \eqn{t}. Note that \eqn{\tau_k^2} is a
#' precision, so that the variance of the \eqn{k}-th component is
#' \eqn{\tau_k^{-2}}.
#'
#' The weights are modelled as the probit regression
#' \eqn{\alpha = \Phi(W^T\theta)}, where \eqn{W} is the matrix of the discrete
#' wavelet transform and \eqn{\theta} a vector of wavelet coefficients. Each
#' sweep of the Gibbs sampler (Algorithm 1 of the paper) draws the component
#' parameters from the conjugate full conditionals
#' \deqn{\mu_k \mid \cdots \sim N(b_k, B_k), \qquad
#'       \tau_k^2 \mid \cdots \sim \Gamma(c_k, C_k),}
#' permutes the labels of the pairs \eqn{(\mu_k, \tau_k^2)} whenever the draw
#' violates the identifying restriction \eqn{\mu_1 < \mu_2}, and draws the
#' allocation variables from their Bernoulli full conditionals. The latent
#' variables of the augmentation of Albert and Chib (1993) are then drawn from
#' \deqn{l_t \mid \cdots \sim N(w_t^T\theta, 1)}
#' truncated to the positive half line if \eqn{z_t = 1} and to the negative one
#' otherwise, which turns the problem into the nonparametric regression
#' \eqn{l = W^T\theta + e}. The coefficients of that regression are estimated by
#' wavelet shrinkage under the spike and slab prior
#' \deqn{\theta_{jk} \sim (1 - \pi_j)\delta_0(\theta_{jk}) + \pi_j\gamma(\theta_{jk}),}
#' where \eqn{\delta_0} is a point mass at zero and \eqn{\gamma} is the slab,
#' either a \eqn{N(0, v_j^2)} distribution or a Laplace distribution with scale
#' \eqn{a_j}. Unlike \command{\link{bayesthresh}}, whose hyperparameters are
#' estimated from the data, here the whole prior is sampled: the sparsity
#' parameters are given the prior \eqn{\pi_j \sim \mathrm{Beta}(\zeta, \rho)},
#' and the slab parameters the priors \eqn{v_j^{-2} \sim \Gamma(\kappa, \xi)}
#' or \eqn{a_j \sim \Gamma(\kappa, \xi)}. The scaling coefficient of the
#' coarsest level carries a diffuse prior and is never shrunk. This is the main
#' difference from \command{\link{bwmixreg}}, which estimates the weights by
#' regularizing the transformed responses \eqn{(y_t - \mu_1)/(\mu_2 - \mu_1)}
#' instead of a latent probit regression.
#'
#' The default priors of the component parameters are the ones used in Section 4
#' of the paper, namely \eqn{\mu_1 \sim N(q_1, s^2)} and
#' \eqn{\mu_2 \sim N(q_3, s^2)}, where \eqn{q_1} and \eqn{q_3} are the quartiles
#' and \eqn{s^2} the sample variance of the data, and
#' \eqn{\tau_k^2 \sim \Gamma(0.01, 0.01)}. The defaults of \code{shrinkage} are
#' the values of the reference implementation, \eqn{\zeta = \rho = \kappa = 1}
#' and \eqn{\xi = 100}; Section 4 of the paper describes the gamma priors as
#' \eqn{\Gamma(0.01, 0.01)}, which corresponds to
#' \code{shrinkage = list(kappa = 0.01, xi = 0.01)}. Following the paper, the
#' point estimates are posterior medians, which are optimal under the absolute
#' loss, and they are reported together with highest posterior density intervals
#' computed by \command{\link{hpdi}}.
#'
#' The entry \code{cut} of \code{shrinkage} excludes outright every coefficient
#' whose posterior probability of being non-null falls below it. The default
#' \code{cut = 0} keeps every coefficient the model selects, and is the method
#' as described in the paper and as implemented in the reference code of its
#' applications.
#'
#' A positive cut is nevertheless worth trying, and \code{cut = 0.05} is what
#' the reference code of the Monte Carlo experiments of the paper uses. The
#' reason it helps is that the marginal of a coefficient of pure noise is close
#' to the marginal of the spike whenever the slab is estimated to be as narrow
#' as the noise, which is what the default \eqn{\xi} produces at the finest
#' levels. The Bayes factor of expression (16) is then about one, the beta full
#' conditional (14) of \eqn{\pi_j} raises the sparsity parameter by about one
#' coefficient at every sweep, and the levels are progressively left unshrunk
#' until the weights interpolate the allocation variables. With \code{cut = 0}
#' the estimated weights are therefore closer to zero and one, as in Figures 7
#' and 9 of the paper, whereas a positive cut recovers the smooth weight
#' functions and the mean squared errors of Tables 2 to 5. Two other settings
#' address the same behaviour without excluding anything by hand, a wider slab
#' (say \code{xi = 1000}) and a sparsity favouring prior for \eqn{\pi_j} (say
#' \code{rho = 50}). The choice does not change the regimes identified in the
#' example below.
#'
#' The whole sampler runs in compiled code, including the wavelet transforms and
#' the shrinkage, so that no memory is allocated during the sweeps and only two
#' wavelet transforms are performed by sweep. It uses the random number
#' generator of R, and the chains are therefore reproducible through
#' \command{\link{set.seed}}.
#'
#' Two remarks on the assumptions. First, the discrete wavelet transform
#' requires equally spaced observations, and the sample size to be a power of
#' two; the latter is dealt with by \code{padding}, and the former is assumed.
#' Second, the estimates of \eqn{\alpha_t} live on the grid of the observations,
#' and not on a wavelet expansion, so \command{\link{predict.bwregime}}
#' interpolates them.
#'
#' This function was implemented based on the reference R code of the authors,
#' available at
#' \url{https://github.com/flaviamotta/Regime-switches-and-Bayesian-wavelet-estimation},
#' which uses the packages \pkg{wavethresh} and \pkg{EbayesThresh}. Some
#' deliberate differences are worth recording. The truncated normal draws and
#' the ratios of densities are evaluated on the logarithmic scale, which keeps
#' them accurate where the reference implementation guards divisions with a
#' small constant, that is, precisely on the plateaus where a regime is
#' identified with near certainty. The scaling coefficient of the coarsest level
#' is drawn from its full conditional under the diffuse prior, and not held at
#' the observed coefficient. The resolution levels are visited from the coarsest
#' to the finest one, the order in which the transform lays them out. The
#' retained sweeps are those of indices \code{burn + lag},
#' \code{burn + 2*lag}, and so on. And the wavelet filters of this package and
#' those of \pkg{wavethresh} follow different conventions, so that the chains
#' produced here are statistically equivalent to, but not numerically identical
#' with, the ones of the reference implementation.
#'
#' @return An object of class \code{"bwregime"}, that is, a list with components
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
#' \item{inclusion}{The posterior probability that each wavelet coefficient is
#'   non-null, in the order of the transform. Its first entry is the scaling
#'   coefficient, which is always included.}
#' \item{draws}{List with the retained draws of \code{mu} and \code{tau2}
#'   (\code{nchain} by 2 matrices), of \code{alpha} (a \code{nchain} by
#'   \eqn{n} matrix, whose columns are ready for \command{\link{hpdi}}), and of
#'   the level hyperparameters \code{pi} and \code{slab}, the latter holding
#'   \eqn{v_j^2} or \eqn{a_j} according to the slab used.}
#' \item{x, y}{The observation times and the observed mixture.}
#' \item{nobs, npad, padding}{The sample size, the padded sample size and the
#'   padding scheme.}
#' \item{nchain, burn, lag, level}{The settings of the sampler.}
#' \item{slab, link}{The prior of the wavelet coefficients and the link used.}
#' \item{prior, init, shrinkage}{The priors, the initial values and the
#'   hyperparameters actually used.}
#' \item{family, filter.size, wavelet.filter}{The wavelet basis used.}
#'
#' @references
#' Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. \emph{Journal of the American Statistical
#' Association}, 88(422), 669--679, \doi{10.1080/01621459.1993.10476321}.
#'
#' Johnstone, I. M. and Silverman, B. W. (2005). Empirical Bayes selection of
#' wavelet thresholds. \emph{The Annals of Statistics}, 33(4), 1700--1752,
#' \doi{10.1214/009053605000000345}.
#'
#' Lai, W. R., Johnson, M. D., Kucherlapati, R. and Park, P. J. (2005).
#' Comparative analysis of algorithms for identifying amplifications and
#' deletions in array CGH data. \emph{Bioinformatics}, 21(19), 3763--3770,
#' \doi{10.1093/bioinformatics/bti611}.
#'
#' Motta, F. C. and Montoril, M. H. (2026a). A Bayesian estimation approach for
#' the wavelet-based mixture regression. \emph{Communications in Statistics -
#' Simulation and Computation}, 55(6), 2426--2434,
#' \doi{10.1080/03610918.2025.2470797}.
#'
#' Motta, F. C. and Montoril, M. H. (2026b). Identifying regime switches through
#' Bayesian wavelet estimation: application to environmental and genetic data.
#' \emph{Journal of Applied Statistics}, \doi{10.1080/02664763.2025.2612551}.
#'
#' @seealso \command{\link{bwmixreg}}, \command{\link{wmixreg}},
#'   \command{\link{bayesthresh}}, \command{\link{hpdi}},
#'   \command{\link{plot.bwregime}}, \command{\link{predict.bwregime}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # Application of Section 5.2 of Motta and Montoril (2026b).
#' #
#' # The data consist of array Comparative Genomic Hybridization values from a
#' # glioblastoma multiforme study (Lai et al., 2005). The regions where the
#' # estimated weights are close to one are the ones where a DNA copy number
#' # alteration is present.
#' #
#' if(requireNamespace("changepoint", quietly = TRUE)){
#'
#'   y <- changepoint::Lai2005fig4$GBM29
#'
#'   # A lighter thinning than the one of the paper, to keep the example quick.
#'   set.seed(123)
#'   fit <- bwregime(y, nchain = 1000, burn = 1000, lag = 5)
#'
#'   fit
#'
#'   # The three regions of copy number alteration reported in the paper.
#'   which(fit$alpha > 0.5)
#'
#'   # The observations behind them.
#'   plot(fit, data = TRUE, band = FALSE)
#'
#'   \donttest{
#'   # The settings of the paper: 51,000 sweeps, of which 1,000 are discarded
#'   # and every 50th of the rest is kept.
#'   set.seed(123)
#'   paper <- bwregime(y, nchain = 1000, burn = 1000, lag = 50)
#'   coef(paper)
#'   }
#' }
#'
#' #
#' # A simulated example with a known sinusoidal weight function, which is one
#' # of the scenarios of Section 4 of the paper. A smooth weight function is
#' # where the cut of the inclusion probabilities discussed in the Details
#' # makes a difference.
#' #
#' set.seed(321)
#' n <- 512
#' t <- (1:n)/n
#' a <- 0.4*cos(2*pi*(t + pi)) + 0.5
#' z <- rbinom(n, size = 1, prob = a)
#' y <- z*rnorm(n, mean = 2, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)
#'
#' fit <- bwregime(y, nchain = 500, burn = 500, lag = 5, plot = FALSE)
#' cut <- bwregime(y, nchain = 500, burn = 500, lag = 5,
#'                 shrinkage = list(cut = 0.05), plot = FALSE)
#'
#' plot(cut, band = FALSE)
#' lines(t, fit$alpha, col = "grey60")
#' lines(t, a, col = 2, lwd = 2, lty = 2)
#' legend("topright", bty = "n", col = c(1, "grey60", 2), lty = c(1, 1, 2),
#'        lwd = c(1, 1, 2),
#'        legend = c("cut = 0.05", "cut = 0 (the default)", "Real function"))
#'
#' c(default = mean((fit$alpha - a)^2), cut = mean((cut$alpha - a)^2))
#'
#' @keywords smooth
#' @importFrom stats quantile var median
#' @export
bwregime <- function(y, x = seq_along(y)/length(y),
                     nchain = 1000, burn = 1000, lag = 50,
                     slab = c("laplace", "gaussian"), link = "probit",
                     prior = list(), init = list(), shrinkage = list(),
                     family = "Daublets", filter.size = 20,
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

  # --- The model ----------------------------------------------------------
  slab <- match.arg(tolower(slab[1L]), c("laplace", "gaussian"))
  link <- match.arg(tolower(link[1L]), "probit")

  model <- c(link = switch(link, probit = 1L),
             slab = switch(slab, gaussian = 1L, laplace = 2L),
             cprior = 1L)

  # --- The dyadic extension of the sample ---------------------------------
  padding <- match.arg(tolower(padding[1L]), c("reflect", "periodic", "none"))
  npad <- 2^ceiling(log2(n))

  if(padding == "none" && npad != n)
    stop("With padding = \"none\", the sample size must be a power of two.")

  # Zero-based positions of y, which C gathers to place the allocation
  # variables on the padded grid. The reflection is the whole-point mirror of
  # period 2n, which reproduces the tail reversal of the Remark of Section 5.2
  # and keeps working when more than n points have to be added.
  idx <- seq_len(npad) - 1L
  padidx <- switch(padding,
    none = idx,
    periodic = idx %% n,
    reflect = ifelse(idx %% (2L*n) < n, idx %% (2L*n),
                     2L*n - 1L - (idx %% (2L*n))))

  # --- Priors, initial values and shrinkage settings ----------------------
  q <- stats::quantile(y, probs = c(0.25, 0.75), names = FALSE)
  s2 <- stats::var(y)

  prior <- .wb_merge(prior, list(mean = q, var = c(s2, s2),
                                 shape = c(0.01, 0.01), rate = c(0.01, 0.01)),
                     "prior")
  init <- .wb_merge(init, list(mu = prior$mean, tau2 = c(1/s2, 1/s2)), "init")
  shrinkage <- .wb_merge(shrinkage,
                         list(zeta = 1, rho = 1, kappa = 1, xi = 100,
                              cut = 0),
                         "shrinkage")

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

  for(nm in c("zeta", "rho", "kappa", "xi"))
    if(length(shrinkage[[nm]]) != 1L || !is.finite(shrinkage[[nm]]) ||
       shrinkage[[nm]] <= 0)
      stop("The component '", nm, "' of 'shrinkage' must be a single positive value.")

  if(length(shrinkage$cut) != 1L || !is.finite(shrinkage$cut) ||
     shrinkage$cut < 0 || shrinkage$cut >= 1)
    stop("The component 'cut' of 'shrinkage' must be a single value in [0, 1).")

  fam <- .wb_family_code(family, wavelet.filter)

  nverb <- if(is.logical(verbose)) (if(isTRUE(verbose)) 1000L else 0L)
           else as.integer(verbose)
  if(is.na(nverb) || nverb < 0)
    stop("'verbose' must be TRUE, FALSE or a non-negative integer.")

  # --- The Gibbs sampler ---------------------------------------------------
  draws <- .Call("_WaveBased_C_BayesRegime", y, as.integer(padidx),
                 list(family = fam$fam,
                      filter.size = as.integer(filter.size),
                      filter = fam$filter),
                 model,
                 lapply(list(mean = prior$mean, var = prior$var,
                             shape = prior$shape, rate = prior$rate,
                             zeta = shrinkage$zeta, rho = shrinkage$rho,
                             kappa = shrinkage$kappa, xi = shrinkage$xi,
                             cut = shrinkage$cut),
                        as.double),
                 lapply(init[c("mu", "tau2")], as.double),
                 c(nchain, burn, lag, nverb))

  J <- as.integer(round(log2(npad)))
  dimnames(draws$mu) <- list(NULL, c("mu1", "mu2"))
  dimnames(draws$tau2) <- list(NULL, c("tau2.1", "tau2.2"))
  dimnames(draws$pi) <- list(NULL, paste0("j", seq_len(J) - 1L))
  dimnames(draws$slab) <- dimnames(draws$pi)

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
              inclusion = draws$inclusion,
              draws = draws[c("mu", "tau2", "alpha", "pi", "slab")],
              x = x, y = y,
              nobs = n, npad = npad, padding = padding,
              nchain = nchain, burn = burn, lag = lag, level = level,
              slab = slab, link = link,
              prior = prior, init = init, shrinkage = shrinkage,
              family = fam$name, fam = fam$fam, filter.size = filter.size,
              wavelet.filter = if(fam$fam == 4L) wavelet.filter else NULL)

  names(out$mu) <- c("mu1", "mu2")
  names(out$tau2) <- c("tau2.1", "tau2.2")

  class(out) <- "bwregime"

  if(plot)
    plot(out, ...)

  out
}


#' @rdname bwregime
#' @param digits Significant digits in the printout.
#' @export
print.bwregime <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nRegime switches through Bayesian wavelet estimation (bwregime)\n\n")
  cat("Call:", deparse(x$call), "\n\n", sep = " ")
  cat("  Wavelet basis  :",
      if(is.null(x$wavelet.filter))
        paste0(x$family, " (filter size ", x$filter.size, ")")
      else
        paste0("own filter (size ", length(x$wavelet.filter), ")"),
      "\n")
  cat("  Weights        :", x$link, "link with a", x$slab, "slab\n")
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
      paste(signif(range(x$alpha), digits), collapse = ", "), "]\n", sep = "")
  cat("  Sparsity       :",
      signif(sum(x$inclusion[-1L]), digits), "of", x$npad - 1L,
      "detail coefficients non-null, on average\n\n")

  invisible(x)
}


#' @title Plot method for the Bayesian regime switching model
#'
#' @description Plots the estimated dynamic mixture weights of a
#' \code{"bwregime"} object, together with their highest posterior density
#' band, as in Figures 7 and 9 of Motta and Montoril (2026b).
#'
#' @param x A fitted object of class \code{"bwregime"}.
#' @param band Logical. Should the highest posterior density band be drawn?
#' @param data Logical. Should the observed mixture be added to the plot, on a
#'   secondary axis? It is useful to relate the regimes to the observations, as
#'   in Figures 8 and 9 of the paper.
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
#' \command{\link{bwregime}}. It is the plot method of
#' \command{\link{bwmixreg}}, applied to the estimates of this model.
#'
#' @return Nothing. The function is called for its side effect.
#'
#' @seealso \command{\link{bwregime}}, \command{\link{plot.bwmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(321)
#' n <- 256
#' t <- (1:n)/n
#' a <- 0.05 + 0.9*(t > 0.3 & t < 0.5)
#' z <- rbinom(n, size = 1, prob = a)
#' y <- z*rnorm(n, mean = 3, sd = 0.5) + (1 - z)*rnorm(n, mean = 0, sd = 0.5)
#'
#' fit <- bwregime(y, nchain = 300, burn = 300, lag = 5, plot = FALSE)
#'
#' plot(fit)
#' plot(fit, data = TRUE, band = FALSE)
#'
#' @keywords smooth
#' @method plot bwregime
#' @export
plot.bwregime <- function(x, band = TRUE, data = FALSE,
                          estimate = c("median", "mean"), main = "",
                          xlab = "Index", ylab = "Estimated probability",
                          type = "l", col = 1, band.col = "grey80", ...){

  # The estimated weights of the two samplers are summarized in the same way,
  # so they are drawn by the same code.
  plot.bwmixreg(x, band = band, data = data, estimate = estimate, main = main,
                xlab = xlab, ylab = ylab, type = type, col = col,
                band.col = band.col, ...)
}


#' @title Predictions and estimates from a Bayesian regime switching model
#'
#' @description Evaluates the estimated mixture weights at arbitrary points, and
#' extracts the estimated component parameters.
#'
#' @param object A fitted object of class \code{"bwregime"}.
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
#' The mixture weights are estimated on the grid of the observations, because
#' that is where the discrete wavelet transform delivers them. Values at other
#' points are therefore obtained by linear interpolation, and points outside the
#' range of the data are not extrapolated.
#'
#' The function \code{coef} returns the posterior medians of the component
#' parameters, and \code{fitted} returns the estimated mixture weights at the
#' observation times. The three methods are the ones of
#' \command{\link{bwmixreg}}, applied to the estimates of this model.
#'
#' @return A vector of estimated probabilities, or a matrix with the columns
#'   requested by \code{what}. For \code{coef}, a named vector with the
#'   estimates of \eqn{\mu_1}, \eqn{\tau_1^2}, \eqn{\mu_2} and \eqn{\tau_2^2}.
#'
#' @seealso \command{\link{bwregime}}, \command{\link{plot.bwregime}},
#'   \command{\link{predict.bwmixreg}}
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
#' fit <- bwregime(y, nchain = 300, burn = 300, lag = 5, plot = FALSE)
#'
#' predict(fit, newx = c(0.1, 0.5, 0.9))
#' predict(fit, newx = c(0.1, 0.5, 0.9), what = "all")
#' coef(fit)
#'
#' @keywords smooth
#' @method predict bwregime
#' @importFrom stats predict
#' @export
predict.bwregime <- function(object, newx = NULL,
                             what = c("estimate", "lower", "upper", "all"),
                             estimate = c("median", "mean"), ...){

  predict.bwmixreg(object, newx = newx, what = what, estimate = estimate, ...)
}


#' @rdname predict.bwregime
#' @method coef bwregime
#' @importFrom stats coef
#' @export
coef.bwregime <- function(object, estimate = c("median", "mean"), ...){

  coef.bwmixreg(object, estimate = estimate, ...)
}


#' @rdname predict.bwregime
#' @method fitted bwregime
#' @importFrom stats fitted
#' @export
fitted.bwregime <- function(object, estimate = c("median", "mean"), ...){

  fitted.bwmixreg(object, estimate = estimate, ...)
}
