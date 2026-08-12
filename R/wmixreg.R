#' @title Wavelet-based estimators for mixture regression
#'
#' @description Estimates the time-varying mixing probability of a two-component
#' mixture observed over time, using the wavelet-based estimators of Montoril,
#' Pinheiro and Vidakovic (2019).
#'
#' @param y Vector with the observed mixture. At each time point, the process
#'   behaves as a random variable \eqn{U} with probability \eqn{f(t)}, or as a
#'   random variable \eqn{V} with probability \eqn{1 - f(t)}.
#' @param x Vector of observation times, of the same length as \code{y} and with
#'   values in \eqn{[0, 1]}. It defaults to the equally spaced grid
#'   \eqn{i/n, i = 1, \ldots, n}. The times are treated as order statistics, so
#'   \code{x} is sorted internally (and \code{y} is permuted accordingly).
#' @param mu.U,mu.V The means \eqn{\mu_U} and \eqn{\mu_V} of the two components.
#'   The theory assumes them to be known; in practice they may be replaced by
#'   consistent estimates. If they are not supplied, they are estimated from the
#'   sample means of the observations above and below \code{cut}.
#' @param cut Threshold used to estimate \code{mu.U} and \code{mu.V} when these
#'   are not supplied. Observations greater than \code{cut} are attributed to
#'   \eqn{U}, and the remaining ones to \eqn{V}.
#' @param J The resolution level of the multiresolution space onto which
#'   \eqn{f} is projected. It defaults to \eqn{\lceil 0.75 \log_2 n \rceil},
#'   the choice made in the application of the paper.
#' @param s The constant \eqn{s} defining the intervals
#'   \eqn{[T_{i-s/2}, T_{i+s/2}]}. It must be a positive even integer, and the
#'   value \code{s = 2} suggested in the paper is the default. \bold{Caution:}
#'   the accompanying MATLAB implementation uses a constant equal to
#'   \code{s/2}.
#' @param estimator Which raw estimator of the scaling coefficients to compute.
#'   One of \code{"both"} (the default), \code{"phi"} or \code{"integral"}. See
#'   Details.
#' @param regularization Character vector with the regularization methods to
#'   apply, any subset of \code{"hard"}, \code{"soft"}, \code{"efromovich"},
#'   \code{"loclin"} and \code{"none"}. See Details.
#' @param family The family of wavelets to use, as in \command{\link{PHI}}. The
#'   default follows the paper, which uses a Symmlet basis.
#' @param filter.size The size of the wavelet filter. The default is the 10-tap
#'   Daubechies least asymmetric filter used throughout the paper.
#' @param prec.wavelet The number of iterations of the Daubechies-Lagarias
#'   algorithm.
#' @param wavelet.filter Use this to provide your own filter, together with
#'   \code{family = "Own"}.
#' @param wavelet.table An optional object created by \command{\link{wtable}}.
#'   When provided, the scaling functions are evaluated by table lookup with
#'   linear interpolation instead of the Daubechies-Lagarias algorithm.
#' @param integral The method used to integrate the scaling functions in the
#'   estimator (9). \code{"primitive"} (the default) evaluates the integrals
#'   through the primitive of the scaling function, and \code{"trapezoid"}
#'   applies the trapezoidal rule of the original implementation. See Details.
#' @param delta,min.pts Length of the sub-intervals and minimum number of points
#'   of the trapezoidal rule. Only used when \code{integral = "trapezoid"}.
#' @param ngrid Number of sub-intervals per unit interval of the dyadic grid on
#'   which the primitive of the scaling function is built. It must be a power of
#'   two, and it is only used when \code{integral = "primitive"}.
#' @param j0 The coarsest resolution level of the decomposition used by the
#'   thresholding methods. Efromovich's shrinkage always uses \code{j0 = 0},
#'   because its coefficient reordering requires it.
#' @param clip Logical. Since \eqn{f} is a probability, should the estimates be
#'   truncated to the interval \eqn{[0, 1]}? Defaults to \code{TRUE}, as in the
#'   paper.
#' @param length.obs Number of equally spaced points at which the estimated
#'   curves are evaluated, when \code{newx} is not supplied.
#' @param newx Optional vector of points at which the estimated curves are
#'   evaluated.
#' @param bandwidth Optional bandwidth of the local linear regression. If it is
#'   not supplied, the rule of Bowman and Azzalini (1997) is used.
#' @param plot Logical. Should the estimated curves be plotted?
#' @param ... Further arguments passed to \command{\link{plot.wmixreg}}.
#'
#' @details
#' Let \eqn{Z_t} be a Bernoulli random variable with parameter \eqn{f(t)}, and
#' let \eqn{U_t} and \eqn{V_t} be independent of \eqn{Z_t}, with known means
#' \eqn{\mu_U} and \eqn{\mu_V}. The observed process is
#' \deqn{Y_t = Z_t U_t + (1 - Z_t) V_t,}
#' so that \eqn{Y_t} is either \eqn{U_t} or \eqn{V_t} with probabilities
#' \eqn{f(t)} and \eqn{1 - f(t)}. The interest lies in estimating the mixing
#' probability \eqn{f}.
#'
#' Defining \eqn{W_i = (Y_i - \mu_V)/(\mu_U - \mu_V)}, one has
#' \eqn{E(W_i \mid T_i = t) = f(t)}, which turns the mixture problem into the
#' heteroscedastic nonparametric regression \eqn{W_i = f(T_i) + \epsilon(T_i, f)}.
#' The function \eqn{f} is then approximated by its projection onto the
#' multiresolution space \eqn{V_J} generated by periodized scaling functions,
#' \deqn{f_J(t) = \sum_{k=1}^{2^J} c_{Jk}\phi_{Jk}(t),}
#' Two raw estimators of the coefficients are available through the argument
#' \code{estimator}:
#' \describe{
#' \item{\code{"phi"}}{The estimator (8) of the paper, which \emph{evaluates}
#'   the scaling functions at the observation times,
#'   \deqn{\hat c_{Jk} = \frac{1}{s}\sum_{i=1}^{n}(T_{i+s/2} - T_{i-s/2})\phi_{Jk}(T_i)W_i.}
#'   It reads as a rectangle rule of numerical integration.}
#' \item{\code{"integral"}}{The estimator (9) of the paper, which
#'   \emph{integrates} the scaling functions over the same intervals,
#'   \deqn{\tilde c_{Jk} = \frac{1}{s}\sum_{i=1}^{n}W_i\int_{T_{i-s/2}}^{T_{i+s/2}}\phi_{Jk}(t)\,\mathrm{d}t.}
#'   It is the refinement of the rectangle rule above that the approximation
#'   leading to (9) suggests.}
#' \item{\code{"both"}}{The default, which computes the two of them in a single
#'   pass over the data.}
#' }
#' In both cases the conventions are \eqn{T_{i-s/2} = T_{\max(0, i - s/2)}},
#' \eqn{T_{i+s/2} = T_{\min(i + s/2, n+1)}}, \eqn{T_0 = 0} and \eqn{T_{n+1} = 1}.
#' The two estimators share the same rates of convergence (Theorems 1 and 2 of
#' the paper), and the components of the fitted object are named after these
#' two labels.
#'
#' The integrals of the second estimator are obtained from the primitive
#' \eqn{\Phi(u) = \int_{-\infty}^{u}\phi} through
#' \deqn{\int_a^b \phi_{Jk}(t)\,\mathrm{d}t = 2^{-J/2}[\Phi(2^Jb - k) - \Phi(2^Ja - k)],}
#' which requires two evaluations of \eqn{\Phi} per non-null coefficient. This
#' is both faster and more accurate than the trapezoidal rule used in the
#' original implementation, which needs one evaluation of \eqn{\phi} per
#' sub-interval of the quadrature. The option \code{integral = "trapezoid"}
#' reproduces that rule and is kept mostly for comparison purposes. Neither
#' method builds the \eqn{n \times 2^J} matrix of basis functions.
#'
#' Four regularization methods are available:
#' \describe{
#' \item{\code{"hard"}, \code{"soft"}}{Thresholding in the wavelet domain. The
#'   raw coefficients are decomposed by the Mallat algorithm, the detail
#'   coefficients are thresholded with the universal policy
#'   \eqn{\lambda = \hat\sigma\sqrt{2J\log 2}}, and the result is reconstructed.
#'   The noise level \eqn{\hat\sigma} is estimated by the median absolute
#'   deviation of the detail coefficients from the finest level. Only
#'   \code{"hard"} is discussed in the paper.}
#' \item{\code{"efromovich"}}{The wavelet version of the shrinkage of Efromovich
#'   and Pinsker (1996). It selects automatically the resolution level that
#'   minimizes an empirical risk among the candidates \eqn{0, \ldots, J}, and
#'   shrinks the coefficients by \eqn{\hat\Theta_j/(\hat\Theta_j + \hat\sigma^2)}.
#'   The selected level is reported by \command{\link{print.wmixreg}}, and the
#'   whole risk sequence is kept in the fitted object. Note that Step 1 of
#'   Section 3.3.2 of the paper describes the raw coefficients as being computed
#'   at level \eqn{J + 1}; as in the implementation that accompanies the paper,
#'   they are computed here at level \eqn{J}, so that the candidate levels are
#'   \eqn{0, \ldots, J}. Whenever the empirical risk keeps decreasing, the
#'   finest candidate is selected and the method reduces to a pure shrinkage.
#'   This reordering of the coefficients into a single index requires the
#'   decomposition to start at level zero, so \code{j0} is ignored here.}
#' \item{\code{"loclin"}}{Regularization in the time domain. The raw estimate is
#'   evaluated at the observation times and smoothed by a local linear
#'   regression with a Gaussian kernel, using the bandwidth of Bowman and
#'   Azzalini (1997).}
#' \item{\code{"none"}}{No regularization, that is, the raw estimate itself.}
#' }
#'
#' @return An object of class \code{"wmixreg"}, that is, a list with components
#' \item{call}{The matched call.}
#' \item{coefficients}{List with the raw coefficient estimates of each
#'   estimator requested, each of length \eqn{2^J}.}
#' \item{fits}{Nested list indexed by estimator and by regularization method,
#'   holding the regularized coefficients and the quantities selected by each
#'   method (the noise level, the threshold, the level chosen by Efromovich's
#'   shrinkage and the bandwidth of the local linear regression).}
#' \item{fitted.values}{Matrix with the estimated curves evaluated at
#'   \code{newx}, with one column per combination of estimator and
#'   regularization method.}
#' \item{raw.fitted}{Matrix with the raw estimates evaluated at the observation
#'   times, with one column per estimator.}
#' \item{newx}{The points at which \code{fitted.values} was evaluated.}
#' \item{x, y, w}{The sorted observation times, the observed mixture and the
#'   transformed responses \eqn{W_i}.}
#' \item{mu.U, mu.V}{The means of the two components, supplied or estimated.}
#' \item{J, j0, s}{The resolution levels and the constant \eqn{s} used.}
#' \item{family, filter.size, prec.wavelet, wavelet.filter}{The wavelet basis
#'   used.}
#'
#' @references
#' Bowman, A. W. and Azzalini, A. (1997). \emph{Applied Smoothing Techniques for
#' Data Analysis}. Clarendon Press, Oxford.
#'
#' Donoho, D. L. and Johnstone, I. M. (1994). Ideal spatial adaptation by
#' wavelet shrinkage. \emph{Biometrika}, 81(3), 425--455,
#' \doi{10.1093/biomet/81.3.425}.
#'
#' Efromovich, S. and Pinsker, M. (1996). Sharp-optimal and adaptive estimation
#' for heteroscedastic nonparametric regression. \emph{Statistica Sinica}, 6(4),
#' 925--942.
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
#' @seealso \command{\link{plot.wmixreg}}, \command{\link{predict.wmixreg}},
#'   \command{\link{PHI}}, \command{\link{wtable}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # Application of Section 5 of Montoril, Pinheiro and Vidakovic (2019).
#' #
#' # The data consist of array Comparative Genomic Hybridization values from a
#' # glioblastoma multiforme study (Lai et al., 2005). Large values indicate
#' # chromosomal aberrations, and the interest lies in estimating the
#' # proportion of aberrations along the genome.
#' #
#' if(requireNamespace("changepoint", quietly = TRUE)){
#'
#'   y <- changepoint::Lai2005fig4$GBM29
#'   n <- length(y)
#'
#'   # The means of the two components are estimated by the sample means of the
#'   # observations above and below 2, which gives 4.594 and 0.249.
#'   fit <- wmixreg(y, cut = 2, J = ceiling(0.75*log2(n)),
#'                  family = "Symmlets", filter.size = 10)
#'
#'   fit
#'   plot(fit)
#'
#'   # The three amplifications reported in the paper are detected by the two
#'   # methods that regularize in the wavelet domain.
#'   round(coef(fit, estimator = "phi", regularization = "hard"), 4)
#' }
#'
#' #
#' # A simulated example with a known mixing probability.
#' #
#' set.seed(123)
#' n <- 500
#' t <- (1:n)/n
#' f <- 0.4*cos(2*pi*(t + pi)) + 0.5     # function f3 of the paper
#' z <- rbinom(n, size = 1, prob = f)
#' y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)
#'
#' fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 5, plot = FALSE)
#'
#' # A single estimator and a single regularization method draw one panel, so
#' # the true function can be superimposed on it.
#' plot(fit, estimator = "integral", regularization = "efromovich")
#' lines(t, f, col = 2, lwd = 2, lty = 2)
#' legend("topright", legend = c("Estimate", "Real function"), col = 1:2,
#'        lty = c(1, 2), lwd = c(1, 2), bty = "n")
#'
#' @keywords smooth
#' @importFrom stats mad
#' @importFrom graphics par lines
#' @export
wmixreg <- function(y, x = seq_along(y)/length(y), mu.U, mu.V, cut,
                    J = ceiling(0.75*log2(length(y))), s = 2,
                    estimator = c("both", "phi", "integral"),
                    regularization = c("hard", "efromovich", "loclin"),
                    family = "Symmlets", filter.size = 10, prec.wavelet = 30,
                    wavelet.filter, wavelet.table = NULL,
                    integral = c("primitive", "trapezoid"),
                    delta = 1e-4, min.pts = 2, ngrid = 2^13, j0 = 0,
                    clip = TRUE, length.obs = 250, newx = NULL,
                    bandwidth = NULL, plot = TRUE, ...){

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
  if(min(x) < 0 || max(x) > 1)
    stop("'x' must belong to the interval [0, 1].")

  # The theory is built upon the order statistics of the observation times.
  if(is.unsorted(x)){
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
  }

  estimator <- match.arg(tolower(estimator[1L]), c("both", "phi", "integral"))
  integral <- match.arg(tolower(integral[1L]), c("primitive", "trapezoid"))
  regularization <- unique(match.arg(tolower(regularization),
                                     c("hard", "efromovich", "loclin", "soft",
                                       "none"),
                                     several.ok = TRUE))

  if(s <= 0 || s != round(s) || s %% 2 != 0)
    stop("'s' must be a positive even integer.")

  J <- as.integer(J)
  j0 <- as.integer(j0)
  if(J < 1)
    stop("'J' must be a positive integer.")
  if(j0 < 0 || j0 >= J)
    stop("'j0' must be a non-negative integer smaller than 'J'.")

  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))

  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")

  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }

  # --- The means of the two components ------------------------------------
  if(missing(mu.U) || missing(mu.V)){
    if(missing(cut))
      stop("Provide 'mu.U' and 'mu.V', or a 'cut' from which they can be estimated.")
    if(!any(y > cut) || !any(y < cut))
      stop("'cut' does not split the data into two non-empty groups.")
    if(missing(mu.U))
      mu.U <- mean(y[y > cut])
    if(missing(mu.V))
      mu.V <- mean(y[y < cut])
  }

  if(!is.finite(mu.U) || !is.finite(mu.V) || mu.U == mu.V)
    stop("'mu.U' and 'mu.V' must be finite and different from each other.")

  # --- The mixture problem written as a regression problem ----------------
  w <- (y - mu.V)/(mu.U - mu.V)

  sn <- s/2
  if(sn >= n)
    stop("'s' is too large for the sample size.")
  lower <- c(rep(0, sn), x[seq_len(n - sn)])
  upper <- c(x[(sn + 1):n], rep(1, sn))

  # --- Wavelet tables -----------------------------------------------------
  # The table supplied by the user drives the pointwise evaluation of phi. A
  # table built here only to support the primitive must not silently replace
  # the exact evaluation, so the two are kept apart.
  evaltab <- NULL
  if(!is.null(wavelet.table))
    evaltab <- .match_wavelet_table(wavelet.table, fam, filter.size,
                                    wavelet.filter)[[1L]]

  # The primitive lives on a dyadic grid built by the cascade algorithm, so
  # 'ngrid' is passed to C as its base-two logarithm.
  mgrid <- log2(ngrid)
  if(!is.finite(mgrid) || mgrid != round(mgrid) || mgrid < 1 || mgrid > 24)
    stop("'ngrid' must be a power of two between 2 and 2^24.")

  what <- switch(estimator, phi = 1L, integral = 2L, both = 3L)

  raw <- .Call("_WaveBased_C_MixCoef", as.double(w), as.double(x),
               as.double(lower), as.double(upper), as.double(s),
               as.integer(J), as.integer(fam), as.integer(filter.size),
               as.integer(prec.wavelet), as.double(wavelet.filter),
               what, if(integral == "primitive") 0L else 1L,
               as.double(c(delta, min.pts)), evaltab, as.integer(mgrid))

  ests <- switch(estimator, phi = "phi", integral = "integral",
                 both = c("phi", "integral"))
  coefs <- stats::setNames(raw[c(if("phi" %in% ests) 1L,
                                 if("integral" %in% ests) 2L)], ests)

  # --- Regularization -----------------------------------------------------
  wargs <- if(fam == 4L)
             list(family = "Own", wavelet.filter = wavelet.filter)
           else
             list(family = family, filter.size = filter.size)

  dec <- function(v, lev) do.call(wavedec, c(list(x = v, j0 = lev), wargs))
  rec <- function(v, lev) do.call(waverec, c(list(x = v, j0 = lev), wargs))

  if(is.null(newx))
    newx <- seq(min(x), max(x), length.out = length.obs)
  newx <- as.double(newx)

  evalc <- function(coef, lev)
    .Call("_WaveBased_C_MixEval", as.double(coef), newx, as.integer(lev),
          as.integer(fam), as.integer(filter.size), as.integer(prec.wavelet),
          as.double(wavelet.filter), evaltab)

  evalx <- function(coef, lev)
    .Call("_WaveBased_C_MixEval", as.double(coef), as.double(x),
          as.integer(lev), as.integer(fam), as.integer(filter.size),
          as.integer(prec.wavelet), as.double(wavelet.filter), evaltab)

  raw.fitted <- vapply(ests, function(e) evalx(coefs[[e]], J), numeric(n))
  dimnames(raw.fitted) <- list(NULL, ests)

  fits <- stats::setNames(vector("list", length(ests)), ests)
  fitted.values <- matrix(NA_real_, length(newx),
                          length(ests)*length(regularization))
  cn <- character(0)
  col <- 0L

  for(e in ests){

    fits[[e]] <- stats::setNames(vector("list", length(regularization)),
                                 regularization)

    for(r in regularization){

      fi <- switch(r,
        none = list(coef = coefs[[e]], J = J),
        hard = ,
        soft = .wmix_thresh(coefs[[e]], j0 = j0, type = r, dec = dec,
                            rec = rec),
        efromovich = .wmix_efromovich(coefs[[e]], dec = dec, rec = rec),
        loclin = .wmix_loclin(x, raw.fitted[, e], bandwidth))

      fits[[e]][[r]] <- fi

      col <- col + 1L
      fitted.values[, col] <- if(r == "loclin")
                                .wmix_smooth(x, raw.fitted[, e], newx, fi$bw)
                              else
                                evalc(fi$coef, fi$J)
      cn <- c(cn, paste(e, r, sep = "."))
    }
  }

  if(clip)
    fitted.values[] <- pmin(pmax(fitted.values, 0), 1)

  dimnames(fitted.values) <- list(NULL, cn)

  out <- list(call = cl,
              coefficients = coefs,
              fits = fits,
              fitted.values = fitted.values,
              raw.fitted = raw.fitted,
              newx = newx,
              x = x, y = y, w = w,
              mu.U = mu.U, mu.V = mu.V,
              nobs = n, J = J, j0 = j0, s = s,
              estimator = ests, regularization = regularization,
              integral = integral, clip = clip,
              family = c("Daublets", "Symmlets", "Coiflets", "Own")[fam],
              fam = fam, filter.size = filter.size,
              prec.wavelet = prec.wavelet,
              wavelet.filter = if(fam == 4L) wavelet.filter else NULL,
              wavelet.table = wavelet.table)

  class(out) <- "wmixreg"

  if(plot)
    plot(out, ...)

  out
}


#' @rdname wmixreg
#' @param digits Significant digits in the printout.
#' @export
print.wmixreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("\nWavelet-based estimator for mixture regression (wmixreg)\n\n")
  cat("Call:", deparse(x$call), "\n\n", sep = " ")
  cat("  Wavelet basis  :",
      if(is.null(x$wavelet.filter))
        paste0(x$family, " (filter size ", x$filter.size, ")")
      else
        paste0("own filter (size ", length(x$wavelet.filter), ")"),
      "\n")
  cat("  Levels         : j0 =", x$j0, "and J =", x$J, "\n")
  cat("  Data           : n =", x$nobs, "observations, s =", x$s, "\n")
  cat("  Component means: mu.U =", signif(x$mu.U, digits),
      "and mu.V =", signif(x$mu.V, digits), "\n")
  cat("  Raw estimator  :",
      paste(c(phi = "(8) scaling functions",
              integral = "(9) integrated")[x$estimator], collapse = ", "),
      if("integral" %in% x$estimator) paste0("[", x$integral, "]") else "",
      "\n\n")

  cat("  Regularization:\n")
  for(e in x$estimator){
    for(r in x$regularization){
      fi <- x$fits[[e]][[r]]
      cat("    ", format(paste(e, r, sep = "/"), width = 20), sep = "")
      cat(switch(r,
                 none = "raw estimate",
                 hard = ,
                 soft = paste0("sigma = ", signif(fi$sigma, digits),
                               ", lambda = ", signif(fi$lambda, digits)),
                 efromovich = paste0("sigma = ", signif(fi$sigma, digits),
                                     ", selected level = ", fi$J),
                 loclin = paste0("bandwidth = ", signif(fi$bw, digits))),
          "\n", sep = "")
    }
  }
  cat("\n")

  invisible(x)
}


#' @title Plot method for wavelet-based mixture regression
#'
#' @description Plots the estimated mixing probability curves of a
#' \code{"wmixreg"} object, arranging the regularization methods along the rows
#' and the raw estimators along the columns, as in Figure 12 of Montoril,
#' Pinheiro and Vidakovic (2019).
#'
#' @param x A fitted object of class \code{"wmixreg"}.
#' @param estimator Optional subset of the raw estimators to plot.
#' @param regularization Optional subset of the regularization methods to plot.
#' @param data Logical. Should the transformed responses be added to the plot?
#'   It is only meaningful when a single panel is drawn.
#' @param main,xlab,ylab,type Graphical parameters with sensible defaults.
#' @param ... Further graphical parameters passed to \command{\link{plot}}.
#'
#' @details
#' When more than one panel is drawn, the function sets \code{par(mfrow)} and
#' restores it on exit. Superimposing curves afterwards with
#' \command{\link{lines}} or \command{\link{points}} therefore only works when a
#' single panel is requested, that is, when both \code{estimator} and
#' \code{regularization} select exactly one element. Otherwise the low-level
#' calls would be drawn on the restored layout instead of on the last panel.
#'
#' @return Nothing. The function is called for its side effect.
#'
#' @seealso \command{\link{wmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(123)
#' n <- 400
#' t <- (1:n)/n
#' f <- 0.8*t + 0.1
#' z <- rbinom(n, size = 1, prob = f)
#' y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)
#'
#' fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)
#' plot(fit)
#' plot(fit, estimator = "phi", regularization = "loclin")
#'
#' @keywords smooth
#' @method plot wmixreg
#' @importFrom graphics par points
#' @export
plot.wmixreg <- function(x, estimator = NULL, regularization = NULL,
                         data = FALSE, main, xlab = "Index",
                         ylab = "Estimated probability", type = "l", ...){

  ests <- if(is.null(estimator)) x$estimator
          else intersect(x$estimator, tolower(estimator))
  regs <- if(is.null(regularization)) x$regularization
          else intersect(x$regularization, tolower(regularization))

  if(length(ests) == 0 || length(regs) == 0)
    stop("No estimator/regularization combination available to plot.")

  labs <- c(none = "raw", hard = "hard threshold", soft = "soft threshold",
            efromovich = "Efromovich's shrinkage",
            loclin = "local linear regression")
  elab <- c(phi = "estimator (8)", integral = "estimator (9)")

  single <- (length(ests)*length(regs) == 1L)

  if(!single){
    op <- par(mfrow = c(length(regs), length(ests)))
    on.exit(par(op))
  }

  for(r in regs){
    for(e in ests){
      cname <- paste(e, r, sep = ".")
      mn <- if(missing(main)) paste0(elab[[e]], ", ", labs[[r]]) else main
      plot(x$newx, x$fitted.values[, cname], type = type, main = mn,
           xlab = xlab, ylab = ylab, ylim = c(0, 1), ...)
      if(data && single)
        points(x$x, x$w, col = "grey60", cex = 0.6)
    }
  }

  invisible(NULL)
}


#' @title Predictions and coefficients from a wavelet-based mixture regression
#'
#' @description Evaluates the estimated mixing probability at arbitrary time
#' points, and extracts the estimated coefficients.
#'
#' @param object A fitted object of class \code{"wmixreg"}.
#' @param newx Vector of time points at which the estimates are evaluated. It
#'   defaults to the grid used by \command{\link{wmixreg}}.
#' @param estimator Optional subset of the raw estimators. It defaults to all
#'   the estimators available in \code{object}.
#' @param regularization Optional subset of the regularization methods. It
#'   defaults to all the methods available in \code{object}.
#' @param clip Logical. Should the estimates be truncated to \eqn{[0, 1]}? It
#'   defaults to the choice made when the object was fitted.
#' @param ... Not used.
#'
#' @details
#' For the regularization methods that act in the wavelet domain, the estimates
#' are obtained by evaluating the scaling function expansion at \code{newx}. For
#' the local linear regression, which is a regularization in the time domain,
#' the raw estimate at the observation times is smoothed again at the requested
#' points, using the bandwidth selected when the object was fitted.
#'
#' The function \code{coef} returns the coefficients of a single combination of
#' estimator and regularization method. Note that the coefficients regularized
#' by Efromovich's shrinkage have length \eqn{2^{\hat J}}, where \eqn{\hat J} is
#' the level selected by the method, which is generally smaller than \eqn{J}.
#'
#' @return A matrix of estimated probabilities, with one row per element of
#'   \code{newx} and one column per combination of estimator and regularization
#'   method. For \code{coef}, a vector of coefficients.
#'
#' @seealso \command{\link{wmixreg}}, \command{\link{plot.wmixreg}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(123)
#' n <- 400
#' t <- (1:n)/n
#' f <- 0.8*t + 0.1
#' z <- rbinom(n, size = 1, prob = f)
#' y <- z*rnorm(n, mean = 2) + (1 - z)*rnorm(n, mean = 0)
#'
#' fit <- wmixreg(y, mu.U = 2, mu.V = 0, J = 4, plot = FALSE)
#' predict(fit, newx = c(0.1, 0.5, 0.9))
#' coef(fit, estimator = "phi", regularization = "hard")
#'
#' @keywords smooth
#' @method predict wmixreg
#' @importFrom stats predict
#' @export
predict.wmixreg <- function(object, newx = NULL, estimator = NULL,
                            regularization = NULL, clip = object$clip, ...){

  ests <- if(is.null(estimator)) object$estimator
          else intersect(object$estimator, tolower(estimator))
  regs <- if(is.null(regularization)) object$regularization
          else intersect(object$regularization, tolower(regularization))

  if(length(ests) == 0 || length(regs) == 0)
    stop("No estimator/regularization combination available.")

  if(is.null(newx))
    return(object$fitted.values[, paste(rep(ests, each = length(regs)), regs,
                                        sep = "."), drop = FALSE])

  newx <- as.double(newx)

  out <- matrix(NA_real_, length(newx), length(ests)*length(regs))
  cn <- character(0)
  col <- 0L

  for(e in ests){
    for(r in regs){
      fi <- object$fits[[e]][[r]]
      col <- col + 1L
      out[, col] <- if(r == "loclin")
                      .wmix_smooth(object$x, object$raw.fitted[, e], newx,
                                   fi$bw)
                    else
                      .Call("_WaveBased_C_MixEval", as.double(fi$coef), newx,
                            as.integer(fi$J), as.integer(object$fam),
                            as.integer(object$filter.size),
                            as.integer(object$prec.wavelet),
                            as.double(if(is.null(object$wavelet.filter)) 0
                                      else object$wavelet.filter),
                            if(is.null(object$wavelet.table)) NULL
                            else object$wavelet.table$phi)
      cn <- c(cn, paste(e, r, sep = "."))
    }
  }

  if(clip)
    out[] <- pmin(pmax(out, 0), 1)

  dimnames(out) <- list(NULL, cn)

  out
}


#' @rdname predict.wmixreg
#' @method coef wmixreg
#' @importFrom stats coef
#' @export
coef.wmixreg <- function(object, estimator = NULL, regularization = NULL, ...){

  e <- if(is.null(estimator)) object$estimator[1L] else tolower(estimator)[1L]
  r <- if(is.null(regularization)) object$regularization[1L]
       else tolower(regularization)[1L]

  if(!e %in% object$estimator)
    stop("Estimator '", e, "' is not available in this object.")
  if(!r %in% object$regularization)
    stop("Regularization '", r, "' is not available in this object.")

  if(r == "loclin")
    stop("The local linear regularization acts in the time domain and has no wavelet coefficients. Use predict() instead.")

  object$fits[[e]][[r]]$coef
}


#' @rdname predict.wmixreg
#' @method fitted wmixreg
#' @importFrom stats fitted
#' @export
fitted.wmixreg <- function(object, ...){
  object$fitted.values
}


#==============================================================================
# Internal helpers
#==============================================================================

#' @keywords internal
#' @noRd
.wmix_thresh <- function(coefs, j0, type, dec, rec){

  m <- length(coefs)
  J <- as.integer(round(log2(m)))

  wc <- dec(coefs, j0)
  sig <- mad(wc[(m/2 + 1):m])
  lambda <- sig*sqrt(2*J*log(2))

  idx <- -seq_len(2^j0)
  d <- wc[idx]
  wc[idx] <- if(type == "hard")
               d*(abs(d) > lambda)
             else
               sign(d)*pmax(abs(d) - lambda, 0)

  list(coef = rec(wc, j0), J = J, sigma = sig, lambda = lambda)
}


#' @keywords internal
#' @noRd
.wmix_efromovich <- function(coefs, dec, rec){

  m <- length(coefs)
  J <- as.integer(round(log2(m)))

  if(J < 1)
    stop("Efromovich's shrinkage needs at least two coefficients.")

  # The reordering of the coefficients into a single index requires the
  # decomposition to start at the coarsest possible level.
  wc <- dec(coefs, 0L)
  sig <- mad(wc[(m/2 + 1):m])

  theta <- pmax(wc^2 - sig^2, 0)
  emprisk <- seq_len(m)*sig^2 + (sum(theta) - cumsum(theta))

  risk <- numeric(J + 1L)
  risk[1:2] <- emprisk[1:2]
  if(J > 1)
    for(kk in 2:J)
      risk[kk + 1L] <- mean(emprisk[(2^(kk - 1) + 1):(2^kk)])

  Jest <- which.min(risk) - 1L

  ind <- seq_len(2^Jest)
  th <- theta[ind]
  shrunk <- if(sig == 0) wc[ind] else (th/(th + sig^2))*wc[ind]

  coef <- if(Jest == 0L) shrunk else rec(shrunk, 0L)

  list(coef = coef, J = Jest, sigma = sig, risk = risk)
}


#' @keywords internal
#' @noRd
.wmix_loclin <- function(x, fraw, bandwidth){

  n <- length(x)

  if(is.null(bandwidth)){
    hx <- mad(x)*(4/3/n)^0.2
    hy <- mad(fraw)*(4/3/n)^0.2
    bandwidth <- sqrt(hx*hy)
  }

  if(!is.finite(bandwidth) || bandwidth <= 0)
    stop("The bandwidth of the local linear regression is not positive. Provide 'bandwidth' explicitly.")

  list(coef = NULL, J = NA_integer_, bw = bandwidth)
}


#' @keywords internal
#' @noRd
.wmix_smooth <- function(x, fraw, xout, bw)
  .Call("_WaveBased_C_LocLin", as.double(x), as.double(fraw), as.double(xout),
        as.double(bw))
