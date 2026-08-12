#' @title Bayesian wavelet thresholding
#'
#' @description Regularizes a signal by the Bayesian wavelet thresholding of
#' Abramovich, Sapatinas and Silverman (1998), which places a spike and slab
#' prior on the detail wavelet coefficients and shrinks them towards the
#' posterior median. It is the regularization step of \command{\link{bwmixreg}}.
#'
#' @param x Vector with the signal to be regularized. Its length must be a
#'   power of two.
#' @param j0 The coarsest resolution level whose detail coefficients are
#'   thresholded. The coefficients of the levels below \code{j0} are left
#'   untouched, as is the usual practice. The default \code{j0 = 3} is the one
#'   of \code{wavethresh::threshold} and the one used by Motta and Montoril
#'   (2026).
#' @param alpha,beta The parameters \eqn{\alpha} and \eqn{\beta} of the prior.
#'   They must be non-negative. The defaults \eqn{\alpha = 0.5} and
#'   \eqn{\beta = 1} are the robust choice recommended by Abramovich, Sapatinas
#'   and Silverman (1998).
#' @param C1,C2 The hyperparameters \eqn{C_1} and \eqn{C_2} of the prior. When
#'   they are left as \code{NA}, the default, they are estimated from the data.
#'   See Details.
#' @param C1.start Sets the search interval \eqn{[0, 50\sqrt{C_{1,start}}]} of
#'   the estimation of \eqn{C_1}. It is only used when \code{C1} is \code{NA}.
#' @param family The family of wavelets to use, as in \command{\link{wavedec}}.
#' @param filter.size The size of the wavelet filter.
#' @param wavelet.filter Use this to provide your own filter, together with
#'   \code{family = "Own"}.
#'
#' @details
#' The signal is decomposed by the discrete wavelet transform down to the level
#' zero, and each detail coefficient \eqn{d_{jk}} of the levels
#' \eqn{j_0, \ldots, J - 1} is assigned the prior
#' \deqn{\pi_j N(0, v_j^2) + (1 - \pi_j)\delta_0(\theta_{jk}), \qquad
#'       v_j^2 = 2^{-\alpha j}C_1, \qquad \pi_j = \min(1, 2^{-\beta j}C_2),}
#' where \eqn{\delta_0} is a point mass at zero. The coefficients are then
#' replaced by the median of the resulting posterior distribution, which
#' thresholds them, and the signal is reconstructed.
#'
#' The noise level is estimated by the median absolute deviation of the detail
#' coefficients of the finest level. The hyperparameter \eqn{C_1} is estimated
#' by maximizing the marginal likelihood of the coefficients that survive a
#' universal hard threshold, and \eqn{C_2} follows from it in closed form.
#'
#' This function was implemented based on the \code{"BayesThresh"} policy of
#' \code{wavethresh::threshold.wd}, of G. P. Nason, which is the implementation
#' used by Motta and Montoril (2026). The estimation of \eqn{C_1} relies on a
#' one-dimensional minimizer implemented based on \code{Brent_fmin}, the engine
#' of \command{\link{optimize}}, itself a translation of the algorithm of Brent
#' (1973). Note that, since the wavelet filters of this package and those of
#' \pkg{wavethresh} follow different conventions, the coefficients of the two
#' implementations are indexed differently, and the regularized signals agree
#' only up to that convention.
#'
#' @return A list with components
#' \item{fit}{The regularized signal, of the same length as \code{x}.}
#' \item{coefficients}{The thresholded wavelet coefficients, decomposed down to
#'   the level zero.}
#' \item{sigma}{The estimated noise level.}
#' \item{C1, C2}{The hyperparameters of the prior, supplied or estimated. They
#'   are \code{NA} when the noise level is degenerate, in which case no
#'   thresholding is performed.}
#'
#' @references
#' Abramovich, F., Sapatinas, T. and Silverman, B. W. (1998). Wavelet
#' thresholding via a Bayesian approach. \emph{Journal of the Royal Statistical
#' Society Series B}, 60(4), 725--749,
#' \doi{10.1111/1467-9868.00151}.
#'
#' Brent, R. P. (1973). \emph{Algorithms for Minimization without Derivatives}.
#' Prentice-Hall, Englewood Cliffs.
#'
#' Motta, F. C. and Montoril, M. H. (2026). A Bayesian estimation approach for
#' the wavelet-based mixture regression. \emph{Communications in Statistics -
#' Simulation and Computation}, 55(6), 2426--2434,
#' \doi{10.1080/03610918.2025.2470797}.
#'
#' @seealso \command{\link{bwmixreg}}, \command{\link{wavedec}},
#'   \command{\link{waverec}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(123)
#' n <- 256
#' t <- (1:n)/n
#' f <- sin(4*pi*t)
#' y <- f + rnorm(n, sd = 0.3)
#'
#' bt <- bayesthresh(y, family = "Coiflets", filter.size = 30)
#'
#' plot(t, y, col = "grey60", cex = 0.6, xlab = "t", ylab = "y")
#' lines(t, f, col = 2, lwd = 2, lty = 2)
#' lines(t, bt$fit, lwd = 2)
#'
#' c(sigma = bt$sigma, C1 = bt$C1, C2 = bt$C2)
#'
#' @keywords smooth
#' @export
bayesthresh <- function(x, j0 = 3, alpha = 0.5, beta = 1, C1 = NA, C2 = NA,
                        C1.start = 100, family = "Coiflets", filter.size = 30,
                        wavelet.filter = NULL){

  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

  x <- as.double(x)
  n <- length(x)

  if(any(!is.finite(x)))
    stop("'x' must be a finite numeric vector.")

  J <- log2(n)
  if(J != round(J) || J < 1)
    stop("The length of 'x' must be a power of two, and at least two.")

  j0 <- as.integer(j0)
  if(is.na(j0) || j0 < 0 || j0 >= J)
    stop("'j0' must be a non-negative integer smaller than log2(length(x)).")

  if(!is.finite(alpha) || alpha < 0 || !is.finite(beta) || beta < 0)
    stop("'alpha' and 'beta' must be non-negative.")
  if(!is.na(C1) && (!is.finite(C1) || C1 < 0))
    stop("'C1' must be non-negative, or NA to be estimated.")
  if(!is.na(C2) && (!is.finite(C2) || C2 < 0))
    stop("'C2' must be non-negative, or NA to be estimated.")
  if(!is.finite(C1.start) || C1.start < 0)
    stop("'C1.start' must be non-negative.")

  fam <- .wb_family_code(family, wavelet.filter)

  out <- .Call("_WaveBased_C_BayesThresh", x, fam$fam, as.integer(filter.size),
               fam$filter, j0,
               as.double(c(alpha, beta, C1, C2, C1.start)))

  list(fit = out[[1L]], coefficients = out[[2L]],
       sigma = out[[3L]][["sigma"]], C1 = out[[3L]][["C1"]],
       C2 = out[[3L]][["C2"]])
}


#==============================================================================
# Internal helpers
#==============================================================================

# Resolves the 'family' and 'wavelet.filter' arguments into the integer code
# and the filter vector expected by the C entry points, with the same
# conventions and messages of the other exported functions of the package.
#' @keywords internal
#' @noRd
.wb_family_code <- function(family, wavelet.filter = NULL){

  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))

  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")

  if(is.null(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }

  list(fam = as.integer(fam), filter = as.double(wavelet.filter),
       name = c("Daublets", "Symmlets", "Coiflets", "Own")[fam])
}
