#' @title Highest posterior density interval
#'
#' @description Computes the highest posterior density interval (HPDI) of MCMC
#' draws, that is, the shortest contiguous interval containing a proportion
#' \code{prob} of the draws. Unlike an equal-tailed (quantile) interval, the
#' HPDI is the shortest such interval, and it may be asymmetric for skewed
#' posteriors.
#'
#' @param data A numeric vector (a single chain of length \eqn{n}) or a numeric
#'   matrix in which each row is one MCMC draw and each column an independent
#'   chain, for instance one column per time point of a trajectory. This is the
#'   orientation of the \code{draws} component returned by
#'   \command{\link{bwmixreg}}. The values must be finite.
#' @param prob Probability mass contained in the interval, a single value
#'   strictly between 0 and 1. It defaults to \code{0.9}.
#'
#' @details
#' The computation is performed in C: each chain is sorted once and, over all
#' the windows of fixed span \eqn{m = \lfloor prob \cdot n \rfloor}, the window
#' of smallest width \eqn{x_{(i + m)} - x_{(i)}} is selected.
#'
#' This function was brought over from the \pkg{pdm} package of the same author,
#' where it supports the posterior summaries of the polynomial dynamic models.
#'
#' @return For a vector input, a named numeric vector of length two, with
#'   elements \code{lower} and \code{upper}. For a matrix input, a numeric
#'   matrix with one row per column of \code{data} and columns \code{lower} and
#'   \code{upper}.
#'
#' @seealso \command{\link{bwmixreg}}, \command{\link{bwregime}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' set.seed(1)
#'
#' # A single chain of a skewed posterior.
#' chain <- rgamma(2000, shape = 2, rate = 1)
#' hpdi(chain, prob = 0.9)
#'
#' # Several chains at once (rows are draws, columns are time points).
#' draws <- matrix(rnorm(2000*5), nrow = 2000, ncol = 5)
#' hpdi(draws, prob = 0.95)
#'
#' @keywords univar
#' @export
hpdi <- function(data, prob = 0.9){

  if(!is.numeric(data))
    stop("'data' must be a numeric vector or matrix.", call. = FALSE)

  if(length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1)
    stop("'prob' must be a single value strictly between 0 and 1.",
         call. = FALSE)

  # Ensure REALSXP without dropping the dim attribute of a matrix.
  storage.mode(data) <- "double"

  .Call("_WaveBased_C_hpdi", data, as.numeric(prob))
}
