#' @title Scaling Function Basis
#'
#' @description Calculates a matrix of scaling functions \eqn{\phi_{Jk}}, where
#' the lines are related to the elements of the data set \eqn{x} and the columns
#' are related to the values of \eqn{k} for which \eqn{\phi_{Jk}} is non-null.
#'
#' @param x Vector containing the data. Its length does not need to be power of 2.
#' @param J The resolution level at which the scaling function basis must be
#'   calculated.
#' @param family The family of wavelets to use. It can be \emph{Daublets},
#'   \emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter
#'   is provided. See Details for more information.
#' @param filter.size The size of the filter. The available values depend on the
#'   chosen family. See Details for more information.
#' @param prec.wavelet The number of iterations to be performed in the
#'   Daubechies-Lagarias algorithm, which is used to evaluate the scaling
#'   functions of the specified wavelet basis at the data points.
#' @param wavelet.filter Use this to provide your own filter. To use this
#'   argument, you must specify \code{family = "Own"}. Do not use it, if you are
#'   not sure about what you are doing.
#' @param boundary The boundary treatment of the basis. One of
#'   \code{"periodic"} (the default, the periodized basis), \code{"none"}
#'   (the raw translates) or \code{"interval"} (the boundary-corrected
#'   orthonormal basis of Cohen, Daubechies and Vial, 1993, which requires
#'   the data to lie in [0, 1], a Daublet or Symmlet filter and a
#'   sufficiently large resolution level -- an informative error states the
#'   exact minimum). See \command{\link{wbasis}} for a detailed description
#'   of the interval basis.
#' @param wavelet.table An optional object created by \command{\link{wtable}}.
#'   When provided, the scaling functions are evaluated by table lookup with
#'   linear interpolation instead of the Daubechies-Lagarias algorithm, which
#'   is considerably faster for large data sets. The table must have been built
#'   for the same \code{family} and \code{filter.size} requested here; the
#'   argument \code{prec.wavelet} is then ignored. It is not used with
#'   \code{boundary = "interval"} (the boundary evaluation is always exact).
#'   See \command{\link{wtable}} for accuracy considerations.
#'
#' @details
#' The scaling function \eqn{\phi} is obtained according to a wavelet filter with
#' specific size, which may be provided by the specification of \code{family} and
#' \code{filter.size}, or by \code{wavelet.filter}. Such a function has an
#' associated multiresolution analysis
#' \deqn{\ldots \subset V_{J-1} \subset V_J \subset V_{J+1} \subset \ldots} and
#' generates the orthonormal basis \eqn{\{\phi_{Jk}(\cdot)\}_k} for the
#' multiresolution space \eqn{V_J}, where
#' \deqn{\phi_{Jk}(x) = 2^{J/2} \phi(2^J x - k).} This function calculates a
#' matrix of values \eqn{\phi_{Jk}(x)}, where the \eqn{i}-th line is related to
#' the \eqn{i}-th element of the data set \eqn{x} and the columns are related to
#' the values of \eqn{k} for which \eqn{\phi_{Jk}} is non-null.
#'
#' If \code{boundary = "none"}, the first column corresponds to the minimum
#' \eqn{k} where \eqn{\phi_{Jk}} is non-null for at least one element of \eqn{x}.
#' Analogously, the last column is related to the maximum reasonable value of
#' \eqn{k}. If \code{boundary = "periodic"}, the columns will correspond to
#' the values of \eqn{k} from \eqn{0} to \eqn{2^J - 1}.
#'
#' The scaling functions are evaluated at the data points efficiently, using the
#' Daubechies-Lagarias algorithm (Daubechies and Lagarias, 1992). Coded kindly
#' by Brani Vidakovic. Part of the code used here is based on the C-function
#' \code{phi}, in the \pkg{wavethresh} package.
#'
#' The families available here are \emph{Daublets}, \emph{Symmlets} and
#' \emph{Coiflets}. The first and the second correspond to the families of
#' Daubechies' extremal phase and Daubechies' least asymmetric wavelets,
#' respectively. The argument \code{family} here is not case sensitive, and it
#' accepts the first letter of the available families. For example, it is allowed
#' to use \code{family = "d"} for the \emph{Daublet} family. The filters used
#' here were taken from the Python package \pkg{PyWavelets} (Lee et al., 2019),
#' which provides extended precision for the families \emph{Daublets} and
#' \emph{Coiflets}. Furhtermore, \pkg{PyWavelets} provides a wider range of
#' filters for the three families.
#'
#' The filter sizes available depend on the chosen family:
#' \describe{
#' \item{Daublets}{\code{filter.size = 2} \code{(Haar), 4, 6, ..., 74, 76}.}
#' \item{Symmlets}{\code{filter.size = 4, 6, 8, ..., 38, 40}.}
#' \item{Coiflets}{\code{filter.size = 6, 12, 18, ..., 96, 102}.}
#' }
#'
#' The argument \code{filter.size} is associated with the number of vanishing
#' moments (\eqn{N}) of each wavelet family. For \emph{Daublets} and
#' \emph{Symmlets}, one has \eqn{\int x^n \psi(x) dx = 0, n = 0, \ldots, N - 1},
#' satisfying the relation \code{N = filter.size/2}. In the case of
#' \emph{Coiflets}, both scaling function \eqn{\phi} and wavelet \eqn{\psi} have
#' \eqn{N} vanishing moments satisfying the relation \code{N = filter.size/6}.
#'
#' This function accepects wavelet filters provided by the user. In this case,
#' the filter is given in the argument \code{wavelet filter}, while setting
#' \code{family = "Own"}. If the filter is provided and \code{family} is not
#' \code{"Other"}, a warning is released and the function will still consider the
#' wavelet filter given. \bold{Caution: the function does not check if the filter
#' provided is in fact a wavelet filter. Therefore, do not use
#' \code{wavelet.filter}, unless you know what you are doing.}
#'
#' This function works only with real valued observations. If, for some reason,
#' the data set is composed by complex numbers, the function will be applied only
#' to the real part of the data and a warning message will be released.
#' Furthermore, eventually, some observations of the data may not be finite
#' (\code{NA}, \code{NaN}, \code{Inf} or \code{-Inf}). If this happens, the
#' associated line to this observation in the matrix of the scaling function
#' basis will be composed by \code{NA} values, and a warning message will be
#' released.
#'
#' @return A matrix in which each line \eqn{i} corresponds to
#'   \eqn{\phi_{Jk}(x_i)}, for different values of \eqn{k} (associated to the
#'   columns of the matrix).
#'
#' @references
#' Daubechies, I. and Lagarias, J.C. (1992). Two-Scale Difference Equations II.
#' Local Regularity, Infinite Products of Matrices and Fractals. \emph{SIAM
#' Journal on Mathematical Analysis}, 24(4), 1031--1079,
#' \url{https://epubs.siam.org/doi/10.1137/0523059}.
#'
#' Lee, G., Gommers, R., Waselewski, F., Wohlfahrt, K. and O'Leary, A. (2019).
#' PyWavelets: A Python package for wavelet analysis. \emph{Journal of Open
#' Source Software}, 4(36), 1237, \url{https://doi.org/10.21105/joss.01237}.
#'
#' @seealso \code{\link{PSI}}, \code{\link{wbasis}}, \code{\link{wtable}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # Least squares estimation of the scaling function coefficients
#' # for the non-linear regression model
#' # Y = f(X) + e
#' #
#' set.seed(123)
#' n <- 100
#' x <- sort(runif(n))
#' f <- function(x) sin(2*pi*x)
#' e <- rnorm(n, mean = 0, sd = 0.25)
#' y <- f(x) + e
#'
#' # Let's apply the data to a periodized Daublet, built with a 18-tap filter.
#' # Consider J = 3 as the resolution level:
#' w <- PHI(x, J = 3, family = "daublets", filter.size = 18, prec.wavelet = 30,
#'          boundary = "periodic")
#'
#' # Estimating the scaling function coefficients
#' mod <- lm(y ~ w - 1)
#' alpha <- coef(mod) # estimates
#'
#' # Calculating estimates of f(x)
#' new.obs <- 0:(n-1)/n
#' myphi <- PHI(new.obs, J = 3, family = "daublets", filter.size = 18,
#'              prec.wavelet = 30, boundary = "periodic")
#' f.est <- drop(myphi %*% coef(mod)) # estimates of f(x)
#'
#' # Let's see the result
#' plot(x, y, main = "Regression model")
#' plot(f, 0, 1, lwd = 2, add = TRUE)
#' points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
#' legend("topright", legend = c("Real function", "Estimate"), col = 1:2, lty = 1)
#' \dontrun{#
#' #
#' # Next example
#' # ------------
#' #
#' # Comparison between wavethresh::denproj and PHI.
#' # Simulate data from the claw density and find the
#' # empirical scaling function coefficients
#' library(wavethresh)
#' data <- rclaw(100)
#' datahr <- denproj(data, J = 8, filter.number = 4, family = "DaubLeAsymm",
#'                   nT = 60)
#' coefdenproj <- round(datahr$coef, 8)
#'
#' # In our function, the analogous case is
#' matPHI <- PHI(data, J = 8, family = "symmlets", filter.size = 8,
#'               boundary = "none", prec.wavelet = 60)
#' coefPHI <- round(apply(matPHI, 2, mean), 8)
#'
#' identical(coefdenproj, coefPHI) # Are the coefficient estimates exact?
#' sum(abs(coefdenproj-coefPHI)) # Are they close with a good precision?
#' }
#'
#' @keywords smooth
#' @export
PHI <- function(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                wavelet.filter, boundary = c("periodic", "none", "interval"),
                wavelet.table = NULL){

  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

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

  bcode <- .wb_boundary_code(boundary)

  cdv <- NULL
  if(bcode == 2L)
    cdv <- .cdv_prepare(as.double(x), fam, filter.size, wavelet.filter,
                        level = J, what = "basis")

  wtab <- NULL
  if(!is.null(wavelet.table)){
    if(bcode == 2L)
      warning("wavelet.table is not used with boundary = \"interval\" and was ignored.")
    else
      wtab <- .match_wavelet_table(wavelet.table, fam, filter.size, wavelet.filter)
  }

  PHI <- .Call("_WaveBased_C_PHImat", as.double(x), as.integer(J), as.integer(fam),
               as.integer(filter.size), as.integer(prec.wavelet),
               bcode, as.double(wavelet.filter), wtab, cdv)
  return(PHI)

}
