#' @title Decomposed Wavelet Basis
#'
#' @description This function calculates, for each observation of the data, a
#' matrix of wavelet basis using the Daubechies-Lagarias algorithm. The first few
#' columns correspond to the scaling functions \eqn{\phi_{Jk}}, followed by the
#' wavelets at the coarsest level up to the wavelets at the finest resolution
#' level, i.e., \eqn{\psi_{jk}, j_0 \le j \le J-1}.
#'
#' @param x Vector containing the data. Its length does not need to be power of 2.
#' @param j0 The coarsest level at which the decomposed wavelet basis must be
#'   calculated. Default is \code{j0 = 0}.
#' @param J The resolution level at which the wavelet basis is decomposed to be
#'   calculated. It is associated to the multiresolution space of the calculated
#'   basis.
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
#'   orthonormal basis of Cohen, Daubechies and Vial, 1993). See Details for
#'   the requirements of \code{boundary = "interval"}.
#' @param wavelet.table An optional object created by \command{\link{wtable}}.
#'   When provided, the scaling functions and wavelets are evaluated by table
#'   lookup with linear interpolation instead of the Daubechies-Lagarias
#'   algorithm, which is considerably faster for large data sets. The same
#'   table serves any pair of levels \code{j0} and \code{J}. It must have been
#'   built for the same \code{family} and \code{filter.size} requested here;
#'   the argument \code{prec.wavelet} is then ignored. It is not used with
#'   \code{boundary = "interval"} (the boundary evaluation is always exact).
#'   See \command{\link{wtable}} for accuracy considerations.
#'
#' @details
#' The scaling function \eqn{\phi} and the wavelet \eqn{\psi} are obtained
#' according to a wavelet filter with specific size, which may be provided by the
#' specification of \code{family} and \code{filter.size}, or by
#' \code{wavelet.filter}. Such a function has an associated multiresolution
#' analysis
#' \deqn{\ldots \subset V_{J-1} \subset V_J \subset V_{J+1} \subset \ldots.}
#' The wavelet \eqn{\psi} generates an orthonormal basis
#' \eqn{\{\psi_{jk}(\cdot)\}_k} of the orthogonal complement \eqn{W_j} of the
#' multiresolution space \eqn{V_j} in \eqn{V_{j+1}}, i.e.,
#' \deqn{W_{j} = V_{j+1} \ominus V_{j},} for any integer \eqn{j}.
#' Therefore, one can see that
#' \deqn{V_J = V_{j_0} \oplus W_{j_0} \oplus \ldots \oplus W_{J-1},}
#' where \eqn{j_0} corresponds to the coarsest level and \eqn{J-1} is known as
#' the finest level of the wavelet decomposition. Thus, one have that
#' \deqn{\{\{\phi_{j_0k}(\cdot); k\} \cup \{\psi_{jk}(\cdot); j_0 \le j \le J-1, k\}\}}
#' compose an orthonormal basis for \eqn{V_J}. The basis above is provided by the
#' present function.
#'
#' Each line \eqn{i} of the resulting matrix is related to each element \eqn{x_i}
#' of \eqn{x}. The columns of the matrix are organized as
#' \deqn{\phi_{j_0,k_{\min}}(\cdot) \ldots \phi_{j_0,k_{\max}}(\cdot)
#' \psi_{j_0,k_{\min}}(\cdot) \ldots \psi_{J-1,k_{\max}}(\cdot).}
#' If \code{boundary = "periodic"}, then \eqn{k_{\min} = 0} and
#' \eqn{k_{\max} = 2^(j-1)-1}, for each level \eqn{j}. If
#' \code{boundary = "none"}, \eqn{k_{\min}}
#' and \eqn{k_{\max}} will be all the integers where their corresponding
#' function, \eqn{\phi_{j_0k}} or \eqn{\psi_{jk}, j_0 \le j \le J - 1}, is
#' non-null for at least one observation of \eqn{x}.
#'
#' With \code{boundary = "interval"}, the basis used is the orthonormal
#' wavelet basis of \eqn{L^2[0,1]} of Cohen, Daubechies and Vial (1993),
#' often called \emph{wavelets on the interval}. At each level \eqn{j}, the
#' basis keeps the \eqn{2^j - L} interior translates (\eqn{L} denotes the
#' filter size) that are entirely supported inside \eqn{[0,1]} and completes
#' them with \eqn{L/2} boundary functions at each endpoint, constructed so
#' that polynomials up to degree \eqn{L/2 - 1} are exactly reproduced on the
#' whole interval -- in particular, no artificial wrap-around discontinuity is
#' introduced at the boundaries, in contrast with the periodized basis. The
#' columns of each level block are organized as left boundary functions,
#' interior translates and right boundary functions, and the resulting matrix
#' has \eqn{2^J} columns, exactly as in the periodic case. The boundary
#' functions agree with the ones tabulated by Cohen, Daubechies and Vial up
#' to an orthogonal rotation within each boundary block, which spans the
#' same multiresolution spaces and therefore yields the same fitted
#' projections.
#'
#' For \emph{Daublets} and \emph{Symmlets} with \code{filter.size} between 4
#' and 24, the boundary coefficient blocks are precomputed: they were derived
#' from the filters in 320-bit arithmetic (see
#' \file{tools/make_cdv_tables.R} in the package sources) and are accurate to
#' the precision of the filter coefficients themselves. For other filters
#' (\code{family = "Own"}, or non-tabulated sizes) the construction is
#' carried out at run time in double precision, which is reliable up to
#' moderate filter sizes and validates itself, raising an informative error
#' when the filter is too long for a stable construction.
#'
#' The interval basis requires: (i) the data to lie in \eqn{[0,1]};
#' (ii) \code{family} to be \emph{Daublets} or \emph{Symmlets} (the
#' construction needs minimal-length filters, so \emph{Coiflets} are not
#' supported); and (iii) the coarsest level to be large enough for the two
#' boundary zones not to interact -- the function reports the exact minimum
#' (roughly \eqn{j_0 \ge \log_2(5 L)}) when the requested level is too small.
#'
#' The scaling functions \eqn{\phi_{j_0k}} and wavelets \eqn{\psi_{jk}} are
#' internally evaluated at the data points efficiently, using the
#' Daubechies-Lagarias algorithm (Daubechies and Lagarias, 1992).
#'
#' With respect to the levels \eqn{j_0} and \eqn{J}, it is important to mention
#' that the function only works for integer values satisfying
#' \eqn{0 \le j_0 \le J}. In the case where \eqn{j_0 = J}, the function will
#' assume the decomposition of the multiresolution space as \eqn{V_J = V_{j_0}}.
#' Therefore, it will take into account an orthonormal basis only with scaling
#' functions \eqn{\{\phi_{j_0k}(\cdot)\}_k}.
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
#'   \deqn{\{\{\phi_{j_0k}(x_i); k\} \cup \{\psi_{jk}(x_i); j_0 \le j \le J-1, k\}\},}
#'   for different values of \eqn{k}.
#'
#' @references
#' Cohen, A., Daubechies, I. and Vial, P. (1993). Wavelets on the Interval and
#' Fast Wavelet Transforms. \emph{Applied and Computational Harmonic
#' Analysis}, 1(1), 54--81,
#' \url{https://doi.org/10.1006/acha.1993.1005}.
#'
#' Daubechies, I. and Lagarias, J.C. (1992). Two-Scale Difference Equations II.
#' Local Regularity, Infinite Products of Matrices and Fractals. \emph{SIAM
#' Journal on Mathematical Analysis}, 24(4), 1031--1079,
#' \url{https://epubs.siam.org/doi/10.1137/0523059}.
#'
#' Lee, G., Gommers, R., Waselewski, F., Wohlfahrt, K. and O'Leary, A. (2019).
#' PyWavelets: A Python package for wavelet analysis. \emph{Journal of Open
#' Source Software}, 4(36), 1237, \url{https://doi.org/10.21105/joss.01237}.
#'
#' @seealso \command{\link{PHI}}, \command{\link{PSI}}, \command{\link{wtable}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # Wavelet estimation of a non-linear regression model
#' # Y = f(X) + e
#' #
#' set.seed(123)
#' n <- 100
#' x <- sort(runif(n))
#' f <- function(x) sin(2*pi*x)
#' e <- rnorm(n, mean = 0, sd = 0.25)
#' y <- f(x) + e
#'
#' w <- wbasis(x, j0 = 0, J = 3, family = "Daublets", filter.size = 18,
#'             prec.wavelet = 30)
#'
#' # Estimating the wavelet coefficients
#' mod <- lm(y~w-1)
#' beta <- coef(mod) # estimates
#'
#' # Calculating estimates of f(x)
#' new.obs <- 0:(n-1)/n
#' myw <- wbasis(new.obs, j0 = 0, J = 3, family = "Daublets", filter.size = 18,
#'              prec.wavelet = 30)
#' f.est <- drop(myw %*% beta) # estimates of f(new.obs)
#'
#' # Let's see the result
#' plot(x, y, main = "Regression model")
#' plot(f, 0, 1, lwd = 2, add = TRUE)
#' points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
#' legend("topright", legend = c("Real function", "Estimate"), col = 1:2, lty = 1)
#' #
#' #
#' # Next example
#' # ------------
#' #
#' # As above, but now with a different value of J and
#' # regularizing the wavelet coefficients by the hard threshold
#' # with lambda = 0.1
#' # Y = f(X) + e
#' #
#' set.seed(123)
#' n <- 100
#' x <- sort(runif(n))
#' f <- function(x) sin(2*pi*x)
#' e <- rnorm(n, mean = 0, sd = 0.25)
#' y <- f(x) + e
#'
#' w <- wbasis(x, j0 = 0, J = 5, family = "Daublets", filter.size = 18,
#'             prec.wavelet = 30)
#'
#' # Estimating the wavelet coefficients
#' mod <- lm(y~w-1)
#' beta <- coef(mod) # estimates
#'
#' # Calculating estimates of f(x)
#' new.obs <- 0:(n-1)/n
#' myw <- wbasis(new.obs, j0 = 0, J = 5, family = "Daublets", filter.size = 18,
#'              prec.wavelet = 30)
#' f.est <- drop(myw %*% beta) # estimates of f(new.obs)
#'
#' # Regularizing beta
#' lambda <- 0.1 # threshold
#' beta.thr <- beta*(abs(beta)>lambda)
#' f.est.thr <- drop(myw %*% beta.thr) # regularized estimates of f(new.obs)
#'
#' # Let's see the result
#' plot(x, y, main = "Regression model")
#' plot(f, 0, 1, lwd = 2, add = TRUE)
#' points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
#' points(new.obs, f.est.thr, col = 3, lwd = 2, type = 'l')
#' legend("topright", legend = c("Real function", "Raw Estimate", "Reg. Estimate"),
#'        col = 1:3, lty = 1)
#' #
#' #
#' # Next example
#' # ------------
#' #
#' # Wavelets on the interval (Cohen-Daubechies-Vial): same regression
#' # problem, but with a function that does not match at the boundaries
#' # (f(0) != f(1)), where the periodized basis suffers from wrap-around
#' # artifacts and the interval basis does not.
#' #
#' set.seed(123)
#' n <- 500
#' x <- sort(runif(n))
#' f <- function(x) x + sin(2*pi*x)
#' y <- f(x) + rnorm(n, sd = 0.2)
#'
#' w <- wbasis(x, j0 = 4, J = 5, family = "Daublets", filter.size = 8,
#'             boundary = "interval")
#' beta <- coef(lm(y~w-1))
#'
#' new.obs <- 0:200/200
#' myw <- wbasis(new.obs, j0 = 4, J = 5, family = "Daublets", filter.size = 8,
#'               boundary = "interval")
#' f.est <- drop(myw %*% beta)
#'
#' plot(x, y, main = "Regression with the interval basis")
#' plot(f, 0, 1, lwd = 2, add = TRUE)
#' points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
#' legend("topleft", legend = c("Real function", "Estimate"), col = 1:2, lty = 1)
#'
#' @keywords smooth
#' @export
wbasis <- function(x, j0 = 0, J, family = "Daublets", filter.size = 20,
                   prec.wavelet = 30, wavelet.filter,
                   boundary = c("periodic", "none", "interval"),
                   wavelet.table = NULL){

  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

  if(j0 < 0)
    stop("j0 must be a non-negative integer.")

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
                        level = j0,
                        what = if (j0 < J) "decompose" else "basis")

  wtab <- NULL
  if(!is.null(wavelet.table)){
    if(bcode == 2L)
      warning("wavelet.table is not used with boundary = \"interval\" and was ignored.")
    else
      wtab <- .match_wavelet_table(wavelet.table, fam, filter.size, wavelet.filter)
  }

  wmat <- .Call("_WaveBased_C_WavBasis", as.double(x), as.integer(j0), as.integer(J),
                as.integer(fam), as.integer(filter.size),
                as.integer(prec.wavelet), bcode,
                as.double(wavelet.filter), wtab, cdv)

  return(wmat)

}
