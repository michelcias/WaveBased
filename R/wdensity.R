#' @title Wavelet-based density estimation for univariate size-biased data
#'
#' @description The \code{wdensity} function calculates wavelet-based density
#' estimates for univariate size-biased data. The method is based on the
#' estimation of the power of the density of interest. It uses the
#' Daubechies-Lagarias algorithm to calculate the periodized version of a
#' compactly supported orthonormal wavelet basis.
#'
#' @param data Vector containing the data from which the density estimate is to
#'   be computed. The length of this vector does not have to be a power of 2.
#' @param from,to Minimum and maximum values of the grid of points at which the
#'   density is to be estimated. The default is \code{from = min(data)} and
#'   \code{to = max(data)}.
#' @param length.obs Number of equally spaced points at which the density is to
#'   be estimated. The default is \code{length.obs = 250}.
#' @param wf Weighting function. It must be a function.
#' @param power.dens Power of the density to be estimated. It must be a value
#'   greater than or equal to 0.5. The Default is \code{power.dens = 0.5}.
#' @param J0 The coarsest resolution level for the wavelet basis. Default is
#'   \code{J0 = 0}.
#' @param J1 The resolution level at which the wavelet basis is decomposed to be
#'   calculated. It is associated to the multiresolution space where the function
#'   is projected.
#' @param family The family of wavelets to use. It can be \emph{Daublets},
#'   \emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter
#'   is provided \code{wavelet.filter}.
#' @param filter.size The size of the wavelet filter. The available values depend
#'   on the chosen wavelet family.
#' @param prec.wavelet The number of iterations to be performed in the
#'   Daubechies-Lagarias algorithm.
#' @param wavelet.filter Use this to provide your own filter, the low-pass
#'   filter of an orthonormal wavelet basis. To use this argument, you must
#'   specify \code{family = "Own"}. The filter is checked: see Details of
#'   \command{\link{wbasis}}.
#' @param check.lawton Logical. Whether a filter given in
#'   \code{wavelet.filter} must satisfy the condition of Lawton, as in
#'   \command{\link{PHI}}.
#' @param wavelet.table An optional object created by \command{\link{wtable}}.
#'   When provided, the wavelet basis is evaluated by table lookup with linear
#'   interpolation instead of the Daubechies-Lagarias algorithm, which is
#'   considerably faster for large data sets. The table must have been built
#'   for the same \code{family} and \code{filter.size} requested here. See
#'   \command{\link{wtable}} for accuracy considerations.
#' @param dens.biased Density estimates of the size-biased data. These estimates
#'   of the (ordered) \code{data} must be provided. By default, a kernel density
#'   estimate is calculated.
#' @param warped Logical. It indicates if the wavelet basis must be warped.
#' @param warp.fun A (cumulative distribution) function. By default, it provides
#'   the empirical distribution function of \code{data}. Only used if
#'   \code{warped = TRUE}.
#' @param dwarp.fun A function. Derivative of the \code{warp.fun}. By default, it
#'   provides a kernel density estimate of \code{data}. Only used if
#'   \code{warped = TRUE}.
#' @param rescale Logical. If TRUE, the \code{data} will be rescaled in the
#'   interval [\code{eps}, 1-\code{eps}], where \code{eps} is another argument.
#'   Only used if \code{warped = FALSE}.
#' @param eps A value in \eqn{[0, 0.5)}. It is used to rescale the \code{data}
#'   in the interval \eqn{[}\code{eps}\eqn{, 1-}\code{eps}\eqn{]}. By default,
#'   \code{eps} \eqn{= 1.9^{-J1}}. A value outside \eqn{[0, 0.5)}, supplied or
#'   given by the default, stops with an error, since the interval would then
#'   be a single point or reversed. Only used if \code{warped = FALSE} and
#'   \code{rescale = TRUE}.
#' @param thresh Indicates if the wavelet-density estimate should be regularized
#'   and how. If \code{thresh = "linear"}, the density estimate will not be
#'   regularized; if \code{thresh = "hard"}, the hard-thresholding regularizer is
#'   applied to the detail coefficient estimates; and if \code{thresh = "soft"},
#'   the soft-thresholding regularizer is applied to the detail coefficient
#'   estimates.
#' @param boundary The boundary treatment of the wavelet basis and of the
#'   thresholding transform. Either \code{"periodic"} (default) or
#'   \code{"interval"} (the boundary-corrected basis of Cohen, Daubechies
#'   and Vial, 1993; see \command{\link{wbasis}}). With
#'   \code{"interval"}, the levels \code{J0} and \code{J1} must satisfy
#'   the minimum reported by the corresponding error message.
#' @param correct Logical. If \code{TRUE}, the density estimate is corrected to
#'   be a bona fide density by the method of Glad, Hjort and Ushakov (2003): it
#'   becomes \eqn{\max(0, \hat{f} - \xi)}, where the constant \eqn{\xi} makes it
#'   integrate to one over the domain of the estimator. See Details. Default is
#'   \code{correct = FALSE}.
#' @param plot Logical. If \code{plot = TRUE}, the density estimate is to be
#'   plotted.
#' @param main,xlab,ylab,type Plotting parameters with useful defaults.
#' @param ... Further plotting parameters.
#'
#' @details
#' This function provides density estimates of size-biased data according to
#' Montoril et al. (2021). It estimates the power of the density of interest by
#' wavelet basis, i.e., \eqn{\hat{f^a}, a \ge 0.5}. After that, the density is
#' finally estimated by \eqn{[\hat{f^a}]^{1/a}}, where the power is taken of the
#' positive part \eqn{\max(0, \hat{f^a})}, except for \eqn{a = 0.5}, whose
#' square is never negative. Periodized version of ordinary
#' orthonormal wavelet bases of families \emph{Daublets}, \emph{Symmlets} and
#' \emph{Coiflets} are allowed, as well as warped bases of these families. It
#' also accepts wavelet filters provided in \code{wavelet.filter} (see the help
#' of \command{\link{wbasis}}).
#'
#' The elements of the \code{data} must all be finite numbers. Otherwise, they
#' will be removed before the density estimation procedure. In this case, a
#' warning message will be released.
#'
#' For stability reasons, the values at which the density estimates are to be
#' computed must satisfy \eqn{\min(data) \le from \le to \le \max(data)}.
#'
#' By default, the argument \code{dens.biased} provides the kernel density
#' estimates of the \code{data}. It is computed by
#' \command{\link[stats]{density.default}}, where the kernel adopted is the
#' Gaussian and the bandwidth selection used is the "\code{SJ}" proposed by
#' Sheather and Jones, (1991).
#'
#' For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
#' \code{prec.wavelet} and \code{wavelet.filter}, see the help of
#' \command{\link{wbasis}}.
#'
#' The arguments \code{warp.fun} and \code{dwarp.fun} are used only when
#' \code{warped = TRUE}. They must be functions representing a cumulative
#' distribution function and a probability density function, respectively. By
#' default, the function uses, respectively, the empirical distribution function
#' (EDF) and the kernel density estimate (KDE) of the \code{data}. The EDF is
#' computed by \command{\link[stats]{ecdf}}. The KDE is computed by
#' \command{\link[stats]{density.default}}, where the kernel adopted is the
#' Gaussian and the bandwidth selection used is the \code{"SJ"} proposed by
#' Sheather and Jones, (1991). Furthermore, when \code{warped = TRUE}, the
#' \code{data} are not rescaled, because they are scaled by the warping function
#' to be used.
#'
#' When \code{warped = FALSE}, the \code{data} will be scaled when
#' \code{rescale = TRUE}. In this case, a linear transformation is applied and
#' the observations will belong to the interval \eqn{[eps, 1-eps]}. By default,
#' \eqn{eps = 1.9^{-J1}}, following the suggestion by Montoril et al. (2019).
#' When \code{rescale = FALSE}, the original \code{data} are to be used.
#'
#' \bold{Caution: as the wavelet basis is periodized in the unit interval, if
#' \eqn{\max(data) > 1} or \eqn{\min(data) < 0}, DO NOT use \code{warped = FALSE}
#' with \code{rescale = FALSE}. Otherwise, the estimates will be unfeasible.}
#'
#' When \code{correct = TRUE}, the estimate \eqn{\hat{f}} is replaced by
#' \eqn{\tilde{f} = \max(0, \hat{f} - \xi)}, the correction of Glad, Hjort and
#' Ushakov (2003), where the constant \eqn{\xi} makes \eqn{\tilde{f}} integrate
#' to one over the domain \eqn{D} of the estimator. With \code{warped = FALSE},
#' \eqn{D} is the interval of the basis, \eqn{[0, 1]}, in the scale of the
#' data: \eqn{[\min(data) - c, \max(data) + c]}, with
#' \eqn{c = eps (\max(data) - \min(data))/(1 - 2 eps)}, if
#' \code{rescale = TRUE}, and \eqn{[0, 1]} itself otherwise. With
#' \code{warped = TRUE}, \eqn{D} is the range of the \code{data}. The domain
#' does not depend on \code{from}, \code{to} and \code{length.obs}, which only
#' set the points where the corrected estimate is evaluated: it integrates to
#' one over \eqn{[}\code{from}\eqn{, }\code{to}\eqn{]} only when this interval
#' covers \eqn{D}.
#'
#' If \eqn{\hat{f}} integrates to more than one over \eqn{D}, then
#' \eqn{\xi > 0}: the estimate is lowered by \eqn{\xi}, and set to zero where
#' it falls below \eqn{\xi}. This never increases the integrated squared error
#' over \eqn{D} (Theorem 1 of Glad et al., 2003). Otherwise, \eqn{\xi \le 0}
#' and \eqn{\tilde{f} = \hat{f} - \xi} on \eqn{D}: the missing mass is spread
#' evenly over \eqn{D}, the counterpart on a bounded domain of their second
#' correction. It does not increase the integrated squared error either when
#' the density of interest vanishes outside \eqn{D}. In both cases,
#' \eqn{\tilde{f}} is the projection of \eqn{\hat{f}} on the densities on
#' \eqn{D}, in \eqn{L_2(D)}. With \code{power.dens = 0.5} and
#' \code{warped = FALSE}, the estimate is already a density on \eqn{D}, since
#' its coefficients are normalized, and \eqn{\xi = 0}. The integral is computed
#' numerically, on the dyadic points of \eqn{[0, 1]}, 16 per interval of the
#' level \code{J1} (fewer when \code{J1 > 18}), where the expansion is
#' evaluated by a convolution of its coefficients with the scaling function.
#'
#' @return A list containing the following components.
#'   \item{y}{The grid of \code{length.obs} observations at which the densities
#'   are estimated.}
#'   \item{dens}{The estimated density values.}
#'   \item{xi}{Only when \code{correct = TRUE}: the constant \eqn{\xi} of the
#'   correction (see Details).}
#'
#' @references
#' Glad, I. K., Hjort, N. L. and Ushakov, N. G. (2003). Correction of density
#' estimators that are not densities. \emph{Scandinavian Journal of
#' Statistics}, 30(2), 415-427, \doi{10.1111/1467-9469.00339}.
#'
#' Montoril, M. H., Chang, W. and Vidakovic, B. (2019). Wavelet-Based Estimation
#' of Generalized Discriminant Functions. \emph{Sankhya B}, 81(2), 318-349,
#' \doi{10.1007/s13571-018-0158-1}.
#'
#' Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2021). Wavelet-based
#' estimation of power densities of size-biased data. \emph{arXiv:2112.12895
#' [math, stat]}. \url{https://arxiv.org/abs/2112.12895}.
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation. \emph{Journal of the Royal
#' Statistical Society series B}, 53, 683-690.
#' \url{http://www.jstor.org/stable/2345597}.
#'
#' @seealso \command{\link{PHI}}, \command{\link{bac}}, \command{\link{wtable}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' # Density estimation of the Blood alcohol concentration of drivers
#' # involved in fatal accidents in the USA in 2019.
#' # This example reproduces the application in Montoril et al. (2021).
#'
#' # Loading the data
#' data(bac)
#'
#' # Let's consider four different possibilities of estimation:
#' # m1 - non-warped wavelets, non-rescaled data and power.density
#' #      equals to half
#' m1 = wdensity(data = bac, wf = function(x) 0.1+.9*x, power.dens = 0.5,
#'               J1 = ceiling(0.45*log2(length(bac))), family = "s",
#'               filter.size = 20, warped = FALSE, rescale = FALSE,
#'               thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
#'               col = 1, lty = 1, main = expression(m[1]))
#'
#' # m2 - non-warped wavelets, non-rescaled data and power.density
#' #      equals to one
#' m2 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 1,
#'               J1 = ceiling(0.45*log2(length(bac))), family = "s",
#'               filter.size = 20, warped = FALSE, rescale = FALSE,
#'               thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
#'               col = 1, lty = 1, main = expression(m[2]))
#'
#' # m3 - warped wavelets and power.density equals to half
#' m3 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 0.5,
#'               J1 = ceiling(0.95*log2(length(bac))), family = "s",
#'               filter.size = 20, warped = TRUE, rescale = FALSE,
#'               thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
#'               col = 1, lty = 1, main = expression(m[3]))
#'
#' # m4 - warped wavelets and power.density equals to one
#' m4 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 1,
#'               J1 = ceiling(0.95*log2(length(bac))), family = "s",
#'               filter.size = 20, warped = TRUE, rescale = FALSE,
#'               thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
#'               col = 1, lty = 1, main = expression(m[4]))
#'
#' @keywords smooth
#' @importFrom stats approx density.default mad
#' @export
wdensity <- function(data, from, to, length.obs = 250, wf = NULL, power.dens,
                     J0 = 0, J1, family = "Daublets", filter.size = 10,
                     prec.wavelet = 30, wavelet.filter, wavelet.table = NULL,
                     dens.biased, warped = TRUE, warp.fun, dwarp.fun,
                     rescale = TRUE, eps = 1.9^(-J1), thresh = "hard",
                     plot = TRUE, main, xlab, ylab, type,
                     boundary = c("periodic", "interval"),
                     check.lawton = TRUE, correct = FALSE, ...){

  boundary <- match.arg(tolower(boundary[1L]), c("periodic", "interval"))
  if(boundary == "interval" && !is.null(wavelet.table))
    warning("wavelet.table is not used with boundary = \"interval\" and was ignored.")

  # taken before data is changed below, for the default xlab
  data.name <- deparse(substitute(data))

  if (!all(is.finite(data))){
    data <- data[is.finite(data)]
    if(!length(data))
      stop("'data' has no finite element.")
    warning("At least one element of 'data' is not finite. Non-finite elements were removed.")
  }

  # eps is only used by the rescaling of the unwarped estimator.
  if(!warped && rescale){
    if(length(eps) != 1L)
      stop("'eps' must be a single value.")
    .check_eps(eps)
  }

  thresh <- tolower(substring(thresh, 1, 1))
  if(!thresh %in% c("l","h","s"))
    stop("The options available for thresh are 'linear', 'hard' and 'soft'.")

  if(missing(power.dens))
    power.dens <- 0.5
  if(power.dens < 0.5)
    stop("power.dens must be greater than or equal to 0.5.")

  ry <- range(data)
  sy <- sort(data)
  n <- length(data)

  if(missing(from))
    from <- sy[1]
  if(missing(to))
    to <- sy[n]

  if(from > to)
    stop("'from' must be smaller than or equal to 'to'.")
  if(missing(dens.biased)){
    dens.bias <- density.default(x = sy, bw = "sj", from = sy[1], to = sy[n])
    dens.biased <- approx(x = dens.bias$x, y = dens.bias$y, xout = sy)$y
  }

  y.vals <- seq(from = from, to = to, length.out = length.obs)

  # The warping function, for the domain of the correction
  Hfun <- NULL
  if(warped){
    scale <- 1
    if(missing(warp.fun)){
      Hyc <- seq_len(length.out = n)/n
      hy <- dens.biased
      Hy.vals <- approx(x = sy, y = Hyc, xout = y.vals)$y
      Hfun <- function(y) approx(x = sy, y = Hyc, xout = y, ties = mean)$y
    }
    else if(!is.function(warp.fun))
      stop("warp.fun must be a function.")
    else if(missing(dwarp.fun))
      stop("dwarp.fun must be a function.")
    else if(!is.function(dwarp.fun))
      stop("dwarp.fun must be a function.")
    else{
      Hyc <- warp.fun(sy)
      hy <- dwarp.fun(sy)
      Hy.vals <- warp.fun(y.vals)
      Hfun <- warp.fun
    }
  }
  else if(rescale){
    a <- eps*diff(ry)/(1 - 2*eps)
    location <- ry[1] - a
    scale <- 2*a + diff(ry)
    Hyc <- (sy - location)/scale
    hy <- rep_len(x = 1/scale, length.out = n)
    Hy.vals <- (y.vals - location)/scale
  }
  else{
    Hyc <- sy
    scale <- 1
    hy <- rep_len(x = 1, length.out = n)
    Hy.vals <- y.vals
  }

  mbasis <- PHI(x = c(Hyc, Hy.vals), J = J1, family = family,
                filter.size = filter.size, prec.wavelet = prec.wavelet,
                boundary = boundary, wavelet.filter = wavelet.filter,
                wavelet.table = if(boundary == "interval") NULL else wavelet.table,
                check.lawton = check.lawton)

  matvals <- mbasis[seq_len(n),]*dens.biased^(power.dens-1)*hy/wf(sy)^power.dens
  coefs <- (scale/mean(1/wf(sy)))^(power.dens)*colMeans(matvals)

  # The thresholding transforms use the filter of the projection. A filter
  # given by the user goes as "Own": PHI has already checked it, and warned
  # if the family was a built-in one.
  wargs <- if(missing(wavelet.filter))
             list(family = family, filter.size = filter.size)
           else
             list(family = "Own", wavelet.filter = wavelet.filter)

  dec <- function(v) .wb_quiet_filter(
    do.call(wavedec, c(list(x = v, j0 = J0, boundary = boundary), wargs)))
  rec <- function(v) .wb_quiet_filter(
    do.call(waverec, c(list(x = v, j0 = J0, boundary = boundary), wargs)))

  if(thresh == "h"){
    coefs <- dec(coefs)
    sig <- mad(coefs[-seq_len(length.out = 2^(J1-1))])
    lambda <- sig*sqrt(2*log(2^(J1-1)))
    d <- coefs[-seq_len(length.out = 2^J0)]
    coefs[-seq_len(length.out = 2^J0)] <- d*(abs(d) > lambda)
    coefs <- rec(coefs)
  }
  if(thresh == "s"){
    coefs <- dec(coefs)
    sig <- mad(coefs[-seq_len(length.out = 2^(J1-1))])
    lambda <- sig*sqrt(2*log(2^(J1-1)))
    d <- coefs[-seq_len(length.out = 2^J0)]
    coefs[-seq_len(length.out = 2^J0)] <- sign(d)*(abs(d) - lambda)*(abs(d) > lambda)
    coefs <- rec(coefs)
  }

  if(!warped & power.dens == 0.5)
    coefs <- coefs/sqrt(sum(coefs^2))

  # The power of the positive part: a negative estimate of f^a has no real
  # power 1/power.dens (except at power.dens = 1, where it stays negative),
  # and gave NaN. At 0.5, the square is never negative and needs no clipping.
  fpow <- function(g)
    if(power.dens == 0.5) g^2 else pmax(g, 0)^(1/power.dens)

  dens.vals <- fpow(drop(mbasis[-seq_len(n),]%*%coefs))/scale

  if(correct){
    xi <- if(!warped && power.dens == 0.5) 0 else
      .wdensity_xi(coefs = coefs, J1 = J1, fpow = fpow, scale = scale,
                   Hfun = Hfun, sy = sy, family = family,
                   filter.size = filter.size, prec.wavelet = prec.wavelet,
                   wavelet.filter = if(missing(wavelet.filter)) NULL
                                    else wavelet.filter,
                   boundary = boundary)
    dens.vals <- pmax(dens.vals - xi, 0)
  }

  if(plot){
    if(missing(main))
      main <- paste("Density estimates using power =", power.dens)
    if(missing(xlab))
      xlab <- data.name
    if(missing(ylab))
      ylab <- "Density"
    if(missing(type))
      type <- "l"

    plot(x = y.vals, y = dens.vals, type = type, main = main,
         xlab = xlab, ylab = ylab, ...)
  }

  out <- list(y = y.vals, dens = dens.vals)
  if(correct)
    out$xi <- xi

  return(invisible(out))

}

# The constant xi of the correction of Glad, Hjort and Ushakov (2003) for the
# estimate fpow(g)/scale, g = sum_k coefs[k] phi_{J1,k-1}, over the domain of
# the estimator: [0, 1] in the scale of the basis without warping (Hfun =
# NULL), and the range of the sorted data sy with the warping Hfun. The
# expansion is evaluated on the dyadic grid of [0, 1] with 2^r points per
# interval of level J1. Without warping, the grid is the quadrature (in u,
# times the jacobian scale). With warping, the integral is in the scale of the
# data: the nodes are a uniform grid of the range, the data and the preimages
# of the dyadic grid, so that they follow the resolution of the warped basis,
# and g is interpolated between the dyadic points.
.wdensity_xi <- function(coefs, J1, fpow, scale, Hfun, sy, family,
                         filter.size, prec.wavelet, wavelet.filter, boundary){

  r <- min(4L, max(0L, 22L - J1))
  M <- 2^(J1 + r)

  fam <- if(is.null(wavelet.filter))
           which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
         else 4L
  gd <- .wb_fine_grid(coefs, J = J1, r = r, fam = fam,
                      filter.size = filter.size, prec.wavelet = prec.wavelet,
                      wavelet.filter = if(is.null(wavelet.filter)) 0
                                       else wavelet.filter,
                      boundary = boundary)

  if(is.null(Hfun)){
    fq <- fpow(gd)/scale
    wq <- c(0.5, rep_len(1, M - 1), 0.5)*scale/M
  }
  else{
    n <- length(sy)
    ud <- (0:M)/M
    yk <- sort(unique(c(seq(from = sy[1], to = sy[n], length.out = M + 1),
                        sy)))
    Hk <- Hfun(yk)
    ui <- ud[ud > Hk[1] & ud < Hk[length(Hk)]]
    yq <- sort(unique(c(yk, approx(x = Hk, y = yk, xout = ui,
                                   ties = mean)$y)))
    fq <- fpow(approx(x = ud, y = gd, xout = Hfun(yq))$y)
    dq <- diff(yq)
    wq <- (c(dq, 0) + c(0, dq))/2
  }

  .wb_bonafide_xi(fq, wq)
}

# Values of sum_k coefs[k] phi_{J,k-1}(u), in the basis of PHI(), at the
# dyadic points u = l/2^(J + r), l = 0, ..., 2^(J + r), computed in C without
# a basis matrix. fam is the code of the family (4 for wavelet.filter).
.wb_fine_grid <- function(coefs, J, r, fam, filter.size, prec.wavelet,
                          wavelet.filter = 0, boundary = "periodic"){
  bcode <- .wb_boundary_code(boundary)
  cdv <- if(bcode == 2L)
           .cdv_prepare(c(0, 1), fam, filter.size, level = J, what = "basis")
  .Call("_WaveBased_C_WavFineGrid", as.double(coefs), as.integer(J),
        as.integer(r), as.integer(fam), as.integer(filter.size),
        as.integer(prec.wavelet), as.double(wavelet.filter), bcode, cdv)
}

# The xi solving sum(w * pmax(f - xi, 0)) == 1, computed in C.
.wb_bonafide_xi <- function(f, w)
  .Call("_WaveBased_C_BonaFideXi", as.double(f), as.double(w))
