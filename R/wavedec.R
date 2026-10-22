#' @title Wavelet decomposition
#'
#' @description Decomposes a signal using the Mallat's pyramidal algorithm
#' (Mallat, 1989).
#'
#' @param x Vector containing the data to be decomposed, or a matrix whose rows
#'   are the signals to be decomposed. The length of the vector, or the number
#'   of columns of the matrix, must be power of 2. A matrix with a single
#'   column is taken as a vector.
#' @param j0 The coarsest resolution level where the signal is to be decomposed.
#' @param family The family of wavelets to use. It can be \emph{Daublets},
#'   \emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter
#'   is provided \code{wavelet.filter}.
#' @param filter.size The size of the filter. The available values depend on the
#'   chosen family. See Details for more information.
#' @param wavelet.filter Use this to provide your own filter, the low-pass
#'   filter of an orthonormal wavelet basis. To use this argument, you must
#'   specify \code{family = "Own"}. The filter is checked: see Details of
#'   \command{\link{wbasis}}.
#' @param boundary The boundary treatment of the transform. Either
#'   \code{"periodic"} (the default, the usual periodic wavelet transform)
#'   or \code{"interval"} (the boundary-corrected transform of Cohen,
#'   Daubechies and Vial, 1993, which requires a Daublet or Symmlet filter
#'   and a sufficiently large coarsest level -- an informative error states
#'   the exact minimum). See \command{\link{wbasis}} for a detailed
#'   description of the interval construction.
#'
#' @details
#' This function applies the discrete wavelet transform to the data \code{x},
#' whose length must be a power of 2. As result, the function returns
#' \deqn{c_{j_0,0} \ldots c_{j_0,2^{j_0}-1} d_{j_0,0} \ldots d_{j_0,2^{j_0}-1}
#' \ldots d_{J-1,0} \ldots d_{J-1,2^{J-1}-1},}
#' where \eqn{c_{j_0,k}} represent the scale coefficients in the coarsest level
#' and \eqn{d_{j,k}} represent the detail coefficients in level \eqn{j}.
#'
#' If \code{x} is a matrix, each row is decomposed, as by
#' \code{t(apply(x, 1, wavedec, ...))}, but in a single call to the compiled
#' code, which sets up the filter (and, with \code{boundary = "interval"}, the
#' boundary blocks) only once. This is not the two-dimensional transform of an
#' image.
#'
#' For more details about \code{family}, \code{filter.size},
#' \code{wavelet.filter} and \code{boundary}, see the help of
#' \command{\link{wbasis}}.
#'
#' @return A vector of wavelet coefficients decomposed from \code{x}, where the
#'   coarsest resolution level is \code{j0}. If \code{x} is a matrix, a matrix
#'   of the same dimensions, whose rows hold the coefficients of the rows of
#'   \code{x}, with the same row names.
#'
#' @references
#' Mallat, S.G. (1989). A theory for multiresolution signal decomposition:
#' the wavelet representation. \emph{IEEE Trans. Patt. Anal. and Mach.
#' Intell.}, 11, 674-693, \doi{10.1109/34.192463}.
#'
#' @seealso \code{\link{waverec}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' #
#' # A simple example
#' #
#' set.seed(123)
#' n <- 64
#' x <- sort(runif(n))
#'
#'
#' # Let's apply the wavelet decomposition to x using a 18-tap
#' # filter of a Daulet.
#' # Consider j0 = 3 as the coarsest resolution level:
#' wdx <- wavedec(x, j0 = 3, family = "daublets", filter.size = 18)
#'
#' # The rows of a matrix are decomposed in a single call:
#' X <- rbind(x, rev(x))
#' WX <- wavedec(X, j0 = 3, family = "daublets", filter.size = 18)
#' all.equal(WX[1, ], wdx)
#'
#' @keywords smooth
#' @export
wavedec <- function(x, j0 = 0, family = "Daublets", filter.size = 20,
                    wavelet.filter, boundary = c("periodic", "interval")){

  if(is.list(x))
    stop("'x' must be a numeric vector or matrix.")

  # a matrix is decomposed by rows; a matrix with a single column, or an
  # array with a single dimension of extent greater than one, is a vector
  rows <- is.matrix(x) && ncol(x) > 1L
  if(!rows){
    if(sum(dim(x) > 1L) > 1L)
      stop("'x' must be a vector or a matrix, not an array.")
    x <- as.vector(x)
  }

  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

  if(any(!is.finite(x)))
    stop("'x' must be a numeric vector.")

  n <- if(rows) ncol(x) else length(x)
  if(n < 2L || bitwAnd(n, n - 1L) != 0L)
    stop("The length of 'x', or its number of columns, must be a power of ",
         "two, and at least two (got ", n, ").")

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
    .wb_check_filter(wavelet.filter)
  }

  bcode <- .wb_boundary_code(boundary)
  if(bcode == 0L)
    stop("boundary must be \"periodic\" or \"interval\" for the wavelet transform.")

  xd <- if(rows) matrix(as.double(x), nrow(x), ncol(x)) else as.double(x)

  cdv <- NULL
  if(bcode == 2L)
    cdv <- .cdv_prepare(xd, fam, filter.size,
                        level = j0, what = "decompose", check.range = FALSE)

  wdec <- .Call("_WaveBased_C_WaveDec", xd, as.integer(fam),
                as.integer(filter.size), as.integer(j0),
                as.double(wavelet.filter), bcode, cdv)
  if(rows)
    rownames(wdec) <- rownames(x)

  return(wdec)

}
