#' @title Wavelet reconstruction
#'
#' @description Reconstructs a signal using the Mallat's pyramidal algorithm in
#' the reconstruction stage (Mallat, 1989).
#'
#' @param x Vector containing the decomposed data, or a matrix whose rows are
#'   decomposed signals. The length of the vector, or the number of columns of
#'   the matrix, must be power of 2. A matrix with a single column is taken as
#'   a vector.
#' @param j0 The coarsest resolution level where the signal is decomposed.
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
#' This function applies the inverse discrete wavelet transform to the
#' (decomposed) data \code{x}, whose length must be a power of 2. In this case,
#' it is assumed that the data is decomposed with the coarsest level \code{j0}.
#'
#' If \code{x} is a matrix, each row is reconstructed, as by
#' \code{t(apply(x, 1, waverec, ...))}, but in a single call to the compiled
#' code, which sets up the filter (and, with \code{boundary = "interval"}, the
#' boundary blocks) only once.
#'
#' For the moment, the function deal with boundary conditions by applying
#' periodic extensions to the data. In the future, a reflective extesion may be
#' taken into account.
#'
#' For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
#' \code{prec.wavelet} and \code{wavelet.filter}, see the help of
#' \command{\link{wbasis}}.
#'
#' @return A vector containing a reconstruction of the decomposed data \code{x},
#'   whose coarsest resolution level is \code{j0}. If \code{x} is a matrix, a
#'   matrix of the same dimensions, whose rows hold the reconstruction of the
#'   rows of \code{x}, with the same row names.
#'
#' @references
#' Mallat, S.G. (1989). A theory for multiresolution signal decomposition:
#' the wavelet representation. \emph{IEEE Trans. Patt. Anal. and Mach.
#' Intell.}, 11, 674-693, \doi{10.1109/34.192463}.
#'
#' @seealso \code{\link{wavedec}}
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
#' # Let's apply the wavelet decomposition to x using a 18-tap
#' # filter of a Daulet.
#' # Consider j0 = 3 as the coarsest resolution level:
#' wdx <- wavedec(x, j0 = 3, family = "daublets", filter.size = 18)
#'
#' # Now, let's apply the wavelet reconstrution to wdx.
#' wrx <- waverec(wdx, j0 = 3, family = "daublets", filter.size = 18)
#'
#' # The reconstructed vector should be equal to x (or, at least,
#' # as close as possible)
#' sum(abs(wrx-x))
#'
#' # The rows of a matrix are reconstructed in a single call:
#' X <- rbind(x, rev(x))
#' WX <- wavedec(X, j0 = 3, family = "daublets", filter.size = 18)
#' all.equal(waverec(WX, j0 = 3, family = "daublets", filter.size = 18), X)
#'
#' @keywords smooth
#' @export
waverec <- function(x, j0 = 0, family = "Daublets", filter.size = 20,
                    wavelet.filter, boundary = c("periodic", "interval")){

  if(is.list(x))
    stop("'x' must be a numeric vector or matrix.")

  # a matrix is reconstructed by rows; a matrix with a single column, or an
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

  wrec <- .Call("_WaveBased_C_WaveRec", xd, as.integer(fam),
                as.integer(filter.size), as.integer(j0),
                as.double(wavelet.filter), bcode, cdv)
  if(rows)
    rownames(wrec) <- rownames(x)

  return(wrec)

}
