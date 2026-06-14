#' @title Wavelet decomposition
#'
#' @description Decomposes a signal using the Mallat's pyramidal algorithm
#' (Mallat, 1989).
#'
#' @param x Vector containing the data to be decomposed. Its length must be
#'   power of 2.
#' @param j0 The coarsest resolution level where the signal is to be decomposed.
#' @param family The family of wavelets to use. It can be \emph{Daublets},
#'   \emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter
#'   is provided \code{wavelet.filter}.
#' @param filter.size The size of the filter. The available values depend on the
#'   chosen family. See Details for more information.
#' @param wavelet.filter Use this to provide your own filter. To use this
#'   argument, you must specify \code{family = "Own"}. Do not use it, if you are
#'   not sure about what you are doing.
#'
#' @details
#' This function applies the discrete wavelet transform to the data \code{x},
#' whose length must be a power of 2. As result, the function returns
#' \deqn{c_{j_0,0} \ldots c_{j_0,2^{j_0}-1} d_{j_0,0} \ldots d_{j_0,2^{j_0}-1}
#' \ldots d_{J-1,0} \ldots d_{J-1,2^{J-1}-1},}
#' where \eqn{c_{j_0,k}} represent the scale coefficients in the coarsest level
#' and \eqn{d_{j,k}} represent the detail coefficients in level \eqn{j}.
#'
#' For the moment, the function deal with boundary conditions by applying
#' periodic extensions to the data. In the future, a reflective extesion may be
#' taken into account.
#'
#' For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
#' \code{prec.wavelet} and \code{wavelet.filter}, see the help of
#' \command{\link{wbasis}}.
#'
#' @return A vector of wavelet coefficients decomposed from \code{x}, where the
#'   coarsest resolution level is \code{j0}.
#'
#' @references
#' Mallat, S.G. (1989). A theory for multiresolution signal decomposition:
#' the wavelet representation. \emph{IEEE Trans. Patt. Anal. and Mach.
#' Intell.}, 11, 674-693, \url{https://ieeexplore.ieee.org/document/192463}.
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
#' @keywords smooth
#' @export
wavedec <- function(x, j0 = 0, family = "Daublets", filter.size = 20,
                    wavelet.filter){

  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }

  if(any(!is.finite(x)))
    stop("'x' must be a numeric vector.")

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

  wdec <- .Call("_WaveBased_C_WaveDec", as.double(x), as.integer(fam),
                as.integer(filter.size), as.integer(j0),
                as.double(wavelet.filter))

  return(wdec)

}
