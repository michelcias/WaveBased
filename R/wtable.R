#' @title Interpolation Tables for the Scaling and Wavelet Functions
#'
#' @description Precomputes the scaling function \eqn{\phi} and the mother
#' wavelet \eqn{\psi} on a fine regular grid, so that \code{\link{PHI}},
#' \code{\link{PSI}}, \code{\link{wbasis}} and \code{\link{wdensity}} can
#' replace the Daubechies-Lagarias algorithm by a fast table lookup with
#' linear interpolation.
#'
#' @param family The family of wavelets to use. It can be \emph{Daublets},
#'   \emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter
#'   is provided in \code{wavelet.filter}.
#' @param filter.size The size of the filter. The available values depend on
#'   the chosen family. See the help of \command{\link{wbasis}} for more
#'   information.
#' @param prec.wavelet The number of iterations of the Daubechies-Lagarias
#'   algorithm used to compute the exact values stored in the table.
#' @param ngrid Number of grid intervals in the unit interval. The grid step is
#'   \code{1/ngrid}. The default, \code{2^13}, keeps the interpolation error of
#'   smooth wavelets (e.g. \code{filter.size >= 8}) below \code{1e-5}.
#' @param wavelet.filter Use this to provide your own filter. To use this
#'   argument, you must specify \code{family = "Own"}. Do not use it, if you
#'   are not sure about what you are doing.
#' @param check Logical. If \code{TRUE} (default), the interpolation error is
#'   measured against the exact Daubechies-Lagarias values at off-grid points
#'   and stored in the returned object.
#' @param check.points Number of off-grid points used when \code{check = TRUE}.
#' @param check.tol Tolerance for the measured error. A warning is released
#'   when the measured error exceeds it.
#'
#' @details
#' The Daubechies-Lagarias algorithm evaluates \eqn{\phi} and \eqn{\psi} at a
#' point through \code{prec.wavelet} products of matrices of dimension
#' \eqn{(N-1)\times(N-1)}, where \eqn{N} is the filter size. Both functions,
#' however, depend on the evaluation point only through its fractional part.
#' This function tabulates the exact values on the regular grid
#' \eqn{w = g/\code{ngrid}}, \eqn{g = 0, \ldots, \code{ngrid}}, producing two
#' \eqn{(N-1)\times(\code{ngrid}+1)} matrices that are \strong{independent of
#' the resolution level} \eqn{J}: a single table serves \code{\link{PHI}},
#' \code{\link{PSI}} and \code{\link{wbasis}} at any level.
#'
#' The table for a 20-tap filter with the default grid occupies about 2.5 MB
#' and is generated in a fraction of a second, so it can be created once per
#' session and reused. If desired, it can be stored on disk with
#' \code{\link{saveRDS}} and restored with \code{\link{readRDS}}.
#'
#' The accuracy of the linear interpolation depends on the regularity of the
#' wavelet, which increases with the filter size. For \code{filter.size >= 8}
#' the error of the interpolated \eqn{\phi} and \eqn{\psi} values with the
#' default grid stays below \code{1e-5} (five decimal places). Note that the
#' entries of the basis matrices returned by \code{\link{PHI}},
#' \code{\link{PSI}} and \code{\link{wbasis}} are scaled by \eqn{2^{J/2}},
#' which amplifies this error accordingly at high resolution levels. Very
#' short filters (e.g. Haar or the 4-tap Daublet) generate functions of low
#' regularity, for which interpolation is not recommended; with
#' \code{check = TRUE}, such cases are reported by a warning.
#'
#' @return An object of class \code{"wavelet_table"}: a list with components
#'   \item{phi}{Matrix \eqn{(N-1)\times(\code{ngrid}+1)} with the tabulated
#'   scaling function vectors.}
#'   \item{psi}{Matrix \eqn{(N-1)\times(\code{ngrid}+1)} with the tabulated
#'   wavelet vectors.}
#'   \item{family}{Canonical family name.}
#'   \item{fam}{Internal family code.}
#'   \item{filter.size}{The size of the filter.}
#'   \item{prec.wavelet}{Precision used to build the table.}
#'   \item{ngrid}{Number of grid intervals.}
#'   \item{wavelet.filter}{The custom filter, when \code{family = "Own"}.}
#'   \item{max.error}{Measured interpolation error (\code{NA} when
#'   \code{check = FALSE}).}
#'
#' @seealso \command{\link{PHI}}, \command{\link{PSI}},
#'   \command{\link{wbasis}}, \command{\link{wdensity}}
#'
#' @author Michel H. Montoril \email{michel@@ufscar.br}
#'
#' @examples
#' # Build a table for the 20-tap Symmlet and reuse it at different levels
#' tab <- wtable(family = "symmlets", filter.size = 20)
#' tab
#'
#' set.seed(123)
#' x <- sort(runif(100))
#'
#' # Interpolated versions: same call, extra argument
#' w5 <- wbasis(x, j0 = 0, J = 5, family = "symmlets", filter.size = 20,
#'              wavelet.table = tab)
#' w8 <- wbasis(x, j0 = 0, J = 8, family = "symmlets", filter.size = 20,
#'              wavelet.table = tab)
#'
#' # Agreement with the exact Daubechies-Lagarias evaluation
#' w5.exact <- wbasis(x, j0 = 0, J = 5, family = "symmlets", filter.size = 20)
#' max(abs(w5 - w5.exact))
#'
#' @keywords smooth
#' @export
wtable <- function(family = "Daublets", filter.size = 20, prec.wavelet = 30,
                   ngrid = 2^13, wavelet.filter, check = TRUE,
                   check.points = 1000, check.tol = 1e-5){

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

  if(ngrid < 2)
    stop("'ngrid' must be an integer greater than or equal to 2.")

  N <- if(fam == 4) length(wavelet.filter) else filter.size
  if(N < 6)
    warning("Interpolation tables are not recommended for very short filters (low regularity). Consider the exact evaluation instead.")

  tabs <- .Call("_WaveBased_C_WavTable", as.integer(fam),
                as.integer(filter.size), as.integer(prec.wavelet),
                as.double(wavelet.filter), as.integer(ngrid))

  out <- structure(
    list(phi = tabs[[1]],
         psi = tabs[[2]],
         family = c("Daublets", "Symmlets", "Coiflets", "Own")[fam],
         fam = fam,
         filter.size = if(fam == 4) length(wavelet.filter) else filter.size,
         prec.wavelet = prec.wavelet,
         ngrid = ngrid,
         wavelet.filter = if(fam == 4) wavelet.filter else NULL,
         max.error = NA_real_),
    class = "wavelet_table")

  if(check){
    # Deterministic off-grid points, so the user's RNG state is not touched.
    # The comparison uses J = 0 and periodic = FALSE, so the measured error
    # refers to the function values themselves (no sqrt(2^J) scaling).
    xs <- ((seq_len(check.points) - 0.5)/check.points + pi/10^4) %% 1
    payload <- list(out$phi, out$psi)

    exact.phi <- .Call("_WaveBased_C_PHImat", as.double(xs), 0L,
                       as.integer(fam), as.integer(filter.size),
                       as.integer(prec.wavelet), 0L,
                       as.double(wavelet.filter), NULL)
    apx.phi <- .Call("_WaveBased_C_PHImat", as.double(xs), 0L,
                     as.integer(fam), as.integer(filter.size),
                     as.integer(prec.wavelet), 0L,
                     as.double(wavelet.filter), payload)

    exact.psi <- .Call("_WaveBased_C_PSImat", as.double(xs), 0L,
                       as.integer(fam), as.integer(filter.size),
                       as.integer(prec.wavelet), 0L,
                       as.double(wavelet.filter), NULL)
    apx.psi <- .Call("_WaveBased_C_PSImat", as.double(xs), 0L,
                     as.integer(fam), as.integer(filter.size),
                     as.integer(prec.wavelet), 0L,
                     as.double(wavelet.filter), payload)

    out$max.error <- max(abs(exact.phi - apx.phi), abs(exact.psi - apx.psi))

    if(out$max.error > check.tol)
      warning("The measured interpolation error (", format(out$max.error, digits = 3),
              ") exceeds ", format(check.tol, digits = 3),
              ". Increase 'ngrid' or use the exact evaluation.")
  }

  out
}

#' @rdname wtable
#' @param x An object of class \code{"wavelet_table"}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.wavelet_table <- function(x, ...){
  cat("Wavelet interpolation table\n")
  cat("  Family       :", x$family, "\n")
  cat("  Filter size  :", x$filter.size, "\n")
  cat("  Grid step    : 1/", x$ngrid, "\n", sep = "")
  cat("  Memory       : ",
      format((length(x$phi) + length(x$psi))*8/2^20, digits = 3), " MB\n",
      sep = "")
  if(is.na(x$max.error))
    cat("  Max. error   : not measured (check = FALSE)\n")
  else
    cat("  Max. error   : ", format(x$max.error, digits = 3),
        " (vs Daubechies-Lagarias)\n", sep = "")
  invisible(x)
}

#' @rdname wtable
#' @importFrom graphics abline par
#' @export
plot.wavelet_table <- function(x, ...){

  N <- nrow(x$phi) + 1
  G <- x$ngrid

  # phi(u), u in [0, N - 1]: the tabulated entry T[m, g] holds
  # phi(g/G + N - 2 - m), so offset c = N - 2 - m corresponds to row N - 1 - c
  u.phi <- seq(0, N - 1, by = 1/G)
  v.phi <- c(unlist(lapply(0:(N - 2),
                           function(cc) x$phi[N - 1 - cc, seq_len(G)])),
             x$phi[1, G + 1])

  # psi(u), u in [-(N/2 - 1), N/2]: T[m, g] holds psi(g/G + N/2 - 1 - m),
  # so offset c = N/2 - 1 - m corresponds to row N/2 - c
  u.psi <- seq(-(N/2 - 1), N/2, by = 1/G)
  v.psi <- c(unlist(lapply((-N/2 + 1):(N/2 - 1),
                           function(cc) x$psi[N/2 - cc, seq_len(G)])),
             x$psi[1, G + 1])

  old.par <- par(mfrow = c(1, 2))
  on.exit(par(old.par))

  plot(u.phi, v.phi, type = "l",
       main = expression(paste("Scaling function ", phi)),
       xlab = "x", ylab = expression(phi(x)), ...)
  abline(h = 0, lty = 3)

  plot(u.psi, v.psi, type = "l",
       main = expression(paste("Wavelet ", psi)),
       xlab = "x", ylab = expression(psi(x)), ...)
  abline(h = 0, lty = 3)

  invisible(x)
}

# Validates a user-supplied wavelet table against the family/filter requested
# in the calling function and returns the payload passed to the C routines.
.match_wavelet_table <- function(wavelet.table, fam, filter.size, wavelet.filter){

  if(!inherits(wavelet.table, "wavelet_table"))
    stop("'wavelet.table' must be an object created by wtable().")

  if(wavelet.table$fam != fam)
    stop("'wavelet.table' was built for family '", wavelet.table$family,
         "', which differs from the requested family.")

  if(fam != 4 && wavelet.table$filter.size != filter.size)
    stop("'wavelet.table' was built for filter.size = ",
         wavelet.table$filter.size,
         ", which differs from the requested filter.size = ", filter.size, ".")

  if(fam == 4 && !isTRUE(all.equal(wavelet.table$wavelet.filter,
                                   wavelet.filter)))
    stop("'wavelet.table' was built with a different wavelet.filter.")

  list(wavelet.table$phi, wavelet.table$psi)
}
