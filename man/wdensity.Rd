\name{wdensity}
\alias{wdensity}
\title{
Wavelet-based density estimation for univariate size-biased data
}
\description{
The \code{wdensity} function calculates wavelet-based density estimates for
univariate size-biased data. The method is based on the estimation of the
power of the density of interest. It uses the Daubechies-Lagarias algorithm
to calculate the periodized version of a compactly supported orthonormal
wavelet basis.
}
\usage{
wdensity(data, from, to, length.obs = 250, wf = NULL, power.dens, J0 = 0, J1,
         family = "Daublets", filter.size = 10, prec.wavelet = 30,
         wavelet.filter, dens.biased, warped = TRUE, warp.fun, dwarp.fun,
         rescale = TRUE, eps = 1.9^(-J1), thresh = "hard", plot = TRUE, main,
         xlab, ylab, type, ...)
}
\arguments{
\item{data}{Vector containing the data from which the density estimate is to be
computed. The length of this vector does not have to be a power of 2.}
\item{from, to}{Minimum and maximum values of the grid of points at which the
density is to be estimated. The default is \code{from = min(data)} and
\code{to = max(data)}.}
\item{length.obs}{Number of equally spaced points at which the density is to be
estimated. The default is \code{length.obs = 250}.}
\item{wf}{Weighting function. It must be a function.}
\item{power.dens}{Power of the density to be estimated. It must be a value
greater than or equal to 0.5. The Default is \code{power.dens = 0.5}.}
\item{J0}{The coarsest resolution level for the wavelet basis. Default is
\code{J0 = 0}.}
\item{J1}{The resolution level at which the wavelet basis is decomposed to be
calculated. It is associated to the multiresolution space where the function is
projected.}
\item{family}{The family of wavelets to use. It can be \emph{Daublets},
\emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter is
provided \code{wavelet.filter}.}
\item{filter.size}{The size of the wavelet filter. The available values depend
on the chosen wavelet family.}
\item{prec.wavelet}{The number of iterations to be performed in the Daubechies-
Lagarias algorithm.}
\item{wavelet.filter}{Use this to provide your own filter. To use this argument,
you must specify \code{family = "Own"}. Do not use it, if you are not sure about
what you are doing.}
\item{dens.biased}{Density estimates of the size-biased data. These estimates of
the (ordered) \code{data} must be provided. By default, a kernel density
estimate is calculated.}
\item{warped}{Logical. It indicates if the wavelet basis must be warped.}
\item{warp.fun}{A (cumulative distribution) function. By default, it
provides the empirical distribution function of \code{data}. Only used if
\code{warped = TRUE}.}
\item{dwarp.fun}{A function. Derivative of the \code{warp.fun}. By default, it
provides a kernel density estimate of \code{data}. Only used if
\code{warped = TRUE}.}
\item{rescale}{Logical. If TRUE, the \code{data} will be rescaled in the
interval [\code{eps}, 1-\code{eps}], where \code{eps} is another argument. Only
used if \code{warped = FALSE}.}
\item{eps}{A value in the unit interval. It is used to rescale the \code{data}
in the interval \eqn{[}\code{eps}\eqn{, 1-}\code{eps}\eqn{]}. By default,
\code{eps} \eqn{= 1.9^{-J1}}. Only used if \code{warped = FALSE}.}
\item{thresh}{Indicates if the wavelet-density estimate should be regularized
and how. If \code{thresh = "linear"}, the density estimate will not be
regularized; if \code{thresh = "hard"}, the hard-thresholding regularizer is
applied to the detail coefficient estimates; and if \code{thresh = "soft"}, the
soft-thresholding regularizer is applied to the detail coefficient estimates.}
\item{plot}{Logical. If \code{plot = TRUE}, the density estimate is to be
plotted.}
\item{main, xlab, ylab, type}{Plotting parameters with useful defaults.}
\item{...}{Further plotting parameters.}
}
\details{
This function provides density estimates of size-biased data according to
Montoril et al. (2021). It estimates the power of the density of interest by
wavelet basis, i.e., \eqn{\hat{f^a}, a \ge 0.5}. After that, the density is
finally estimated by \eqn{[\hat{f^a}]^{1/a}}. Periodized version of ordinary
orthonormal wavelet bases of families \emph{Daublets}, \emph{Symmlets} and
\emph{Coiflets} are allowed, as well as warped bases of these families. It also
accepts wavelet filters provided in \code{wavelet.filter} (see the help of
\command{\link{wbasis}}).

The elements of the \code{data} must all be finite numbers. Otherwise, they
will be removed before the density estimation procedure. In this case, a
warning message will be released.

For stability reasons, the values at which the density estimates are to be
computed must satisfy \eqn{\min(data) \le from \le to \le \max(data)}.

By default, the argument \code{dens.biased} provides the kernel density
estimates of the \code{data}. It is computed by
\command{\link[stats]{density.default}}, where the kernel adopted is the
Gaussian and the bandwidth selection used is the "\code{SJ}" proposed by
Sheather and Jones, (1991).

For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
\code{prec.wavelet} and \code{wavelet.filter}, see the help of
\command{\link{wbasis}}.

The arguments \code{warp.fun} and \code{dwarp.fun} are used only when
\code{warped = TRUE}. They must be functions representing a cumulative
distribution function and a probability density function, respectively. By
default, the function uses, respectively, the empirical distribution function
(EDF) and the kernel density estimate (KDE) of the \code{data}. The EDF is
computed by \command{\link[stats]{ecdf}}. The KDE is computed by
\command{\link[stats]{density.default}}, where the kernel adopted is the
Gaussian and the bandwidth selection used is the \code{"SJ"} proposed by
Sheather and Jones, (1991). Furthermore, when \code{warped = TRUE}, the
\code{data} are not rescaled, because they are scaled by the warping function
to be used.

When \code{warped = FALSE}, the \code{data} will be scaled when
\code{rescale = TRUE}. In this case, a linear transformation is applied and the
observations will belong to the interval \eqn{[eps, 1-eps]}. By default,
\eqn{eps = 1.9^{-J1}}, following the suggestion by Montoril et al. (2019). When
\code{rescale = FALSE}, the original \code{data} are to be used.

\bold{Caution: as the wavelet basis is periodized in the unit interval, if
\eqn{\max(data) > 1} or \eqn{\min(data) < 0}, DO NOT use \code{warped = FALSE}
with \code{rescale = FALSE}. Otherwise, the estimates will be unfeasible.}

}
\value{
A list containing the following components.
\item{y}{The grid of \code{length.obs} observations at which the densities are
estimated.}
\item{dens}{The estimated density values.}
}
\author{Michel H. Montoril \email{michel@ufscar.br}}
\references{
Montoril, M. H., Chang, W. and Vidakovic, B. (2019). Wavelet-Based Estimation
of Generalized Discriminant Functions. \emph{Sankhya B}, 81(2), 318-349,
\doi{10.1007/s13571-018-0158-1}.

Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2021). Wavelet-based
estimation of power densities of size-biased data. \emph{arXiv:2112.12895
[math, stat]}. \url{https://arxiv.org/abs/2112.12895}.

Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
selection method for kernel density estimation. \emph{Journal of the Royal
Statistical Society series B}, 53, 683-690.
\url{http://www.jstor.org/stable/2345597}.
}
\seealso{\command{\link{PHI}}, \command{\link{bac}}}
\examples{
# Density estimation of the Blood alcohol concentration of drivers
# involved in fatal accidents in the USA in 2019.
# This example reproduces the application in Montoril et al. (2021).

# Loading the data
data(bac)

# Let's consider four different possibilities of estimation:
# m1 - non-warped wavelets, non-rescaled data and power.density
#      equals to half
m1 = wdensity(data = bac, wf = function(x) 0.1+.9*x, power.dens = 0.5,
              J1 = ceiling(0.45*log2(length(bac))), family = "s",
              filter.size = 20, warped = FALSE, rescale = FALSE,
              thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
              col = 1, lty = 1, main = expression(m[1]))

# m2 - non-warped wavelets, non-rescaled data and power.density
#      equals to one
m2 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 1,
              J1 = ceiling(0.45*log2(length(bac))), family = "s",
              filter.size = 20, warped = FALSE, rescale = FALSE,
              thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
              col = 1, lty = 1, main = expression(m[2]))

# m3 - warped wavelets and power.density equals to half
m3 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 0.5,
              J1 = ceiling(0.95*log2(length(bac))), family = "s",
              filter.size = 20, warped = TRUE, rescale = FALSE,
              thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
              col = 1, lty = 1, main = expression(m[3]))

# m4 - warped wavelets and power.density equals to one
m4 = wdensity(data = bac, wf = function(x) 0.1+0.9*x, power.dens = 1,
              J1 = ceiling(0.95*log2(length(bac))), family = "s",
              filter.size = 20, warped = TRUE, rescale = FALSE,
              thresh = "hard", ylim = c(0.00000, 16.37781), lwd = 2,
              col = 1, lty = 1, main = expression(m[4]))
}
\keyword{smooth}
