\name{wavedec}
\alias{wavedec}
\title{
Wavelet decomposition
}
\description{
Decomposes a signal using the Mallat's pyramidal algorithm (Mallat, 1989).
}
\usage{
wavedec(x, j0 = 0, family = "Daublets", filter.size = 20, wavelet.filter)
}
\arguments{
\item{x}{Vector containing the data to be decomposed. Its length must be
be power of 2.}
\item{j0}{The coarsest resolution level where the signal is to be
decomposed.}
\item{family}{The family of wavelets to use. It can be \emph{Daublets},
\emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter is
provided \code{wavelet.filter}.}
\item{filter.size}{The size of the filter. The available values depend on
the chosen family. See Details for more information.}
\item{wavelet.filter}{Use this to provide your own filter. To use this argument,
you must specify \code{family = "Own"}. Do not use it, if you are not sure about
what you are doing.}
}
\details{ 
This function applies the discrete wavelet transform to the data
\code{x}, whose length must be a power of 2. As result, the function returns
\deqn{c_{j_0,0} \ldots c_{j_0,2^{j_0}-1} d_{j_0,0} \ldots d_{j_0,2^{j_0}-1}
\ldots d_{J-1,0} \ldots d_{J-1,2^{J-1}-1},}
where \eqn{c_{j_0,k}} represent the scale coefficients in the coarsest
level and \eqn{d_{j,k}} represent the detail coefficients in level \eqn{j}.

For the moment, the function deal with boundary conditions by applying
periodic extensions to the data. In the future, a reflective extesion may
be taken into account.

For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
\code{prec.wavelet} and \code{wavelet.filter}, see the help of
\command{\link{wbasis}}.

}
\value{
A vector of wavelet coefficients decomposed from \code{x}, where the
coarsest resolution level is \code{j0}.
}
\author{Michel H. Montoril \email{michel@ufscar.br}}
\references{
Mallat, S.G. (1989). A theory for multiresolution signal decomposition:
the wavelet representation. \emph{IEEE Trans. Patt. Anal. and Mach.
Intell.}, 11, 674-693, \url{https://ieeexplore.ieee.org/document/192463}.
}
\seealso{\code{\link{waverec}}}
\examples{
#
# A simple example
#
set.seed(123)
n <- 64
x <- sort(runif(n))


# Let's apply the wavelet decomposition to x using a 18-tap
# filter of a Daulet.
# Consider j0 = 3 as the coarsest resolution level:
wdx <- wavedec(x, j0 = 3, family = "daublets", filter.size = 18)
}
\keyword{smooth}
