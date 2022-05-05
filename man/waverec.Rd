\name{waverec}
\alias{waverec}
\title{
Wavelet reconstruction
}
\description{
Reconstructs a signal using the Mallat's pyramidal algorithm in the
reconstruction stage (Mallat, 1989).
}
\usage{
waverec(x, j0 = 0, family = "Daublets", filter.size = 20, wavelet.filter)
}
\arguments{
\item{x}{Vector containing the decomposed data. Its length must be
be power of 2.}
\item{j0}{The coarsest resolution level where the signal is decomposed.}
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
This function applies the inverse discrete wavelet transform to the (decomposed)
data \code{x}, whose length must be a power of 2. In this case, it is assumed
that the data is decomposed with the coarsest level \code{j0}.

For the moment, the function deal with boundary conditions by applying periodic
extensions to the data. In the future, a reflective extesion may be taken into
account.

For more details about \code{J0}, \code{J1}, \code{family}, \code{filter.size},
\code{prec.wavelet} and \code{wavelet.filter}, see the help of
\command{\link{wbasis}}.

}
\value{
A vector containing a reconstruction of the decomposed data \code{x}, whose
coarsest resolution level is \code{j0}.
}
\author{Michel H. Montoril \email{michel@ufscar.br}}
\references{
Mallat, S.G. (1989). A theory for multiresolution signal decomposition:
the wavelet representation. \emph{IEEE Trans. Patt. Anal. and Mach.
Intell.}, 11, 674-693, \url{https://ieeexplore.ieee.org/document/192463}.
}
\seealso{\code{\link{wavedec}}}
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

# Now, let's apply the wavelet reconstrution to wdx.
wrx <- waverec(wdx, j0 = 3, family = "daublets", filter.size = 18)

# The reconstructed vector should be equal to x (or, at least,
# as close as possible)
sum(abs(wrx-x))
}
\keyword{smooth}
