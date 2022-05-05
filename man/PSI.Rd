\name{PSI}
\alias{PSI}
\title{
Mother Wavelet Basis
}
\description{
Calculates a matrix of wavelet functions \eqn{\psi_{Jk}}, where the lines are
related to the elements of the data set \eqn{x} and the columns are related to
the values of \eqn{k} for which \eqn{\psi_{Jk}} is non-null.
}
\usage{
PSI(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
    periodic = TRUE, wavelet.filter)
}
\arguments{
\item{x}{Vector containing the data. Its length does not need to be power of 2.}
\item{J}{The resolution level at which the (mother) wavelet basis must be
calculated.}
\item{family}{The family of wavelets to use. It can be \emph{Daublets},
\emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter is
provided. See Details for more information.}
\item{filter.size}{The size of the filter. The available values depend on the
chosen family. See Details for more information.}
\item{prec.wavelet}{The number of iterations to be performed in the Daubechies-
Lagarias algorithm, which is used to evaluate the scaling functions of the
specified wavelet basis at the data points.}
\item{periodic}{If it is TRUE (default), the periodic (mother) wavelet basis
will be used. We recommend to use \code{periodic = TRUE}, because it can handle
boundary conditions.}
\item{wavelet.filter}{Use this to provide your own filter. To use this argument,
you must specify \code{family = "Own"}. Do not use it, if you are not sure about
what you are doing.}
}
\details{
The wavelet function \eqn{\psi} is obtained according to a wavelet filter with
specific size, which may be provided by the specification of \code{family} and
\code{filter.size}, or by \code{wavelet.filter}. Such a function generates an
orthonormal basis \eqn{\{\psi_{Jk}(\cdot)\}_k} of the orthogonal complement
\eqn{W_J} of the multiresolution space \eqn{V_J} in \eqn{V_{J+1}},
i.e., \deqn{W_{J} = V_{J+1} \ominus V_{J}} of an associated multiresolution
analysis. This function calculates a matrix of values \eqn{\psi_{Jk}(x)}, where
the \eqn{i}-th line is related to the \eqn{i}-th element of the data set \eqn{x}
and the columns are related to the values of \eqn{k} for which \eqn{\psi_{Jk}}
is non-null. The (mother) wavelet basis is defined as
\deqn{\psi_{Jk}(x) = 2^{J/2} \psi(2^J x - k).}

If \code{periodic = FALSE}, the first column corresponds to the minimum \eqn{k}
where \eqn{\psi_{Jk}} is non-null for at least one element of \eqn{x}.
Analogously, the last column is related to the maximum reasonable value of
\eqn{k}. If \code{periodic = TRUE}, the columns will correspond to the values of
\eqn{k} from \eqn{0} to \eqn{2^J - 1}.

The wavelets are calculated based on the dilation equation
(see, e.g., Vidakovic, 1999, Theorem 3.5.5). The scaling function used is
evaluated at the data points efficiently by the Daubechies-Lagarias algorithm
(Daubechies and Lagarias, 1992). Coded kindly by Brani Vidakovic. More details
about scaling functions in \command{\link{PHI}}.

The families available here are \emph{Daublets}, \emph{Symmlets} and
\emph{Coiflets}. The first and the second correspond to the families of
Daubechies' extremal phase and Daubechies' least asymmetric wavelets,
respectively. The argument \code{family} here is not case sensitive, and it
accepts the first letter of the available families. For example, it is allowed
to use \code{family = "d"} for the \emph{Daublet} family. The filters used here
were taken from the Python package \pkg{PyWavelets} (Lee et al., 2019), which
provides extended precision for the families \emph{Daublets} and
\emph{Coiflets}. Furhtermore, \pkg{PyWavelets} provides a wider range of filters
for the three families.

The filter sizes available depend on the chosen family:
\describe{
\item{Daublets}{\code{filter.size = 2} \code{(Haar), 4, 6, ..., 74, 76}.}
\item{Symmlets}{\code{filter.size = 4, 6, 8, ..., 38, 40}.}
\item{Coiflets}{\code{filter.size = 6, 12, 18, ..., 96, 102}.}
}

The argument \code{filter.size} is associated with the number of vanishing
moments (\eqn{N}) of each wavelet family. For \emph{Daublets} and
\emph{Symmlets}, one has \eqn{\int x^n \psi(x) dx = 0, n = 0, \ldots, N - 1},
satisfying the relation \code{N = filter.size/2}. In the case of
\emph{Coiflets}, both scaling function \eqn{\phi} and wavelet \eqn{\psi} have
\eqn{N} vanishing moments satisfying the relation \code{N = filter.size/6}.

This function works only with real valued observations. If, for some reason, the
data set is composed by complex numbers, the function will be applied only to
the real part of the data and a warning message will be released. Furthermore,
eventually, some observations of the data may not be finite (\code{NA},
\code{NaN}, \code{Inf} or \code{-Inf}). If this happens, the associated line to
this observation in the matrix of the scaling function basis will be composed by
\code{NA} values, and a warning message will be released.
}
\value{
A matrix in which each line \eqn{i} corresponds to \eqn{\psi_{Jk}(x_i)}, for
different values of \eqn{k} (associated to the columns of the matrix).
}
\author{Michel H. Montoril \email{michel@ufscar.br}}
\references{
Daubechies, I. and Lagarias, J.C. (1992). Two-Scale Difference Equations II.
Local Regularity, Infinite Products of Matrices and Fractals. \emph{SIAM Journal
on Mathematical Analysis}, 24(4), 1031--1079,
\url{https://epubs.siam.org/doi/10.1137/0523059}.

Lee, G., Gommers, R., Waselewski, F., Wohlfahrt, K. and O'Leary, A. (2019).
PyWavelets: A Python package for wavelet analysis. \emph{Journal of Open Source
Software}, 4(36), 1237, \url{https://doi.org/10.21105/joss.01237}.

Vidakovic, B. (1999). \emph{Statistical Modeling by Wavelets}. John Wiley, New
York.
}
\seealso{\command{\link{PHI}}, \command{\link{wbasis}}}
\examples{
#
# Calculating the non-periodized and the periodized versions
# of psi
#
set.seed(123)
n <- 10 # setting the size of the vector
x <- sort(runif(n)) # generating some specific vector of size n

psi <- PSI(x, J = 3, family = "Daublets", filter.size = 18, prec.wavelet = 30,
           periodic = FALSE)

psi.per <- PSI(x, J = 3, family = "Daublets", filter.size = 18,
               prec.wavelet = 30, periodic = TRUE)
}
\keyword{smooth}
