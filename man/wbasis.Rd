\name{wbasis}
\alias{wbasis}
\title{
Decomposed Wavelet Basis
}
\description{
This function calculates, for each observation of the data, a matrix of
wavelet basis using the Daubechies-Lagarias algorithm. The first few columns
correspond to the scaling functions \eqn{\phi_{Jk}}, followed by the wavelets at
the coarsest level up to the wavelets at the finest resolution level, i.e.,
\eqn{\psi_{jk}, j_0 \le j \le J-1}.
}
\usage{
wbasis(x, j0 = 0, J, family = "Daublets", filter.size = 20,
       prec.wavelet = 30, periodic = TRUE, wavelet.filter)
}
\arguments{
\item{x}{Vector containing the data. Its length does not need to be power of 2.}
\item{j0}{The coarsest level at which the decomposed wavelet basis must be
calculated. Default is \code{j0 = 0}.}
\item{J}{The resolution level at which the wavelet basis is decomposed to be
calculated. It is associated to the multiresolution space of the calculated
basis.}
\item{family}{The family of wavelets to use. It can be \emph{Daublets},
\emph{Symmlets} or \emph{Coiflets}. It can also be \emph{Own}, if a filter is
provided. See Details for more information.}
\item{filter.size}{The size of the filter. The available values depend on the
chosen family. See Details for more information.}
\item{prec.wavelet}{The number of iterations to be performed in the Daubechies-
Lagarias algorithm, which is used to evaluate the scaling functions of the
specified wavelet basis at the data points.}
\item{periodic}{If it is TRUE (default), the periodic wavelet basis will be
used. We recommend to use \code{periodic = TRUE}, because it can handle boundary
conditions.}
\item{wavelet.filter}{Use this to provide your own filter. To use this argument,
you must specify \code{family = "Own"}. Do not use it, if you are not sure about
what you are doing.}
}
\details{
The scaling function \eqn{\phi} and the wavelet \eqn{\psi} are obtained
according to a wavelet filter with specific size, which may be provided by the
specification of \code{family} and \code{filter.size}, or by
\code{wavelet.filter}. Such a function has an associated multiresolution
analysis
\deqn{\ldots \subset V_{J-1} \subset V_J \subset V_{J+1} \subset \ldots.}
The wavelet \eqn{\psi} generates an orthonormal basis
\eqn{\{\psi_{jk}(\cdot)\}_k} of the orthogonal complement \eqn{W_j} of the
multiresolution space \eqn{V_j} in \eqn{V_{j+1}}, i.e.,
\deqn{W_{j} = V_{j+1} \ominus V_{j},} for any integer \eqn{j}.
Therefore, one can see that
\deqn{V_J = V_{j_0} \oplus W_{j_0} \oplus \ldots \oplus W_{J-1},}
where \eqn{j_0} corresponds to the coarsest level and \eqn{J-1} is known as the
finest level of the wavelet decomposition. Thus, one have that
\deqn{\{\{\phi_{j_0k}(\cdot); k\} \cup \{\psi_{jk}(\cdot); j_0 \le j \le J-1, k\}\}}
compose an orthonormal basis for \eqn{V_J}. The basis above is provided by the
present function.

Each line \eqn{i} of the resulting matrix is related to each element \eqn{x_i}
of \eqn{x}. The columns of the matrix are organized as
\deqn{\phi_{j_0,k_{\min}}(\cdot) \ldots \phi_{j_0,k_{\max}}(\cdot)
\psi_{j_0,k_{\min}}(\cdot) \ldots \psi_{J-1,k_{\max}}(\cdot).}
If \code{periodic = TRUE}, then \eqn{k_{\min} = 0} and
\eqn{k_{\max} = 2^(j-1)-1}, for each level \eqn{j}. Otherwise, \eqn{k_{\min}}
and \eqn{k_{\max}} will be all the integers where their corresponding function,
\eqn{\phi_{j_0k}} or \eqn{\psi_{jk}, j_0 \le j \le J - 1}, is non-null for at
least one observation of \eqn{x}.

The scaling functions \eqn{\phi_{j_0k}} and wavelets \eqn{\psi_{jk}} are
internally evaluated at the data points efficiently, using the
Daubechies-Lagarias algorithm (Daubechies and Lagarias, 1992).% by the commands
%\command{\link{PHI}} and \command{\link{PSI}}, respectively.

With respect to the levels \eqn{j_0} and \eqn{J}, it is important to mention
that the function only works for integer values satisfying
\eqn{0 \le j_0 \le J}. In the case where \eqn{j_0 = J}, the function will assume
the decomposition of the multiresolution space as \eqn{V_J = V_{j_0}}.
Therefore, it will take into account an orthonormal basis only with scaling
functions \eqn{\{\phi_{j_0k}(\cdot)\}_k}.

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

This function accepects wavelet filters provided by the user. In this case, the
filter is given in the argument \code{wavelet filter}, while setting
\code{family = "Own"}. If the filter is provided and \code{family} is not
\code{"Other"}, a warning is released and the function will still consider the
wavelet filter given. \bold{Caution: the function does not check if the filter
provided is in fact a wavelet filter. Therefore, do not use
\code{wavelet.filter}, unless you know what you are doing.}

This function works only with real valued observations. If, for some reason, the
data set is composed by complex numbers, the function will be applied only to
the real part of the data and a warning message will be released. Furthermore,
eventually, some observations of the data may not be finite (\code{NA},
\code{NaN}, \code{Inf} or \code{-Inf}). If this happens, the associated line to
this observation in the matrix of the scaling function basis will be composed by
\code{NA} values, and a warning message will be released.
}
\value{
A matrix in which each line \eqn{i} corresponds to
\deqn{\{\{\phi_{j_0k}(x_i); k\} \cup \{\psi_{jk}(x_i); j_0 \le j \le J-1, k\}\},}
for different values of \eqn{k}.
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
}
\seealso{\command{\link{PHI}}, \command{\link{PSI}}}
\examples{
#
# Wavelet estimation of a non-linear regression model
# Y = f(X) + e
#
set.seed(123)
n <- 100
x <- sort(runif(n))
f <- function(x) sin(2*pi*x)
e <- rnorm(n, mean = 0, sd = 0.25)
y <- f(x) + e

w <- wbasis(x, j0 = 0, J = 3, family = "Daublets", filter.size = 18,
            prec.wavelet = 30)

# Estimating the wavelet coefficients
mod <- lm(y~w-1)
beta <- coef(mod) # estimates

# Calculating estimates of f(x)
new.obs <- 0:(n-1)/n
myw <- wbasis(new.obs, j0 = 0, J = 3, family = "Daublets", filter.size = 18,
             prec.wavelet = 30)
f.est <- drop(myw \%*\% beta) # estimates of f(new.obs)

# Let's see the result
plot(x, y, main = "Regression model")
plot(f, 0, 1, lwd = 2, add = TRUE)
points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
legend("topright", legend = c("Real function", "Estimate"), col = 1:2, lty = 1)
#
#
# Next example
# ------------
# 
# As above, but now with a different value of J and
# regularizing the wavelet coefficients by the hard threshold
# with lambda = 0.1
# Y = f(X) + e
#
set.seed(123)
n <- 100
x <- sort(runif(n))
f <- function(x) sin(2*pi*x)
e <- rnorm(n, mean = 0, sd = 0.25)
y <- f(x) + e

w <- wbasis(x, j0 = 0, J = 5, family = "Daublets", filter.size = 18,
            prec.wavelet = 30)

# Estimating the wavelet coefficients
mod <- lm(y~w-1)
beta <- coef(mod) # estimates

# Calculating estimates of f(x)
new.obs <- 0:(n-1)/n
myw <- wbasis(new.obs, j0 = 0, J = 5, family = "Daublets", filter.size = 18,
             prec.wavelet = 30)
f.est <- drop(myw \%*\% beta) # estimates of f(new.obs)

# Regularizing beta
lambda <- 0.1 # threshold
beta.thr <- beta*(abs(beta)>lambda)
f.est.thr <- drop(myw \%*\% beta.thr) # regularized estimates of f(new.obs)

# Let's see the result
plot(x, y, main = "Regression model")
plot(f, 0, 1, lwd = 2, add = TRUE)
points(new.obs, f.est, col = 2, lwd = 2, type = 'l')
points(new.obs, f.est.thr, col = 3, lwd = 2, type = 'l')
legend("topright", legend = c("Real function", "Raw Estimate", "Reg. Estimate"),
       col = 1:3, lty = 1)
}
\keyword{smooth}
