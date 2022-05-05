\name{bac}
\docType{data}
\alias{bac}
\title{Blood alcohol concentration}
\description{
Blood alcohol concentration of drivers involved in fatal accidents that
occurred in the USA, during the year of 2019.
}
\usage{data(bac)}
\source{
This dataset contains 2,495 blood alcohol concentrations (BAC) of drivers
involved in fatal accidents that occurred in the USA, during the year of
2019. The data was collected from the National Highway Traffic Safety
Administration Department of Transportation (\url{www.nhtsa.dot.gov}). It
is part of The Fatality Analysis Reporting System (FARS).

The FARS contains data of traffic crashes within the 50 States, the
District of Columbia, and Puerto Rico. To be included in FARS, a crash
must involve a motor vehicle traveling on a traffic way customarily open
to the public, and must result in the death of a vehicle occupant or a
nonoccupant within 30 days of the crash. A better description of the data
is available in this
\href{https://www.nhtsa.gov/crash-data-systems/fatality-analysis-reporting-system}{link}.


The \command{bac} data are expressed in grams/100 ml, with three decimal
places (or significant digits). The 2019 FARS/CRSS Coding and Validation Manual
(available \href{https://static.nhtsa.gov/nhtsa/downloads/FARS/FARS-DOC/Coding\%20and\%20Validation\%20Manual/2019\%20FARS\%20CRSS\%20Coding\%20and\%20Validation\%20Manual\%20-\%20DOT\%20HS\%20813\%20010.pdf}{here})
contains details of the variables collected. The \command{bac} data
available in this package consider only fatal accidents where alcohol is
involved (according to the police report); and vehicles classified as
automobiles, automobiles derivatives, utility vehicles and two-wheel
motorcycles. Crashes that are not included in the state highway inventory,
not reported or unknown were discarded.

The \command{bac} data are used as an application in
Montoril et al. (2021). It is reproduced as an example in
\command{\link{wdensity}}.
}
\author{Michel H. Montoril \email{michel@ufscar.br}}
\references{
Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2021). Wavelet-based
estimation of power densities of size-biased data. \emph{arXiv:2112.12895
[math, stat]}. \url{https://arxiv.org/abs/2112.12895}
}
\examples{
data(bac)
hist(bac)
}
\keyword{datasets}
