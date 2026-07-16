# WaveBased

<!-- badges: start -->
[![R-CMD-check](https://github.com/michelcias/WaveBased/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/michelcias/WaveBased/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Wavelet-Based Estimation** — an R package with tools to tackle different
nonparametric problems through compactly supported orthonormal wavelet bases.

Since compactly supported orthonormal wavelets (other than Haar) have no closed
form, `WaveBased` evaluates the scaling and wavelet functions at arbitrary
points using the **Daubechies–Lagarias algorithm** (Daubechies and Lagarias,
1992). On top of that it provides the discrete wavelet transform and a
wavelet-based density estimator for size-biased data.

The package ships with families **Daublets**, **Symmlets** and **Coiflets**
(filters from PyWavelets, Lee et al., 2019), and also accepts user-supplied
filters.

## Installation

```r
# install.packages("remotes")
remotes::install_github("michelcias/WaveBased")
```

## Main functions

| Function     | Purpose                                                                  |
|--------------|--------------------------------------------------------------------------|
| `PHI()`      | Matrix of scaling functions evaluated at the data                        |
| `PSI()`      | Matrix of (mother) wavelets evaluated at the data                        |
| `wbasis()`   | Full decomposed wavelet basis (scaling functions + wavelets)             |
| `wavedec()`  | Discrete wavelet decomposition (Mallat's pyramidal algorithm)            |
| `waverec()`  | Inverse wavelet reconstruction                                           |
| `wdensity()` | Wavelet-based density estimation for univariate size-biased data         |
| `wtable()`   | Precomputed interpolation tables for fast basis evaluation               |

The package also includes the `bac` dataset (blood alcohol concentrations of
drivers in fatal US accidents, 2019), used in Montoril et al. (2021).

## Quick start

Decomposition/reconstruction (perfect reconstruction):

```r
library(WaveBased)

set.seed(123)
x   <- sort(runif(64))                                   # length must be a power of 2
wdx <- wavedec(x, j0 = 3, family = "daublets", filter.size = 18)
wrx <- waverec(wdx, j0 = 3, family = "daublets", filter.size = 18)

sum(abs(wrx - x))   # ~ 0
```

Wavelet regression with `wbasis()`:

```r
set.seed(123)
n <- 100
x <- sort(runif(n))
y <- sin(2 * pi * x) + rnorm(n, sd = 0.25)

w   <- wbasis(x, j0 = 0, J = 3, family = "Daublets", filter.size = 18)
fit <- lm(y ~ w - 1)
```

Density estimation for size-biased data:

```r
data(bac)
est <- wdensity(data = bac, wf = function(x) 0.1 + 0.9 * x,
                power.dens = 0.5, J1 = ceiling(0.95 * log2(length(bac))),
                family = "s", filter.size = 20, warped = TRUE)
```

### Fast evaluation with precomputed tables

The Daubechies–Lagarias evaluation is exact but costly for large samples.
`wtable()` tabulates the scaling function and the wavelet once (a ~2.5 MB
object, independent of the resolution level `J`), and the basis functions are
then evaluated by table lookup with linear interpolation — typically **30 to
500 times faster**, with errors below `1e-5` for filters of size 8 or larger:

```r
tab <- wtable(family = "symmlets", filter.size = 20)   # once per session
tab   # prints the measured interpolation error

x <- sort(runif(1e5))
w <- PHI(x, J = 9, family = "symmlets", filter.size = 20,
         wavelet.table = tab)                          # ~30x faster
```

The same `wavelet.table` argument is accepted by `PSI()`, `wbasis()` and
`wdensity()`, and `wtable()` also accepts user-supplied filters via
`family = "Own"` with `wavelet.filter`.

With `cache = TRUE` the table is additionally saved to the user cache
directory, so it is generated only once *across R sessions* — later calls
load it from disk in milliseconds:

```r
tab <- wtable(family = "symmlets", filter.size = 20, cache = TRUE)

wtable_cache()             # inspect the cached tables
wtable_cache(clear = TRUE) # remove them
```

## Citation

If you use `WaveBased` in your research, please cite the package. You can get an
up-to-date citation (including BibTeX) from R with:

```r
citation("WaveBased")
```

For the size-biased density methodology, please also cite:

> Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2021). Wavelet-based
> estimation of power densities of size-biased data. *arXiv:2112.12895*.

## References

- Daubechies, I. and Lagarias, J. C. (1992). Two-Scale Difference Equations II.
  Local Regularity, Infinite Products of Matrices and Fractals.
  *SIAM Journal on Mathematical Analysis*, 24(4), 1031–1079.
- Mallat, S. G. (1989). A theory for multiresolution signal decomposition: the
  wavelet representation. *IEEE Trans. Pattern Anal. Mach. Intell.*, 11, 674–693.
- Lee, G., Gommers, R., Waselewski, F., Wohlfahrt, K. and O'Leary, A. (2019).
  PyWavelets: A Python package for wavelet analysis. *Journal of Open Source
  Software*, 4(36), 1237.
- Montoril, M. H., Pinheiro, A. and Vidakovic, B. (2021). Wavelet-based
  estimation of power densities of size-biased data. *arXiv:2112.12895*.

## License

GPL (>= 2). Maintainer: Michel H. Montoril.
