# WaveBased 2.10-1

## Performance

* The Daubechies-Lagarias algorithm, which evaluates the scaling function and
  the wavelet wherever no `wavelet.table` is given, applies its product of
  `prec.wavelet` matrices to a vector, by matrix-vector products, instead of
  forming the product: the cost per point falls from O(prec N^3) to
  O(prec N^2), for a filter of size N. Every function that evaluates the basis
  gains: the evaluation is about 15 times faster with Daublets of size 10, 44
  times with Symmlets of size 20 and 70 times with Coiflets of size 30. The
  values change by rounding only, below 1e-13.
* `PSI()`, `wbasis()` and the dense designs of `wall()` and `cv.wall()` no
  longer keep a workspace for every point until they return: `PSI()` at
  200,000 points with Daublets of size 20 reached 1.4 GB for a result of
  12 MB, and took 184 s, now 1 s.
* `wdensity()` no longer forms the matrix of the basis at the data and at the
  evaluation points, with (n + `length.obs`) x 2^`J1` entries, of which at
  most `filter.size - 1` per row are not null. The coefficients and the
  estimate are computed in C point by point, so the memory no longer grows
  with the product of the sample size and 2^`J1`, which made
  `warped = TRUE`, whose `J1` grows with the sample size, impractical for
  large samples: with 10,000 observations and Symmlets of size 20, a fit took
  12.6 s and 1.4 GB, and now takes 0.05 s; with 100,000 (`J1 = 16`), whose
  matrix alone would take 52 GB, it takes 0.8 s. The four `bac` examples of
  its help page fall from 2.4-2.8 s to about 0.013 s each.

## Minor changes

* `wdensity()` with `warped = TRUE` and the default warping no longer warns
  "collapsing to unique 'x' values" on tied data. The warping was already
  computed with `ties = mean`, which is now stated.

# WaveBased 2.10-0

## New features

* `wdensity()` gains the argument `correct`, `FALSE` by default. With
  `correct = TRUE`, the estimate is turned into a bona fide density by the
  correction of Glad, Hjort and Ushakov (2003): it becomes `max(0, f - xi)`,
  where the constant `xi`, also returned, makes it integrate to one over the
  domain of the estimator. That domain is the interval `[0, 1]` of the basis,
  in the scale of the data, with `warped = FALSE`, and the range of the data
  with `warped = TRUE`. It does not depend on `from`, `to` and `length.obs`,
  which only set where the estimate is evaluated. When `xi > 0`, the
  correction never increases the integrated squared error; when `xi < 0`,
  the missing mass is spread evenly over the domain. With
  `power.dens = 0.5` and `warped = FALSE`, the estimate is already a density
  and `xi = 0`. The integral is computed in C, on a dyadic grid of `[0, 1]`
  where the expansion is a convolution of the coefficients with the scaling
  function at the dyadic points, so no basis matrix is formed: in the warped
  `bac` example (`J1 = 11`, Symmlets of size 20), the correction takes
  0.04 s of a fit of 2.6 s.

## Bug fixes

* `wdensity()` returned `NaN` where the estimate of the power of the density
  was negative, with a `power.dens` other than 0.5 and 1: a negative number
  has no real power `1/power.dens`. The power is now taken of the positive
  part of the estimate, so these points are 0, as they already were with
  `power.dens = 1`.

# WaveBased 2.9-1

* Merges the check of `eps` in `wall()`, `cv.wall()` and `wdensity()`,
  released as 2.8-2, into 2.9-0. Both changes are kept: `eps` must lie in
  `[0, 0.5)`, and the argument `check.lawton` is available in these functions.

# WaveBased 2.9-0

## New features

* `PHI()`, `PSI()`, `wbasis()`, `wtable()`, `wdensity()`, `wall()`,
  `cv.wall()` and `wmixreg()` gain the argument `check.lawton`, `TRUE` by
  default. With `check.lawton = FALSE`, a filter given in `wavelet.filter`
  is not checked against the condition of Lawton (1991), nor is the warning
  of a filter that nearly fails it released, which serves the study of
  degenerate filters. Their scaling function is evaluated all the same: that
  of `c(1, 0, 0, 1)/sqrt(2)` is 1/3 on `[0, 3)`, the exact value, but its
  integer translations are not orthonormal. The even length, the sum of
  `sqrt(2)` and the orthogonality to the translations by even shifts are
  still checked.

## Bug fixes

* `wavedec()`, `waverec()`, `bayesthresh()`, `bwmixreg()` and `bwregime()`
  no longer check the condition of Lawton, which 2.8-0 required of them. It
  concerns the scaling function, which these functions do not evaluate: the
  discrete transform is orthogonal, and reconstructs the signal exactly,
  whenever the filter sums to `sqrt(2)` and is orthogonal to its
  translations by even shifts. A degenerate filter as
  `c(1, 0, 0, 1)/sqrt(2)` stopped them, and is now accepted.

# WaveBased 2.8-2

## Bug fixes

* `wall()`, `cv.wall()` and `wdensity()` stop when the rescaling constant
  `eps` lies outside `[0, 0.5)`, whatever its source: supplied by the user,
  or given by the default `1.9^(-J)` (`1.9^(-J1)` in `wdensity()`), as at
  `J = 1`. The value is checked once it is resolved, by a single rule shared
  by the three functions, and the message reports it. At `eps = 0.5` the
  interval `[eps, 1 - eps]` is a single point, and above it the covariates
  were mapped, reversed, into a narrow band around 1/2. In `wall()`, the fit
  ran on that band, but `predict()` truncates the rescaled values to
  `[eps, 1 - eps]`, which then sends every observation to the same point, so
  the predictions were constant. `cv.wall()` scores its folds on the design
  of the fit, without that truncation, so it could report a low
  cross-validated error and select a `J` whose fitted object only predicts a
  single class. Before, the default `eps >= 0.5` only raised a warning (since
  2.6-1), and `wdensity()` accepted any `eps` in `[0, 1]`.
* `wdensity()` checks `eps` only when it is used, i.e., with
  `warped = FALSE` and `rescale = TRUE`.

# WaveBased 2.8-1

## Bug fixes

* `wall()` and `cv.wall()` now warn when `wavelet.filter` is given alongside a
  built-in family, as the other functions do: the filter is the one used, and
  `family` and `filter.size` are ignored. Before, they used the filter
  silently. `cv.wall()` warns once, and not at each fold or candidate `J`.

# WaveBased 2.8-0

## New features

* A user-supplied filter (`family = "Own"`, or a `wavelet.filter` given
  alongside a built-in family) is now checked to be the low-pass filter of an
  orthonormal wavelet basis of compact support, in `PHI()`, `PSI()`,
  `wbasis()`, `wavedec()`, `waverec()`, `wtable()`, `wdensity()`, `wall()`,
  `cv.wall()`, `wmixreg()`, `bayesthresh()`, `bwmixreg()` and `bwregime()`.
  Before, no function checked it, and the documentation said so. The checks
  are the conditions that hold for any such filter: an even length, a sum of
  `sqrt(2)`, the orthogonality to the translations by even shifts, and the
  condition of Lawton (1991), without which the integer translations of the
  scaling function are not orthonormal (`c(1, 0, 0, 1)/sqrt(2)` satisfies the
  other three but not this one). The number of vanishing moments is not
  checked.
* The function stops when the sum or the orthogonality fail by more than
  `1e-6`, with a hint when the filter is the high-pass one, sums to 1 or sums
  to `-sqrt(2)`, and when the condition of Lawton fails. It warns when they
  hold up to an error between `1e-10` and `1e-6` only, as with coefficients
  copied with a few digits, and when the filter nearly fails the condition of
  Lawton, which makes the scaling function very irregular. The filters of
  the built-in families pass silently. `wdensity()`, `wall()`, `cv.wall()`
  and `wmixreg()` warn once, and not at each of their internal calls.
* This changes the behavior: a filter accepted before may now stop. The odd
  length is still refused with a message that says "even", now in all the
  functions above and not only in `wtable()`.

## Documentation

* The parameter `wavelet.filter` and the Details of `PHI()` and `wbasis()`
  describe the checks, in place of the warning that the filter was not
  checked.

# WaveBased 2.7-5

## Bug fixes

* The error for an unavailable `filter.size` of the Coiflets listed the sizes
  from 12, but 6 is available too.

# WaveBased 2.7-4

## Documentation

* The vanishing moments of the Coiflets were wrong in `PHI()`, `PSI()` and
  `wbasis()`: they gave `filter.size/6` for both the scaling function and the
  wavelet. The wavelet has `N = filter.size/3` vanishing moments, and the
  moments of the scaling function of orders 1 to `N` vanish about the point
  `N` of its support `[0, filter.size - 1]`. The paragraph on user-supplied
  filters of `PHI()` and `wbasis()` named the family `"Other"` where the code
  uses `"Own"`.

# WaveBased 2.7-3

## Bug fixes

* `wdensity()` now removes the non-finite elements of `data`, with a
  warning, as it meant to. The test was inverted: it removed them only when
  no element was finite, so a single `NA` made the function stop with
  "missing value where TRUE/FALSE needed". It now stops when no element is
  finite. The default `xlab` of the plot is still the name of `data`.

# WaveBased 2.7-2

## Bug fixes

* `wdensity()` now uses the filter given in `wavelet.filter` in the
  thresholding as well. Before, the projection used it but the transforms of
  `thresh = "hard"` and `thresh = "soft"` used `family` and `filter.size`:
  with `family = "Own"` the default `thresh = "hard"` stopped, and with a
  filter given alongside a built-in family the coefficients were thresholded
  in the basis of that family. `thresh = "linear"` was not affected.

# WaveBased 2.7-1

## Bug fixes

* `boundary = "interval"` now stops when the filter is user-supplied
  (`family = "Own"`, or a `wavelet.filter` given alongside a built-in family),
  in `PHI()`, `PSI()`, `wbasis()`, `wavedec()` and `waverec()`, and through
  them in `wdensity()`, `wall()` and `cv.wall()`. The documentation already
  required a Daublet or Symmlet filter, since the Cohen-Daubechies-Vial
  construction needs the minimal-length filters of Daubechies, but the code
  accepted any filter and derived the boundary blocks from it at run time.
  With the filters of fewer than `L/2` vanishing moments that we tried (a
  Coiflet, and random orthonormal filters of lengths 4 to 12), the
  construction stopped with a message that blamed the length of the filter.
  The message now names the family. The paragraph of `wbasis()` that
  described the run-time construction for `family = "Own"` was corrected.

# WaveBased 2.7-0

## New features

* `wavedec()` and `waverec()` transform a matrix by rows: each row is a
  signal, and the result is a matrix of the same dimensions, with the same row
  names. This replaces the error of 2.6-4. The rows go through a single call
  to the compiled code, which sets up the filter, and the boundary blocks of
  `boundary = "interval"`, only once. The numbers are the same as those of
  `t(apply(x, 1, wavedec, ...))`, from 2 to 13 times faster in our
  measurements (500 to 5000 rows of 256 or 1024 points). The number of
  columns must be a power of two. A matrix with a single column is still taken
  as a vector, and an array with more than one dimension of extent greater
  than one, or a data frame, still stops. The two-dimensional transform of an
  image is not computed by these functions.

# WaveBased 2.6-4

## Bug fixes

* `wavedec()` and `waverec()` now stop when `x` is a matrix, or an array, with
  more than one dimension of extent greater than one. Before, it was flattened
  by columns and transformed as a single vector, which is neither the
  transform of its rows nor that of its columns, and the result was returned
  as a vector without an error. The message points to
  `t(apply(x, 1, wavedec, ...))` for the rows. A matrix with a single row or
  a single column is still taken as a vector, and a data frame now stops with
  an informative message.

# WaveBased 2.6-3

## Bug fixes

* `wavedec()` and `waverec()` now stop when the length of `x` is not a power
  of two. Before, the check in the C code never fired, because `log2(n)` was
  truncated to an integer before it was tested, and `waverec()` had no check
  at all. A length such as 60, 100 or 1000 was accepted without an error: the
  periodic transform wrote past its work buffers, and R could abort at the
  next garbage collection. When `j0` was equal to the truncated level (for
  instance `n = 100` and `j0 = 6`), `x` was returned unchanged, with a warning
  that there was no decomposition. The C routines carry the same check, for
  direct callers of `.Call`. A vector of length one, returned unchanged with
  a warning before, is now rejected too.

# WaveBased 2.6-2

## Documentation

* The help of `wall()` now warns that, when the scaling functions are kept in
  the design (`boundary = "interval"`, or `j0 > 0`), the constant enters once
  per covariate and once more through the intercept: the design has `d`
  linearly dependent directions, and the intercept, the scaling coefficients
  returned by `coef()` and the level of each fitted component are not
  individually identified. The fitted log-odds, the probabilities and the
  classifier are not affected, and nothing changes in the fit. The help of
  `coef()` and of `plot(type = "components")` points to it, and a test checks
  the rank deficiency. With the defaults (periodized basis, `j0 = 0`) the
  scaling columns are dropped and the issue does not arise.
* The help of `plot.wall()` says which of its displays depend on the scaling
  coefficients in that case: the ranking of the covariates, the curves of
  `type = "path"`, the widths of the edges of `type = "network"` and, with
  `newx`, the values of the scaling nodes, of the components and of the
  intercept. It no longer says that only the sum of the components is
  identifiable, which does not hold for the defaults.

# WaveBased 2.6-1

## Changes in defaults

* The default grid of `J` in `cv.wall()` starts at `J = 2` (or at `j0 + 1`,
  if larger), and is now `J1:max(J1 + 1, ceiling(log2(n)/2))` with
  `J1 = max(j0 + 1, 2)`. At `J = 1`, the default boundary constant of the
  periodized basis is `eps = 1.9^(-1) > 0.5`, outside the range of the
  rescaling: the covariates were mapped, reversed, into a narrow band around
  1/2, and the candidate was degenerate, never selected. A grid supplied
  through `J` is used as given.

## Bug fixes

* `wall()` and `cv.wall()` warn when the default `eps` is 0.5 or more, which
  happens at `J = 1` with the periodized basis. The fit is the same as before;
  the rescaling it rests on is degenerate, and used to be applied silently. A
  value of `eps` supplied by the user is still validated, as before.

# WaveBased 2.6-0

This file starts here. Earlier versions are described by the commit history.

## Breaking changes

The arguments that choose and that hold the priors were renamed in
`bwregime()` and in `bwmixreg()`. They came in pairs, one naming a prior and
one holding its hyperparameters, and the old names did not show which pair a
name belonged to. Each pair is now named after its subject:

| subject               | which prior                     | its hyperparameters      |
|-----------------------|---------------------------------|--------------------------|
| wavelet coefficients  | `wavelet.prior` (new)           | `shrinkage`              |
| component scales      | `scale.prior` (was `cprior`)    | `components` (was `prior`) |

A call written against an earlier version now stops with an error that names
the replacement. It used to be swallowed by `...`, silently, and to return the
fit of the defaults; see the fix below.

## New features

* `bwregime()` gained the argument `wavelet.prior`, which chooses the prior of
  the wavelet coefficients of the mixture weights: `"spikeslab"` (the default,
  the prior of Motta and Montoril, 2026b) or `"horseshoe"` (Carvalho, Polson
  and Scott, 2010), written level by level, with a global scale of its own at
  every resolution. Its hyperparameter is the new entry `scale` of
  `shrinkage`, and its full conditionals are closed form under the scale
  mixture of Makalic and Schmidt (2016).

* `bwregime()` gained the argument `link`, which is now `"probit"` (the
  default, the link of the paper) or `"logit"`. The logit link is estimated by
  the Polya-Gamma augmentation of Polson, Scott and Windle (2013), whose latent
  variables are drawn exactly by Devroye's method, and completed to a common
  scale by the orthogonal data augmentation of Ghosh and Clyde (2011), which is
  what keeps the wavelet coefficients conditionally independent. Both links are
  exact and they identify the same regimes; the logit one costs between ten and
  twenty percent more per effective draw on smooth or moderately varying
  weights, and about twice as much on sharply switching ones.

* The level draws that `bwregime()` returns are named after what they hold:
  `pi` and `slab` under the spike and slab prior, `kappa` and `scale` under the
  horseshoe one. The component `inclusion` holds the posterior probability that
  a coefficient is non-null, or the posterior mean of its shrinkage weight,
  which are read on the same scale.

## Bug fixes

* `bwregime()` and `bwmixreg()` stop on an argument they do not have. Such an
  argument used to reach the plot method through `...`, where the graphics
  engine dropped it without a word, so that a misspelled or an outdated
  argument returned the fit of the defaults instead of the one that was asked
  for. Graphical parameters still pass through `...` as before.

## Internals

* The prior of the wavelet coefficients, its slab, the prior of the component
  scales and the link are each an interface with a registry, and each registry
  is read once, before the first sweep. A sweep never compares a model code nor
  searches a table.

* The chains of the spike and slab prior under the probit link are the ones of
  the previous version, bit for bit.
