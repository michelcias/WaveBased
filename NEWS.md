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
  The message now names the family. The paragraph of `wbasis()` that described the run-time
  construction for `family = "Own"` was corrected.

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
