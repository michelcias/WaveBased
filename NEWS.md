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
