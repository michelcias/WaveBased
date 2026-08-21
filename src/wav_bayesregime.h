#ifndef WAVEBASED_WAV_BAYESREGIME_H
#define WAVEBASED_WAV_BAYESREGIME_H

/**
 * @file wav_bayesregime.h
 * @brief Gibbs sampler for the identification of regime switches.
 * @author Michel H. Montoril
 *
 * @details Implements Algorithm 1 of Motta and Montoril (2026b) for the dynamic
 *          two-component Gaussian mixture
 *          \f[ y_t = (1 - z_t)x_{1t} + z_t x_{2t}, \qquad
 *              x_{kt} \sim N(\mu_k, \tau_k^{-2}), \qquad
 *              z_t \sim \mathrm{Bern}(\alpha_t), \f]
 *          whose mixture weights are written as \f$\alpha = \Phi(W^T\theta)\f$,
 *          with W the discrete wavelet transform matrix. The vector of wavelet
 *          coefficients is estimated by the data augmentation of Albert and
 *          Chib (1993) under a spike and slab prior, so that the weights are
 *          recovered without imposing a functional form on them.
 *
 *          The component parameters and the allocation variables are drawn by
 *          WBDrawComponent() and WBDrawAlloc(), shared with the sampler of
 *          C_BayesMixReg(); what is proper to this file is the latent variable
 *          step, the prior of the wavelet coefficients and the link.
 *
 * @section extending Extending the sampler
 *
 * The sampler is written against three interfaces, so that the model can grow
 * without the driver being rewritten.
 *
 * - **The prior of the wavelet coefficients.** WBWavPrior holds the single
 *   operation the sweep asks of it: draw one whole resolution level, that is,
 *   its hyperparameters and its coefficients. The spike and slab prior of the
 *   article and the horseshoe prior are two instances of it, registered in
 *   WBWavPriorGet(). A level is the unit of the interface because a prior is
 *   free to share information among the coefficients of a level, which the
 *   sparsity parameter of the one and the global scale of the other both do.
 *
 * - **The slab.** Within the spike and slab prior, WBSlab collects the three
 *   operations that depend on the slab itself: the draw of the hyperparameters
 *   of a resolution level, the Bayes factor of the slab against the spike, and
 *   the draw of a coefficient. The Gaussian and the Laplace slabs of the
 *   article are two instances of it, registered in WBSlabGet(). It is the inner
 *   interface of one instance of WBWavPrior, and the horseshoe never consults
 *   it.
 *
 * - **The link.** Only the probit link of Albert and Chib is implemented, and
 *   it enters through WBDrawLatentProbit(), which returns the latent variable,
 *   and WBLinkMean(), which maps the linear predictor to a weight. A link whose
 *   augmentation is not homoscedastic (the logistic one under the Polya-Gamma
 *   scheme, say) delivers a working response with an observation-specific
 *   scale; every formula of the slabs is therefore written in terms of a
 *   general scale @c sigma, and never assumes the unit scale that the probit
 *   link happens to impose.
 *
 * - **The component parameters.** WBComponent collects the two operations that
 *   depend on their prior: the draw of the mean and of the precision of one
 *   component, and the refresh of the auxiliary variable that a scale mixture
 *   carries. The normal-gamma prior of the article and the half-t prior of the
 *   component standard deviations are two instances of it, registered in
 *   WBComponentGet(). The two operations are separate because the sweep permutes
 *   the labels between them, and each auxiliary must end the sweep paired with
 *   the precision that occupies its slot, under the prior of that slot. A prior
 *   without an auxiliary variable leaves the refresh returning what it is
 *   given, which is what the normal-gamma instance does. The interface does not
 *   depend on how the weights are built, so it lives in wav_bayesmix.h and is
 *   shared with C_BayesMixReg().
 *
 * The R entry point takes the model codes and the hyperparameters as coded
 * vectors and named lists, and not as positional arguments, so that new options
 * and new hyperparameters do not change its arity nor the registration table.
 *
 * @see Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
 *      polychotomous response data. Journal of the American Statistical
 *      Association, 88(422), 669-679.
 * @see Motta, F. C. and Montoril, M. H. (2026b). Identifying regime switches
 *      through Bayesian wavelet estimation: application to environmental and
 *      genetic data. Journal of Applied Statistics.
 */

#include <R.h>
#include <Rinternals.h>

#include "wav_bayesmix.h"  /* WBCompPrior, WBComponent, the WB_CPRIOR_* codes */

/** @name Model codes, shared with the R interface. */
/**@{*/
#define WB_LINK_PROBIT   1  /**< Probit link, by data augmentation.        */
#define WB_LINK_LOGIT    2  /**< Logit link, by Polya-Gamma augmentation.  */
#define WB_WPRIOR_SS     1  /**< Spike and slab prior of the coefficients.  */
#define WB_WPRIOR_HS     2  /**< Horseshoe prior of the coefficients.       */
#define WB_SLAB_GAUSSIAN 1  /**< Spike and slab with Gaussian slab (SSG).  */
#define WB_SLAB_LAPLACE  2  /**< Spike and slab with Laplace slab (SSL).   */
/**@}*/
/* The codes of the prior of the component parameters, WB_CPRIOR_*, are shared
 * with C_BayesMixReg() and live in wav_bayesmix.h. */

/* The specification holds the interfaces that the model codes select, and the
 * interfaces are written against the specification, so both are declared before
 * either is defined. Every registry is read once, before the first sweep, and
 * what the sweep sees from then on is a pointer: no code is compared and no
 * table is searched while the chain runs. */
typedef struct WBRegimeSpec WBRegimeSpec;
typedef struct WBSlab WBSlab;
typedef struct WBWavPrior WBWavPrior;
typedef struct WBLink WBLink;

/**
 * @brief The slab of the spike and slab prior.
 *
 * @details The three members are the only places where the sampler depends on
 *          the slab. Both draws consume the same number of variates whatever
 *          the value of the inclusion indicator, which keeps the stream of the
 *          sampler independent of the data.
 */
struct WBSlab {
  const char *name;
  /**
   * @brief Draws the hyperparameters of one resolution level.
   * @param[in]  spec  Model specification.
   * @param[in]  theta Current coefficients of the level, of length m.
   * @param[in]  gam   Current inclusion indicators of the level, of length m.
   * @param[in]  m     Number of coefficients of the level.
   * @param[out] pi    Drawn sparsity parameter.
   * @param[out] par   Drawn slab parameter: the variance v_j^2 of the Gaussian
   *                   slab, or the scale a_j of the Laplace slab.
   */
  void (*level_pars)(const WBRegimeSpec *spec, const double *theta,
                     const int *gam, int m, double *pi, double *par);
  /**
   * @brief Log of the ratio between the slab and the spike marginals.
   * @param[in] d     Empirical wavelet coefficient.
   * @param[in] sigma Scale of the latent regression at that coefficient.
   * @param[in] par   Slab parameter of the level.
   * @return log{g(d)/phi_sigma(d)}, the log Bayes factor of expressions (3)
   *         and (4) of the article.
   */
  double (*log_bf)(double d, double sigma, double par);
  /**
   * @brief Draws a wavelet coefficient from its full conditional.
   * @param[in] d     Empirical wavelet coefficient.
   * @param[in] sigma Scale of the latent regression at that coefficient.
   * @param[in] par   Slab parameter of the level.
   * @param[in] gam   Inclusion indicator just drawn.
   * @return The draw, which is zero whenever @p gam is zero.
   */
  double (*draw_theta)(double d, double sigma, double par, int gam);
};

/**
 * @brief The prior of the wavelet coefficients, seen by the sweep.
 *
 * @details The single member is the only place where the sampler depends on how
 *          the coefficients are shrunk. It is handed one resolution level and
 *          leaves it complete: the hyperparameters of the level and every
 *          coefficient of it are drawn from their full conditionals, given the
 *          state the level is handed in.
 *
 *          The two priors keep different states between sweeps, and the
 *          interface carries the union of them. A prior writes into what it
 *          keeps and leaves the rest untouched: the spike and slab prior uses
 *          @p gam and ignores @p loc, and the horseshoe prior does the reverse.
 *          What both must fill is @p w, the per-coefficient weight that the
 *          sweep averages into the 'inclusion' output, so that the summary has
 *          the same reading under either prior: the posterior probability that
 *          a coefficient is in the model, and the posterior mean of the
 *          shrinkage weight it receives.
 */
struct WBWavPrior {
  const char *name;
  /**
   * @brief Draws the hyperparameters and the coefficients of one level.
   * @param[in]     spec  Model specification.
   * @param[in]     d     Empirical wavelet coefficients of the level, length m.
   * @param[in]     sigma Scale of the latent regression at the level.
   * @param[in]     m     Number of coefficients of the level.
   * @param[in,out] theta Coefficients, handed in at their previous draw.
   * @param[in,out] gam   Inclusion indicators of the level.
   * @param[in,out] loc   Local scales of the level, the squared lambda_k of the
   *                      horseshoe prior.
   * @param[out]    w     Weight of each coefficient, in the unit interval.
   * @param[in,out] pi    Sparsity parameter pi_j, or the mean shrinkage weight
   *                      of the level.
   * @param[in,out] par   Slab parameter (the variance v_j^2 of the Gaussian
   *                      slab or the scale a_j of the Laplace one), or the
   *                      squared global scale tau_j^2 of the horseshoe prior.
   *                      It is handed in at its previous draw, which the
   *                      horseshoe prior needs and the spike and slab one
   *                      overwrites.
   */
  void (*draw_level)(const WBRegimeSpec *spec, const double *d, double sigma,
                     int m, double *theta, int *gam, double *loc, double *w,
                     double *pi, double *par);
};

/**
 * @brief The link between the wavelet expansion and the mixture weights.
 *
 * @details The two members are the only places where the sampler depends on the
 *          link, and both are called once per sweep and not once per
 *          observation, so that the choice costs one indirect call and never a
 *          comparison inside a loop.
 *
 *          What separates the two links is homoscedasticity. The probit
 *          augmentation of Albert and Chib (1993) leaves the latent regression
 *          with a common scale, and the wavelet coefficients of the working
 *          response are therefore independent. The logit augmentation does not,
 *          under the Polya-Gamma scheme of Polson, Scott and Windle (2013) or
 *          under any other: the working response carries a precision of its
 *          own at every observation, the posterior precision of the
 *          coefficients stops being diagonal, and none of the coefficient-wise
 *          formulas of the priors apply. The diagonal is restored by the
 *          orthogonal data augmentation of Ghosh and Clyde (2011), described in
 *          WBDrawWorkingLogit(), and the common scale it completes the data to
 *          is what the first member returns.
 */
struct WBLink {
  const char *name;
  /**
   * @brief Draws the working response of one sweep.
   * @param[in]  spec  Model specification.
   * @param[in]  fit   Linear predictor at each position of the padded grid.
   * @param[in]  z     Allocation variables, on the grid of the sample.
   * @param[in]  idx   Padding index, gathering z onto the padded grid.
   * @param[in]  npad  Padded sample size.
   * @param[out] v     Working response.
   * @param[out] omega Scratch of length npad, which a heteroscedastic
   *                   augmentation uses for its weights and the probit link
   *                   leaves untouched.
   * @return The common scale of the working response.
   */
  double (*draw_working)(const WBRegimeSpec *spec, const double *fit,
                         const int *z, const int *idx, int npad, double *v,
                         double *omega);
  /**
   * @brief Maps the linear predictor to the mixture weights of the sample.
   * @param[in]  fit   Linear predictor, of which the first n entries are read.
   * @param[in]  n     Sample size.
   * @param[out] alpha Mixture weights.
   */
  void (*weights)(const double *fit, int n, double *alpha);
};

/**
 * @brief Everything the sweep needs to know about the model being fitted.
 *
 * @details Filled once from the arguments of C_BayesRegime() and passed around
 *          by pointer, so that a new hyperparameter is a new field and not a
 *          new argument of every routine. The interfaces the model codes select
 *          are resolved into it before the first sweep.
 */
struct WBRegimeSpec {
  int n;              /**< Sample size.                                     */
  int npad;           /**< Padded sample size, a power of two.              */
  int J;              /**< Number of resolution levels, log2(npad).         */
  const WBLink *link; /**< The link, resolved.                              */
  const WBSlab *slab; /**< The slab of the spike and slab prior, resolved.  */
  WBCompPrior comp[2];/**< Priors of the parameters of the two components.  */
  double zeta;        /**< First shape of the beta prior of pi_j.           */
  double rho;         /**< Second shape of the beta prior of pi_j.          */
  double kappa;       /**< Shape of the gamma prior of v_j^-2 or of a_j.    */
  double xi;          /**< Rate of the gamma prior of v_j^-2 or of a_j.     */
  double cut;         /**< Inclusion probabilities below it are set to zero.*/
  double hscale;      /**< Scale of the half-Cauchy prior of tau_j.         */
};

/**
 * @brief Runs the Gibbs sampler for the dynamic mixture weights.
 *
 * @details Each sweep draws the component parameters, permutes the labels
 *          whenever the identifying restriction mu_1 < mu_2 is violated,
 *          refreshes the auxiliary variables of their prior, draws
 *          the allocation variables, the latent variables of the augmentation,
 *          the hyperparameters of each resolution level, the inclusion
 *          indicators and the wavelet coefficients, and rebuilds the mixture
 *          weights through the link. Only the sweeps after the burn-in period,
 *          taken every 'lag' sweeps, are stored. No memory is allocated inside
 *          the loop, and exactly two wavelet transforms are performed by sweep.
 *
 * @param[in] y       Real SEXP: observed mixture of length n.
 * @param[in] padidx  Integer SEXP of length npad (a power of two) with
 *                    zero-based positions of y, describing how the sample is
 *                    extended to a dyadic length. The first n entries must be
 *                    0, 1, ..., n - 1.
 * @param[in] wavelet List SEXP with the components 'family' (integer:
 *                    1 = Daublets, 2 = Symmlets, 3 = Coiflets, 4 = custom),
 *                    'filter.size' (integer) and 'filter' (real, used only when
 *                    the family is custom).
 * @param[in] model   Integer SEXP with the components 'link', 'wprior', 'slab'
 *                    and 'cprior', holding the codes of the model. The 'slab'
 *                    code is read only when 'wprior' selects the spike and slab
 *                    prior.
 * @param[in] hyper   List SEXP with the real components 'mean' and 'var', of
 *                    length two, holding the normal priors of the component
 *                    means, 'shape' and 'rate', of length two, holding the gamma
 *                    priors of the component precisions, 'df' and 'scale', of
 *                    length two, holding the half-t priors of the component
 *                    standard deviations, the scalars 'zeta', 'rho', 'kappa'
 *                    and 'xi', holding the priors of the sparsity and of the
 *                    slab parameters, 'cut', the smallest posterior inclusion
 *                    probability that is not rounded down to zero, and
 *                    'hscale', the scale of the half-Cauchy prior of the global
 *                    scale of a level. Only the entries the 'cprior' and
 *                    'wprior' codes select are read.
 * @param[in] init    List SEXP with the real components 'mu' and 'tau2', of
 *                    length two, holding the initial values of the chains.
 * @param[in] mcmc    Integer SEXP of length 4: number of draws to keep,
 *                    burn-in, lag, and the reporting interval (zero for a
 *                    silent run).
 *
 * @return SEXP list with the retained draws of mu (nchain x 2), tau2
 *         (nchain x 2), alpha (nchain x n), pi (nchain x J) and the slab
 *         parameter (nchain x J), together with the posterior means of the
 *         allocation variables and of the weight of each coefficient. The last
 *         three carry the reading that the 'wprior' code gives them: the
 *         sparsity parameter, the slab parameter and the inclusion indicator
 *         under the spike and slab prior, and the mean shrinkage weight of a
 *         level, its squared global scale and the shrinkage weight of a
 *         coefficient under the horseshoe prior.
 */
SEXP C_BayesRegime(SEXP y, SEXP padidx, SEXP wavelet, SEXP model, SEXP hyper,
                   SEXP init, SEXP mcmc);

#endif /* WAVEBASED_WAV_BAYESREGIME_H */
