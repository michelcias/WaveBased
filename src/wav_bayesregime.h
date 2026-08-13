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
 * - **The slab.** WBSlab collects the three operations that depend on the prior
 *   of the wavelet coefficients: the draw of the hyperparameters of a
 *   resolution level, the Bayes factor of the slab against the spike, and the
 *   draw of a coefficient. The Gaussian and the Laplace slabs of the article
 *   are two instances of it, registered in WBSlabTable(). A further slab is a
 *   further instance, and nothing else.
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
 * - **The component parameters.** They are drawn by WBDrawComponent(), selected
 *   by the @c cprior code of WBRegimeSpec. A heavier tailed prior for the
 *   component scales, such as a half-t through its scale mixture
 *   representation, is a sibling of that routine plus one auxiliary variable
 *   per component, held in the spec.
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

/** @name Model codes, shared with the R interface. */
/**@{*/
#define WB_LINK_PROBIT   1  /**< Probit link, by data augmentation.        */
#define WB_SLAB_GAUSSIAN 1  /**< Spike and slab with Gaussian slab (SSG).  */
#define WB_SLAB_LAPLACE  2  /**< Spike and slab with Laplace slab (SSL).   */
#define WB_CPRIOR_GAMMA  1  /**< Normal-gamma prior of the components.     */
/**@}*/

/**
 * @brief Everything the sweep needs to know about the model being fitted.
 *
 * @details Filled once from the arguments of C_BayesRegime() and passed around
 *          by pointer, so that a new hyperparameter is a new field and not a
 *          new argument of every routine.
 */
typedef struct {
  int n;              /**< Sample size.                                     */
  int npad;           /**< Padded sample size, a power of two.              */
  int J;              /**< Number of resolution levels, log2(npad).         */
  int link;           /**< Link code, one of the WB_LINK_* constants.       */
  int slab;           /**< Slab code, one of the WB_SLAB_* constants.       */
  int cprior;         /**< Component prior code, a WB_CPRIOR_* constant.    */
  double b0[2];       /**< Prior means of mu_1 and mu_2.                    */
  double B0[2];       /**< Prior variances of mu_1 and mu_2.                */
  double v0[2];       /**< Prior shapes of tau_1^2 and tau_2^2.             */
  double V0[2];       /**< Prior rates of tau_1^2 and tau_2^2.              */
  double zeta;        /**< First shape of the beta prior of pi_j.           */
  double rho;         /**< Second shape of the beta prior of pi_j.          */
  double kappa;       /**< Shape of the gamma prior of v_j^-2 or of a_j.    */
  double xi;          /**< Rate of the gamma prior of v_j^-2 or of a_j.     */
  double cut;         /**< Inclusion probabilities below it are set to zero.*/
} WBRegimeSpec;

/**
 * @brief The prior of the wavelet coefficients, seen by the sweep.
 *
 * @details The three members are the only places where the sampler depends on
 *          the slab. Both draws consume the same number of variates whatever
 *          the value of the inclusion indicator, which keeps the stream of the
 *          sampler independent of the data.
 */
typedef struct {
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
} WBSlab;

/**
 * @brief Runs the Gibbs sampler for the dynamic mixture weights.
 *
 * @details Each sweep draws the component parameters, permutes the labels
 *          whenever the identifying restriction mu_1 < mu_2 is violated, draws
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
 * @param[in] model   Integer SEXP with the components 'link', 'slab' and
 *                    'cprior', holding the codes of the model.
 * @param[in] hyper   List SEXP with the real components 'mean', 'var', 'shape'
 *                    and 'rate', of length two, holding the priors of the
 *                    component parameters, and the scalars 'zeta', 'rho',
 *                    'kappa' and 'xi', holding the priors of the sparsity and
 *                    of the slab parameters, and 'cut', the smallest posterior
 *                    inclusion probability that is not rounded down to zero.
 * @param[in] init    List SEXP with the real components 'mu' and 'tau2', of
 *                    length two, holding the initial values of the chains.
 * @param[in] mcmc    Integer SEXP of length 4: number of draws to keep,
 *                    burn-in, lag, and the reporting interval (zero for a
 *                    silent run).
 *
 * @return SEXP list with the retained draws of mu (nchain x 2), tau2
 *         (nchain x 2), alpha (nchain x n), pi (nchain x J) and the slab
 *         parameter (nchain x J), together with the posterior means of the
 *         allocation variables and of the inclusion indicators.
 */
SEXP C_BayesRegime(SEXP y, SEXP padidx, SEXP wavelet, SEXP model, SEXP hyper,
                   SEXP init, SEXP mcmc);

#endif /* WAVEBASED_WAV_BAYESREGIME_H */
