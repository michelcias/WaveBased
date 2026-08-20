/**
 * @file wav_bayesregime.c
 * @brief Gibbs sampler for the identification of regime switches.
 * @author Michel H. Montoril
 *
 * @details Implements Algorithm 1 of Motta and Montoril (2026b). The sampler
 *          was implemented based on the reference R code of the article,
 *          available at
 *          https://github.com/flaviamotta/Regime-switches-and-Bayesian-wavelet-estimation,
 *          with the whole sweep moved into C. See wav_bayesregime.h for the
 *          interfaces through which the model can be extended.
 */

#include <math.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "wav_bayesregime.h"
#include "wav_bayesmix.h"
#include "wav_transform.h"
#include "wav_utilities.h"

//==============================================================================
// NUMERICAL HELPERS
//==============================================================================

/**
 * @brief Logarithm of a sum of two exponentials, without overflowing.
 *
 * @param[in] a First logarithm.
 * @param[in] b Second logarithm.
 * @return log(exp(a) + exp(b)).
 */
static double WBLogSumExp(double a, double b){

  return (a > b) ? a + log1p(exp(b - a)) : b + log1p(exp(a - b));
}

/**
 * @brief Draws a normal variate truncated at zero.
 *
 * @details Inverts the distribution function on the logarithmic scale, so that
 *          the draw stays accurate however far the mean lies from the truncation
 *          point. The plain inversion loses every digit once the truncated tail
 *          falls below the resolution of the double precision, which is the
 *          rule and not the exception here: the mean is the linear predictor of
 *          the probit link, and the regimes it describes are precisely the
 *          regions where that predictor is far from zero.
 *
 *          Exactly one uniform variate is consumed by call.
 *
 * @param[in] m        Mean of the underlying normal distribution.
 * @param[in] sigma    Standard deviation of the underlying normal distribution.
 * @param[in] positive Nonzero to truncate to the positive half line, zero to
 *                     truncate to the negative one.
 * @return The truncated draw.
 */
static double WBDrawTruncNorm(double m, double sigma, int positive){

  /* unif_rand() returns a value in the open unit interval, so the logarithm
   * below is always finite. */
  double lp = log1p(-unif_rand());

  if(positive){
    lp += pnorm(m/sigma, 0.0, 1.0, 1, 1);
    return m - sigma*qnorm(lp, 0.0, 1.0, 1, 1);
  }

  lp += pnorm(-m/sigma, 0.0, 1.0, 1, 1);
  return m + sigma*qnorm(lp, 0.0, 1.0, 1, 1);
}

//==============================================================================
// THE LINK
//==============================================================================

/**
 * @brief Draws the latent variable of the probit augmentation.
 *
 * @details This is expression (12) of the article: the latent variable is
 *          normal with mean given by the linear predictor and unit variance,
 *          truncated to the half line indicated by the allocation variable.
 *
 * @param[in] eta Linear predictor at that index.
 * @param[in] z   Allocation variable at that index, 0 or 1.
 * @return The latent variable.
 */
static double WBDrawLatentProbit(double eta, int z){

  return WBDrawTruncNorm(eta, 1.0, z);
}

/**
 * @brief Maps the linear predictor to a mixture weight, as in expression (7).
 *
 * @param[in] eta Linear predictor.
 * @return The mixture weight.
 */
static double WBProbitMean(double eta){

  return pnorm(eta, 0.0, 1.0, 1, 0);
}

//==============================================================================
// THE SPIKE AND SLAB PRIORS
//==============================================================================

/**
 * @brief Rounds a posterior inclusion probability down to zero below a cut.
 *
 * @param[in] p   Posterior inclusion probability.
 * @param[in] cut Smallest probability that is left untouched.
 * @return Zero when the probability lies below the cut, and the probability
 *         itself otherwise. A cut of zero leaves everything untouched.
 */
static double WBCut(double p, double cut){

  return (p < cut) ? 0.0 : p;
}

/**
 * @brief Draws the sparsity parameter of one resolution level.
 *
 * @details Expression (14) of the article, common to every slab: the beta prior
 *          is conjugate to the inclusion indicators of the level.
 *
 * @param[in] spec Model specification.
 * @param[in] n1   Number of non-null coefficients of the level.
 * @param[in] m    Number of coefficients of the level.
 * @return The drawn sparsity parameter.
 */
static double WBDrawSparsity(const WBRegimeSpec *spec, int n1, int m){

  return rbeta(spec->zeta + n1, spec->rho + m - n1);
}

/* --- Gaussian slab (SSG) --------------------------------------------------- */

static void WBLevelParsNormal(const WBRegimeSpec *spec, const double *theta,
                              const int *gam, int m, double *pi, double *par){

  int k, n1 = 0;
  double ss = 0.0;

  for(k = 0; k < m; k++){
    n1 += gam[k];
    ss += theta[k]*theta[k];
  }

  *pi = WBDrawSparsity(spec, n1, m);

  /* Expression (15) draws the precision of the slab, whose variance is the
   * parameter the remaining routines work with. */
  *par = 1.0/rgamma(spec->kappa + 0.5*m, 1.0/(spec->xi + 0.5*ss));
}

static double WBLogBFNormal(double d, double sigma, double par){

  double s2 = sigma*sigma, r = par/s2;

  /* The marginal of the slab is N(0, v^2 + sigma^2), that of the spike is
   * N(0, sigma^2), and this is the logarithm of the ratio of their densities. */
  return 0.5*(d*d*r/(s2*(1.0 + r)) - log1p(r));
}

static double WBDrawThetaNormal(double d, double sigma, double par, int gam){

  double lambda = par/(par + sigma*sigma);
  double x = lambda*d + sigma*sqrt(lambda)*norm_rand();

  return gam ? x : 0.0;
}

/* --- Laplace slab (SSL) ---------------------------------------------------- */

static void WBLevelParsLaplace(const WBRegimeSpec *spec, const double *theta,
                               const int *gam, int m, double *pi, double *par){

  int k, n1 = 0;
  double sa = 0.0;

  for(k = 0; k < m; k++){
    n1 += gam[k];
    sa += fabs(theta[k]);
  }

  *pi = WBDrawSparsity(spec, n1, m);
  *par = rgamma(spec->kappa + m, 1.0/(spec->xi + sa));
}

static double WBLogBFLaplace(double d, double sigma, double par){

  double zm = d/sigma - par*sigma, zp = d/sigma + par*sigma;

  /* The convolution between the Laplace slab and the normal noise, divided by
   * the density of the noise, written with the two Mills ratios on the
   * logarithmic scale. Their asymptote is taken care of by pnorm(), so no tail
   * approximation is needed. */
  return log(0.5*par*sigma)
         + WBLogSumExp(pnorm(zm, 0.0, 1.0, 1, 1) - dnorm(zm, 0.0, 1.0, 1),
                       pnorm(zp, 0.0, 1.0, 0, 1) - dnorm(zp, 0.0, 1.0, 1));
}

static double WBDrawThetaLaplace(double d, double sigma, double par, int gam){

  double shift = par*sigma*sigma, m, x;
  double lup = -par*d + pnorm((d - shift)/sigma, 0.0, 1.0, 1, 1);
  double ldn =  par*d + pnorm((d + shift)/sigma, 0.0, 1.0, 0, 1);
  int positive;

  /* Expression (5): the non-null component of the posterior is a mixture of two
   * truncated normal distributions, whose weight is evaluated on the
   * logarithmic scale. */
  positive = (unif_rand() < 1.0/(1.0 + exp(ldn - lup)));
  m = positive ? d - shift : d + shift;
  x = WBDrawTruncNorm(m, sigma, positive);

  return gam ? x : 0.0;
}

/* --- The registry of slabs ------------------------------------------------- */

/**
 * @brief Returns the slab requested by a model code.
 *
 * @param[in] code One of the WB_SLAB_* constants.
 * @return Pointer to the corresponding interface. The function does not return
 *         when the code is unknown.
 */
static const WBSlab *WBSlabGet(int code){

  static const WBSlab table[] = {
    {"gaussian", WBLevelParsNormal,  WBLogBFNormal,  WBDrawThetaNormal},
    {"laplace",  WBLevelParsLaplace, WBLogBFLaplace, WBDrawThetaLaplace}
  };

  if(code == WB_SLAB_GAUSSIAN) return &table[0];
  if(code == WB_SLAB_LAPLACE)  return &table[1];

  error("Unknown slab code.");
  return NULL;
}

/**
 * @brief Draws one resolution level under the spike and slab prior.
 *
 * @details Expressions (14) to (16) of the article. The hyperparameters of the
 *          level are drawn from its previous state, and each coefficient then
 *          receives its inclusion indicator and its value. The local scales
 *          @p loc are not part of this prior and are left untouched.
 *
 * @copydetails WBWavPrior::draw_level
 */
static void WBDrawLevelSS(const WBRegimeSpec *spec, const double *d,
                          double sigma, int m, double *theta, int *gam,
                          double *loc, double *w, double *pi, double *par){

  const WBSlab *slab = WBSlabGet(spec->slab);
  double lodds, prob;
  int k;

  (void) loc;

  slab->level_pars(spec, theta, gam, m, pi, par);

  lodds = log(*pi) - log1p(-(*pi));

  for(k = 0; k < m; k++){

    /* Expression (16), on the logarithmic scale. Coefficients whose inclusion
     * probability does not reach the cut are excluded outright, which keeps the
     * sparsity parameters from drifting upwards: a coefficient of pure noise
     * carries a Bayes factor of about one, so that the beta full conditional of
     * pi_j would raise it by about one coefficient at every sweep, level after
     * level, until nothing is shrunk any more. The uniform variate is consumed
     * either way, so that the stream of the sampler does not depend on the
     * data. */
    prob = WBCut(1.0/(1.0 + exp(-lodds - slab->log_bf(d[k], sigma, *par))),
                 spec->cut);
    gam[k] = (unif_rand() < prob);
    theta[k] = slab->draw_theta(d[k], sigma, *par, gam[k]);

    /* The indicator itself, and not the probability it was drawn from, so that
     * the summary is the posterior mean of what the sweep visits. */
    w[k] = (double) gam[k];
  }
}

/* --- The horseshoe prior --------------------------------------------------- */

/**
 * @brief Draws one resolution level under the horseshoe prior.
 *
 * @details The prior of Carvalho, Polson and Scott (2010) written level by
 *          level,
 *          \f[ \theta_k \mid \lambda_k, \tau_j \sim
 *              N(0, \sigma^2\lambda_k^2\tau_j^2), \qquad
 *              \lambda_k \sim C^+(0, 1), \qquad
 *              \tau_j \sim C^+(0, A), \f]
 *          so that each level keeps a global scale of its own, as the sparsity
 *          parameter pi_j of the spike and slab prior does, and the amount of
 *          shrinkage is learnt separately at every resolution.
 *
 *          Every full conditional is closed form under the inverse gamma scale
 *          mixture of Makalic and Schmidt (2016), which is the representation
 *          the sampler of the component scales already uses: a half-Cauchy
 *          variable x with scale A satisfies
 *          \f[ x^2 \mid a \sim \mathrm{IG}(1/2, 1/a), \qquad
 *              a \sim \mathrm{IG}(1/2, 1/A^2), \f]
 *          the two shapes being half the single degree of freedom of the
 *          half-Cauchy distribution. The full conditionals of the auxiliary
 *          variables therefore have shape one, that half joined by the half a
 *          single Gaussian observation contributes. The auxiliary variables are
 *          drawn afresh from the scales they belong to, so nothing but the
 *          scales themselves has to be carried between sweeps.
 *
 *          Nothing is set to zero: the shrinkage is continuous, and it is the
 *          weight
 *          \f$\kappa_k = \lambda_k^2\tau_j^2/(1 + \lambda_k^2\tau_j^2)\f$
 *          given to the empirical coefficient that plays the part of the
 *          inclusion indicator. The indicators are held at one for that reason,
 *          and the cut has nothing to act on.
 *
 * @copydetails WBWavPrior::draw_level
 *
 * @see Carvalho, C. M., Polson, N. G. and Scott, J. G. (2010). The horseshoe
 *      estimator for sparse signals. Biometrika, 97(2), 465-480.
 * @see Makalic, E. and Schmidt, D. F. (2016). A simple sampler for the
 *      horseshoe estimator. IEEE Signal Processing Letters, 23(1), 179-182.
 */
static void WBDrawLevelHS(const WBRegimeSpec *spec, const double *d,
                          double sigma, int m, double *theta, int *gam,
                          double *loc, double *w, double *pi, double *par){

  double s2 = sigma*sigma, tau2 = *par, aux, ss = 0.0, kap, sw = 0.0;
  int k;

  /* The coefficients of the level, standardized by their local scales, which
   * is what the rate of the full conditional of the global scale is built
   * from. Its shape is (m + 1)/2: one half from the prior, and one half for
   * each coefficient the level holds. */
  for(k = 0; k < m; k++)
    ss += theta[k]*theta[k]/loc[k];

  aux = 1.0/rgamma(1.0, 1.0/(1.0/(spec->hscale*spec->hscale) + 1.0/tau2));
  tau2 = 1.0/rgamma(0.5*(m + 1.0), 1.0/(1.0/aux + 0.5*ss/s2));

  for(k = 0; k < m; k++){

    /* The local scale of the coefficient, and then the coefficient itself. The
     * full conditional is the one of the Gaussian slab of variance
     * sigma^2*lambda_k^2*tau_j^2, which is where the two priors meet. */
    aux = 1.0/rgamma(1.0, 1.0/(1.0 + 1.0/loc[k]));
    loc[k] = 1.0/rgamma(1.0, 1.0/(1.0/aux + 0.5*theta[k]*theta[k]/(tau2*s2)));

    kap = loc[k]*tau2/(1.0 + loc[k]*tau2);
    theta[k] = kap*d[k] + sigma*sqrt(kap)*norm_rand();

    gam[k] = 1;
    w[k] = kap;
    sw += kap;
  }

  *par = tau2;

  /* The reading that pi_j receives under this prior: the average weight the
   * coefficients of the level give the data, which is what the proportion of
   * non-null coefficients is under the spike and slab prior. */
  *pi = sw/m;
}

/* --- The registry of coefficient priors ------------------------------------ */

/**
 * @brief Returns the prior of the wavelet coefficients requested by a code.
 *
 * @param[in] code One of the WB_WPRIOR_* constants.
 * @return Pointer to the corresponding interface. The function does not return
 *         when the code is unknown.
 */
static const WBWavPrior *WBWavPriorGet(int code){

  static const WBWavPrior table[] = {
    {"spikeslab", WBDrawLevelSS},
    {"horseshoe", WBDrawLevelHS}
  };

  if(code == WB_WPRIOR_SS) return &table[0];
  if(code == WB_WPRIOR_HS) return &table[1];

  error("Unknown code for the prior of the wavelet coefficients.");
  return NULL;
}

//==============================================================================
// READING THE ARGUMENTS
//==============================================================================

/**
 * @brief Position of a named component of a vector or of a list.
 *
 * @param[in] x    Named vector or list.
 * @param[in] name Name to look for.
 * @return Its position. The function does not return when the name is absent.
 */
static int WBWhich(SEXP x, const char *name){

  SEXP nms = getAttrib(x, R_NamesSymbol);
  int i;

  if(nms == R_NilValue)
    error("A named vector or list was expected.");

  for(i = 0; i < length(nms); i++)
    if(strcmp(CHAR(STRING_ELT(nms, i)), name) == 0)
      return i;

  error("The component '%s' is missing.", name);
  return -1;
}

/** @brief The k-th value of a real component of a list. */
static double WBReal(SEXP list, const char *name, int k){

  SEXP el = VECTOR_ELT(list, WBWhich(list, name));

  if(TYPEOF(el) != REALSXP || length(el) <= k)
    error("The component '%s' must be a numeric vector of length %d or more.",
          name, k + 1);

  return REAL(el)[k];
}

/** @brief The value of a named integer component of a vector. */
static int WBInt(SEXP x, const char *name){

  if(TYPEOF(x) != INTSXP)
    error("An integer vector was expected.");

  return INTEGER(x)[WBWhich(x, name)];
}

//==============================================================================
// THE SAMPLER
//==============================================================================

SEXP C_BayesRegime(SEXP y, SEXP padidx, SEXP wavelet, SEXP model, SEXP hyper,
                   SEXP init, SEXP mcmc){

  int i, j, k, t, m, lo, keep, J = 0, N;
  int n = length(y), npad = length(padidx), total;
  int nchain = INTEGER(mcmc)[0], burn = INTEGER(mcmc)[1];
  int lag = INTEGER(mcmc)[2], nverb = INTEGER(mcmc)[3];
  int *idx = INTEGER(padidx), *z, *gam;
  double *ry = REAL(y), *rwfilter, *work, *lat, *dstar, *theta, *fit, *alpha;
  double *pilev, *parlev, *loc, *wgt;
  double *rmu, *rtau2, *ralpha, *rpi, *rpar, *rzbar, *rincl;
  double mu[2], tau2[2], sd[2], aux[2], tmp;
  /* The probit augmentation leaves the latent regression with a unit scale.
   * The priors of the coefficients are written for a general one, which is what
   * a link with a heteroscedastic augmentation would supply. */
  const double sigma = 1.0;
  WBRegimeSpec spec;
  const WBWavPrior *wprior;
  const WBComponent *comp;
  SEXP out, nms, wutils;

  while((1 << (J + 1)) <= npad) J++;

  if((1 << J) != npad)
    error("The padded sample size must be a power of two.");
  if(npad < n)
    error("The padded sample size must not be smaller than the sample size.");
  if(nchain < 1 || burn < 0 || lag < 1)
    error("Invalid MCMC settings.");

  total = burn + lag*nchain;

  /* --- The model specification ------------------------------------------- */
  spec.n = n;
  spec.npad = npad;
  spec.J = J;
  spec.link = WBInt(model, "link");
  spec.wprior = WBInt(model, "wprior");
  spec.slab = WBInt(model, "slab");
  spec.cprior = WBInt(model, "cprior");

  if(spec.link != WB_LINK_PROBIT)
    error("Only the probit link is available.");

  wprior = WBWavPriorGet(spec.wprior);
  comp = WBComponentGet(spec.cprior);

  for(k = 0; k < 2; k++){
    spec.comp[k].b0 = WBReal(hyper, "mean", k);
    spec.comp[k].B0 = WBReal(hyper, "var", k);
    spec.comp[k].v0 = WBReal(hyper, "shape", k);
    spec.comp[k].V0 = WBReal(hyper, "rate", k);
    spec.comp[k].nu0 = WBReal(hyper, "df", k);
    spec.comp[k].A0 = WBReal(hyper, "scale", k);
  }

  spec.zeta = WBReal(hyper, "zeta", 0);
  spec.rho = WBReal(hyper, "rho", 0);
  spec.kappa = WBReal(hyper, "kappa", 0);
  spec.xi = WBReal(hyper, "xi", 0);
  spec.cut = WBReal(hyper, "cut", 0);
  spec.hscale = WBReal(hyper, "hscale", 0);

  if(spec.zeta <= 0.0 || spec.rho <= 0.0 || spec.kappa <= 0.0 || spec.xi <= 0.0)
    error("The hyperparameters of the spike and slab prior must be positive.");
  if(spec.cut < 0.0 || spec.cut >= 1.0)
    error("The cut of the inclusion probabilities must lie in the unit interval.");
  if(spec.hscale <= 0.0)
    error("The scale of the horseshoe prior must be positive.");

  wutils = PROTECT(WavUtilities(VECTOR_ELT(wavelet, WBWhich(wavelet, "family")),
                                VECTOR_ELT(wavelet, WBWhich(wavelet, "filter.size")),
                                VECTOR_ELT(wavelet, WBWhich(wavelet, "filter"))));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  /* --- Working space, allocated once for the whole run ------------------- */
  work   = (double *) R_alloc(WB_WAVE_WORK(npad), sizeof(double));
  lat    = (double *) R_alloc(npad, sizeof(double));
  dstar  = (double *) R_alloc(npad, sizeof(double));
  theta  = (double *) R_alloc(npad, sizeof(double));
  fit    = (double *) R_alloc(npad, sizeof(double));
  alpha  = (double *) R_alloc(n, sizeof(double));
  pilev  = (double *) R_alloc(J, sizeof(double));
  parlev = (double *) R_alloc(J, sizeof(double));
  loc    = (double *) R_alloc(npad, sizeof(double));
  wgt    = (double *) R_alloc(npad, sizeof(double));
  z      = (int *)    R_alloc(n, sizeof(int));
  gam    = (int *)    R_alloc(npad, sizeof(int));

  /* --- The retained draws ------------------------------------------------ */
  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, allocMatrix(REALSXP, nchain, 2));     /* mu       */
  SET_VECTOR_ELT(out, 1, allocMatrix(REALSXP, nchain, 2));     /* tau2     */
  SET_VECTOR_ELT(out, 2, allocMatrix(REALSXP, nchain, n));     /* alpha    */
  SET_VECTOR_ELT(out, 3, allocMatrix(REALSXP, nchain, J));     /* pi       */
  SET_VECTOR_ELT(out, 4, allocMatrix(REALSXP, nchain, J));     /* slab     */
  SET_VECTOR_ELT(out, 5, allocVector(REALSXP, n));             /* z mean   */
  SET_VECTOR_ELT(out, 6, allocVector(REALSXP, npad));          /* gam mean */

  rmu    = REAL(VECTOR_ELT(out, 0));
  rtau2  = REAL(VECTOR_ELT(out, 1));
  ralpha = REAL(VECTOR_ELT(out, 2));
  rpi    = REAL(VECTOR_ELT(out, 3));
  rpar   = REAL(VECTOR_ELT(out, 4));
  rzbar  = REAL(VECTOR_ELT(out, 5));
  rincl  = REAL(VECTOR_ELT(out, 6));

  for(t = 0; t < n; t++)
    rzbar[t] = 0.0;
  for(t = 0; t < npad; t++)
    rincl[t] = 0.0;

  /* --- Initial values ---------------------------------------------------- */
  mu[0] = WBReal(init, "mu", 0);
  mu[1] = WBReal(init, "mu", 1);
  tau2[0] = WBReal(init, "tau2", 0);
  tau2[1] = WBReal(init, "tau2", 1);

  /* The chain of the coefficients starts from the null model, which places the
   * weights at one half. The scaling coefficient carries a diffuse prior, so it
   * is always part of the model and has no inclusion indicator of its own. */
  for(t = 0; t < npad; t++){
    theta[t] = 0.0;
    fit[t] = 0.0;
    gam[t] = 0;
    loc[t] = 1.0;
    wgt[t] = 0.0;
  }
  gam[0] = 1;
  wgt[0] = 1.0;

  /* The global scale of a level starts at one, the median of a standard
   * half-Cauchy variable, which is the state the horseshoe prior is handed at
   * the first sweep. The sparsity parameter beside it is a placeholder: the
   * spike and slab prior draws both before it reads them, and is left
   * untouched by this. */
  for(j = 0; j < J; j++){
    pilev[j] = 0.5;
    parlev[j] = 1.0;
  }

  for(t = 0; t < n; t++)
    alpha[t] = WBProbitMean(0.0);

  GetRNGstate();

  /* The auxiliary variables of the prior of the component parameters start from
   * their own full conditional given the initial precisions, which is the draw
   * the sweep performs anyway. Under the normal-gamma prior there is nothing to
   * draw, and the two calls leave the stream untouched. */
  aux[0] = comp->refresh(&spec.comp[0], tau2[0], 0.0);
  aux[1] = comp->refresh(&spec.comp[1], tau2[1], 0.0);

  sd[0] = 1.0/sqrt(tau2[0]);
  sd[1] = 1.0/sqrt(tau2[1]);
  WBDrawAlloc(ry, alpha, n, mu, sd, z);

  //============================================================================
  // GIBBS SAMPLING
  //============================================================================

  keep = 0;

  for(i = 1; i <= total; i++){

    R_CheckUserInterrupt();

    /* --- Component parameters, and the label switching constraint -------- */
    comp->draw(&spec.comp[0], ry, z, n, 0, aux[0], &mu[0], &tau2[0]);
    comp->draw(&spec.comp[1], ry, z, n, 1, aux[1], &mu[1], &tau2[1]);

    if(mu[1] < mu[0]){
      tmp = mu[0];    mu[0] = mu[1];       mu[1] = tmp;
      tmp = tau2[0];  tau2[0] = tau2[1];   tau2[1] = tmp;
    }

    /* The auxiliary variables are refreshed after the permutation, and not
     * permuted with the precisions, so that each of them is drawn from the
     * precision that occupies its slot, under the prior of that slot. */
    aux[0] = comp->refresh(&spec.comp[0], tau2[0], aux[0]);
    aux[1] = comp->refresh(&spec.comp[1], tau2[1], aux[1]);

    /* --- Allocation variables -------------------------------------------- */
    sd[0] = 1.0/sqrt(tau2[0]);
    sd[1] = 1.0/sqrt(tau2[1]);
    WBDrawAlloc(ry, alpha, n, mu, sd, z);

    /* --- Latent variables of the augmentation ---------------------------- */
    /* The regression lives on the padded grid, and so do the allocation
     * variables, which are gathered through the padding index. */
    for(t = 0; t < npad; t++)
      lat[t] = WBDrawLatentProbit(fit[t], z[idx[t]]);

    WaveDecP(lat, npad, 0, rwfilter, N, dstar, work);

    /* --- The scaling coefficient, under its diffuse prior ---------------- */
    theta[0] = dstar[0] + sigma*norm_rand();

    /* --- Level hyperparameters, inclusion indicators and coefficients ---- */
    /* The levels are visited from the coarsest to the finest one, which is the
     * order in which the transform lays them out. Each of them is complete
     * before the next one starts, and it is the previous state of the level
     * that its hyperparameters are drawn from. */
    for(j = 0; j < J; j++){

      m = 1 << j;
      lo = m;

      wprior->draw_level(&spec, dstar + lo, sigma, m, theta + lo, gam + lo,
                         loc + lo, wgt + lo, &pilev[j], &parlev[j]);
    }

    /* --- The mixture weights --------------------------------------------- */
    WaveRecP(theta, npad, 0, rwfilter, N, fit, work);

    for(t = 0; t < n; t++)
      alpha[t] = WBProbitMean(fit[t]);

    /* --- Storage ---------------------------------------------------------- */
    if(i > burn && (i - burn) % lag == 0){

      rmu[keep] = mu[0];
      rmu[keep + nchain] = mu[1];
      rtau2[keep] = tau2[0];
      rtau2[keep + nchain] = tau2[1];

      for(j = 0; j < J; j++){
        rpi[keep + nchain*j] = pilev[j];
        rpar[keep + nchain*j] = parlev[j];
      }

      for(t = 0; t < n; t++){
        ralpha[keep + nchain*t] = alpha[t];
        rzbar[t] += z[t];
      }

      for(t = 0; t < npad; t++)
        rincl[t] += wgt[t];

      keep++;
    }

    if(nverb > 0 && i % nverb == 0)
      Rprintf("  %d iterations of %d.\n", i, total);
  }

  PutRNGstate();

  for(t = 0; t < n; t++)
    rzbar[t] /= nchain;
  for(t = 0; t < npad; t++)
    rincl[t] /= nchain;

  PROTECT(nms = allocVector(STRSXP, 7));
  SET_STRING_ELT(nms, 0, mkChar("mu"));
  SET_STRING_ELT(nms, 1, mkChar("tau2"));
  SET_STRING_ELT(nms, 2, mkChar("alpha"));
  SET_STRING_ELT(nms, 3, mkChar("pi"));
  SET_STRING_ELT(nms, 4, mkChar("slab"));
  SET_STRING_ELT(nms, 5, mkChar("z"));
  SET_STRING_ELT(nms, 6, mkChar("inclusion"));
  setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(3);
  return out;
}
