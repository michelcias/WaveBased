/**
 * @file wav_bayesthresh.c
 * @brief Bayesian wavelet thresholding (BayesThresh) for the R interface.
 * @author Michel H. Montoril
 *
 * @note Implemented based on the "BayesThresh" policy of the function
 *       threshold.wd() of the R package 'wavethresh' (G. P. Nason), which
 *       implements Abramovich, Sapatinas and Silverman (1998) and is the
 *       implementation used by Motta and Montoril (2026). The univariate
 *       minimizer was implemented based on Brent_fmin() of R, the engine of
 *       optimize(), itself a translation of the FMIN algorithm of Brent (1973).
 *       Both are used under the terms of the GNU General Public License.
 */

#include <float.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "wav_bayesthresh.h"
#include "wav_transform.h"
#include "wav_utilities.h"
#include "utils.h"

/** @brief Largest number of resolution levels a signal of int length may have. */
#define WB_MAXLEV 32

/** @brief Coefficients below this magnitude are not counted as signal. */
#define WB_BT_ZERO 1e-10

//==============================================================================
// UNIVARIATE MINIMIZATION
//==============================================================================

double WBBrentFmin(double ax, double bx, double (*f)(double, void *),
                   void *info, double tol)
{
  /*  c is the squared inverse of the golden ratio */
  const double c = (3. - sqrt(5.)) * .5;

  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  /*  eps is approximately the square root of the relative machine precision. */
  eps = DBL_EPSILON;
  tol1 = eps + 1.;  /* the smallest 1.000... > 1 */
  eps = sqrt(eps);

  a = ax;
  b = bx;
  v = a + c * (b - a);
  w = v;
  x = v;

  d = 0.;
  e = 0.;
  fx = (*f)(x, info);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  for(;;){
    xm = (a + b) * .5;
    tol1 = eps * fabs(x) + tol3;
    t2 = tol1 * 2.;

    /* check stopping criterion */
    if(fabs(x - xm) <= t2 - (b - a) * .5) break;
    p = 0.;
    q = 0.;
    r = 0.;
    if(fabs(e) > tol1){ /* fit parabola */
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if(q > 0.) p = -p; else q = -q;
      r = e;
      e = d;
    }

    if(fabs(p) >= fabs(q * .5 * r) ||
       p <= q * (a - x) || p >= q * (b - x)){ /* a golden-section step */
      if(x < xm) e = b - x; else e = a - x;
      d = c * e;
    }
    else{ /* a parabolic-interpolation step */
      d = p / q;
      u = x + d;

      /* f must not be evaluated too close to ax or bx */
      if(u - a < t2 || b - u < t2){
        d = tol1;
        if(x >= xm) d = -d;
      }
    }

    /* f must not be evaluated too close to x */
    if(fabs(d) >= tol1)
      u = x + d;
    else if(d > 0.)
      u = x + tol1;
    else
      u = x - tol1;

    fu = (*f)(u, info);

    /* update a, b, v, w, and x */
    if(fu <= fx){
      if(u < x) b = x; else a = x;
      v = w;   w = x;   x = u;
      fv = fw; fw = fx; fx = fu;
    }
    else{
      if(u < x) a = u; else b = u;
      if(fu <= fw || w == x){
        v = w; fv = fw;
        w = u; fw = fu;
      }
      else if(fu <= fv || v == x || v == w){
        v = u; fv = fu;
      }
    }
  }

  return x;
}

//==============================================================================
// BAYESIAN THRESHOLDING
//==============================================================================

/** @brief Data needed by the marginal likelihood of the hyperparameter C1. */
typedef struct {
  const int    *nsignal; /**< Number of surviving coefficients per level. */
  const double *sum2;    /**< Sum of their squares per level. */
  const double *v;       /**< The factors 2^(-alpha j) per level. */
  int           nthresh; /**< Number of thresholded levels. */
  double        sigma2;  /**< Squared noise level. */
  double        zuniv;   /**< -sigma*sqrt(2 log n), the universal quantile. */
} WBBTLik;

/**
 * @brief Minus twice the marginal log-likelihood of C, up to a constant.
 *
 * @details This is the objective minimized by threshold.wd() to estimate C1,
 *          restricted to the coefficients that survive a universal hard
 *          threshold. The normal tail is evaluated on the log scale, which is
 *          equivalent to, and safer than, taking the logarithm afterwards.
 *
 * @param[in] C    Candidate value, so that the prior variances are C^2 v_j.
 * @param[in] info Pointer to a WBBTLik structure.
 * @return The value of the objective at C.
 */
static double WBBayesThreshLik(double C, void *info){

  const WBBTLik *o = (const WBBTLik *) info;
  double ans = 0.0, s2v;
  int i;

  for(i = 0; i < o->nthresh; i++){
    s2v = o->sigma2 + C*C*o->v[i];
    ans += o->nsignal[i]*(log(s2v)
                          - 2.0*pnorm(o->zuniv/sqrt(s2v), 0.0, 1.0, 1, 1))
           + o->sum2[i]/s2v;
  }

  return ans;
}

/**
 * @brief Sets to zero every detail coefficient of the thresholded levels.
 *
 * @param[in,out] wc  Wavelet coefficients of length n.
 * @param[in]     n   Length of wc.
 * @param[in]     jlo Coarsest thresholded level.
 */
static void WBNullLevels(double *wc, int n, int jlo){

  int i;

  for(i = 1 << jlo; i < n; i++)
    wc[i] = 0.0;
}

void WBBayesThresh(double *wc, int n, int jlo, double alpha, double beta,
                   double C1, double C2, double C1start, double *work,
                   WBBayesThreshFit *fit){

  int i, j, k, m, J = 0, nthresh, ntot, nsig = 0;
  int nsignal[WB_MAXLEV];
  double v[WB_MAXLEV], sum2[WB_MAXLEV], tau2[WB_MAXLEV];
  double sigma, sigma2, zuniv, thr, d, ad, ww, z, rat, pr, psum;
  WBBTLik lik;

  while((1 << (J + 1)) <= n) J++;   /* J = log2(n) */

  if(fit != NULL){
    fit->sigma = NA_REAL;
    fit->C1 = NA_REAL;
    fit->C2 = NA_REAL;
  }

  if(jlo < 0 || jlo >= J)
    error("The coarsest level to threshold must lie between 0 and log2(n) - 1.");

  nthresh = J - jlo;

  /* --- The noise level, from the finest detail level --------------------- */
  sigma = WBMad(wc + n/2, n/2, work);

  if(!R_FINITE(sigma) || sigma <= 0.0)
    /* A degenerate noise level leaves the posterior median undefined. The
     * coefficients are kept as they are, and the caller is told about it
     * through the NA values already stored in fit. */
    return;

  sigma2 = sigma*sigma;
  zuniv = -sigma*sqrt(2.0*log((double) n));

  for(i = 0; i < nthresh; i++){
    v[i] = R_pow(2.0, -alpha*(jlo + i));
    nsignal[i] = 0;
    sum2[i] = 0.0;
  }

  /* --- Universal hard thresholding, which feeds C1 and C2 ---------------- */
  /* The thresholded levels are contiguous in the coefficient vector, so the
   * pooled detail coefficients are the tail wc[2^jlo], ..., wc[n - 1]. */
  ntot = n - (1 << jlo);
  thr = sqrt(2.0*log((double) ntot))*WBMad(wc + (1 << jlo), ntot, work);

  for(i = 0; i < nthresh; i++){
    m = 1 << (jlo + i);
    for(k = 0; k < m; k++){
      ad = fabs(wc[m + k]);
      if(ad > thr){
        sum2[i] += wc[m + k]*wc[m + k];
        if(ad > WB_BT_ZERO)
          nsignal[i]++;
      }
    }
    if(nsignal[i] == 0)
      sum2[i] = 0.0;
    nsig += nsignal[i];
  }

  /* --- The hyperparameters of the prior ---------------------------------- */
  if(ISNAN(C1)){
    if(nsig == 0){
      /* No coefficient survives the universal threshold: everything is noise. */
      WBNullLevels(wc, n, jlo);
      if(fit != NULL){
        fit->sigma = sigma;
        fit->C1 = 0.0;
        fit->C2 = 0.0;
      }
      return;
    }

    lik.nsignal = nsignal;
    lik.sum2 = sum2;
    lik.v = v;
    lik.nthresh = nthresh;
    lik.sigma2 = sigma2;
    lik.zuniv = zuniv;

    C1 = WBBrentFmin(0.0, 50.0*sqrt(C1start), WBBayesThreshLik, &lik,
                     R_pow(DBL_EPSILON, 0.25));
    C1 = C1*C1;
  }

  if(C1 < 0.0)
    error("The hyperparameter 'C1' must be non-negative.");

  for(i = 0; i < nthresh; i++)
    tau2[i] = C1*v[i];

  if(ISNAN(C2)){
    psum = 0.0;
    for(i = 0; i < nthresh; i++)
      psum += nsignal[i]/(2.0*pnorm(zuniv/sqrt(sigma2 + tau2[i]), 0.0, 1.0,
                                    1, 0));

    C2 = (beta == 1.0) ?
           psum/J :
           (1.0 - R_pow(2.0, 1.0 - beta))/
             (1.0 - R_pow(2.0, (1.0 - beta)*J))*psum;
  }

  if(C2 < 0.0)
    error("The hyperparameter 'C2' must be non-negative.");

  if(fit != NULL){
    fit->sigma = sigma;
    fit->C1 = C1;
    fit->C2 = C2;
  }

  /* --- Shrinkage towards the posterior median ---------------------------- */
  if(C1 == 0.0 || C2 == 0.0){
    WBNullLevels(wc, n, jlo);
    return;
  }

  for(i = 0; i < nthresh; i++){

    j = jlo + i;
    m = 1 << j;
    pr = fmin2(1.0, C2*R_pow(2.0, -beta*j));
    rat = tau2[i]/(sigma2 + tau2[i]);

    for(k = 0; k < m; k++){
      d = wc[m + k];
      ad = fabs(d);
      ww = (1.0 - pr)/pr/sqrt(sigma2*rat/tau2[i])*exp(-rat*d*d/(2.0*sigma2));
      z = 0.5*(1.0 + fmin2(ww, 1.0));
      /* z = 1 gives an infinite quantile, and the coefficient is killed. */
      ad = rat*ad - sigma*sqrt(rat)*qnorm(z, 0.0, 1.0, 1, 0);
      wc[m + k] = (d < 0 ? -1.0 : 1.0)*fmax2(0.0, ad);
    }
  }
}

//==============================================================================
// R INTERFACE
//==============================================================================

SEXP C_BayesThresh(SEXP x, SEXP family, SEXP fs, SEXP waveletfilter, SEXP j0,
                   SEXP hyper){

  int n = length(x), N, jlo = INTEGER(j0)[0], J = 0;
  double *rhyper = REAL(hyper), *work, *rwfilter;
  WBBayesThreshFit fit;
  SEXP out, coef, sig, nms;

  while((1 << (J + 1)) <= n) J++;

  if((1 << J) != n)
    error("The sample size must be a power of two.");
  if(J < 1)
    error("The sample size must be at least two.");

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  work = (double *) R_alloc(WB_WAVE_WORK(n), sizeof(double));

  PROTECT(coef = allocVector(REALSXP, n));
  WaveDecP(REAL(x), n, 0, rwfilter, N, REAL(coef), work);

  WBBayesThresh(REAL(coef), n, jlo, rhyper[0], rhyper[1], rhyper[2], rhyper[3],
                rhyper[4], work, &fit);

  PROTECT(out = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, allocVector(REALSXP, n));
  WaveRecP(REAL(coef), n, 0, rwfilter, N, REAL(VECTOR_ELT(out, 0)), work);
  SET_VECTOR_ELT(out, 1, coef);

  PROTECT(sig = allocVector(REALSXP, 3));
  REAL(sig)[0] = fit.sigma;
  REAL(sig)[1] = fit.C1;
  REAL(sig)[2] = fit.C2;
  PROTECT(nms = allocVector(STRSXP, 3));
  SET_STRING_ELT(nms, 0, mkChar("sigma"));
  SET_STRING_ELT(nms, 1, mkChar("C1"));
  SET_STRING_ELT(nms, 2, mkChar("C2"));
  setAttrib(sig, R_NamesSymbol, nms);
  SET_VECTOR_ELT(out, 2, sig);

  UNPROTECT(5);
  return out;
}
