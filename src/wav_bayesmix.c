/**
 * @file wav_bayesmix.c
 * @brief Gibbs sampler for the wavelet-based Bayesian mixture regression.
 * @author Michel H. Montoril
 *
 * @details Implements Algorithm 1 of Motta and Montoril (2026). The sampler was
 *          implemented based on the reference R code of the article, available
 *          at https://github.com/flaviamotta/Bayesian-wavelet-mixture-regression,
 *          with the whole sweep moved into C.
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#include "wav_bayesmix.h"
#include "wav_bayesthresh.h"
#include "wav_transform.h"
#include "wav_utilities.h"

void WBDrawComponent(const double *y, const int *z, int n, int grp,
                     double b0, double B0, double v0, double V0,
                     double *mu, double *tau2){

  int t, nk = 0;
  double sy = 0.0, ss = 0.0, Bk, bk, m;

  for(t = 0; t < n; t++)
    if(z[t] == grp){
      nk++;
      sy += y[t];
    }

  Bk = 1.0/(nk*(*tau2) + 1.0/B0);
  bk = Bk*(sy*(*tau2) + b0/B0);
  m = bk + sqrt(Bk)*norm_rand();

  for(t = 0; t < n; t++)
    if(z[t] == grp)
      ss += (y[t] - m)*(y[t] - m);

  *mu = m;
  *tau2 = rgamma(v0 + 0.5*nk, 1.0/(V0 + 0.5*ss));
}

void WBDrawAlloc(const double *y, const double *alpha, int n,
                 const double *mu, const double *sd, int *z){

  int t;
  double l0, l1, u;

  for(t = 0; t < n; t++){
    u = unif_rand();
    if(alpha[t] <= 0.0)
      z[t] = 0;
    else if(alpha[t] >= 1.0)
      z[t] = 1;
    else{
      l0 = log1p(-alpha[t]) + dnorm(y[t], mu[0], sd[0], 1);
      l1 = log(alpha[t]) + dnorm(y[t], mu[1], sd[1], 1);
      z[t] = (u < 1.0/(1.0 + exp(l0 - l1))) ? 1 : 0;
    }
  }
}

SEXP C_BayesMixReg(SEXP y, SEXP padidx, SEXP family, SEXP fs,
                   SEXP waveletfilter, SEXP j0, SEXP hyper, SEXP prior,
                   SEXP init, SEXP mcmc){

  int i, t, n = length(y), npad = length(padidx), N, keep, J = 0;
  int jlo = INTEGER(j0)[0], *idx = INTEGER(padidx), *z;
  int nchain = INTEGER(mcmc)[0], burn = INTEGER(mcmc)[1];
  int lag = INTEGER(mcmc)[2], nverb = INTEGER(mcmc)[3];
  int total = burn + lag*nchain, ndeg = 0;
  double *ry = REAL(y), *rhyper = REAL(hyper), *rprior = REAL(prior);
  double *rwfilter, *work, *w, *wc, *fit, *alpha;
  double *rmu, *rtau2, *ralpha, *rsigma, *rC1, *rC2, *rzbar;
  double mu[2], tau2[2], delta, sd[2], tmp;
  WBBayesThreshFit btfit;
  SEXP out, nms;

  while((1 << (J + 1)) <= npad) J++;

  if((1 << J) != npad)
    error("The padded sample size must be a power of two.");
  if(npad < n)
    error("The padded sample size must not be smaller than the sample size.");
  if(nchain < 1 || burn < 0 || lag < 1)
    error("Invalid MCMC settings.");
  /* Checked here, and not only inside WBBayesThresh, so that no error can be
   * raised after the state of the random number generator has been read. */
  if(jlo < 0 || jlo >= J)
    error("The coarsest level to threshold must lie between 0 and log2(n) - 1.");

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  /* --- Working space, allocated once for the whole run ------------------- */
  work  = (double *) R_alloc(WB_WAVE_WORK(npad), sizeof(double));
  w     = (double *) R_alloc(npad, sizeof(double));
  wc    = (double *) R_alloc(npad, sizeof(double));
  fit   = (double *) R_alloc(npad, sizeof(double));
  alpha = (double *) R_alloc(n, sizeof(double));
  z     = (int *)    R_alloc(n, sizeof(int));

  /* --- The retained draws ------------------------------------------------ */
  PROTECT(out = allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, allocMatrix(REALSXP, nchain, 2));  /* mu     */
  SET_VECTOR_ELT(out, 1, allocMatrix(REALSXP, nchain, 2));  /* tau2   */
  SET_VECTOR_ELT(out, 2, allocMatrix(REALSXP, nchain, n));  /* alpha  */
  SET_VECTOR_ELT(out, 3, allocVector(REALSXP, nchain));     /* sigma  */
  SET_VECTOR_ELT(out, 4, allocVector(REALSXP, nchain));     /* C1     */
  SET_VECTOR_ELT(out, 5, allocVector(REALSXP, nchain));     /* C2     */
  SET_VECTOR_ELT(out, 6, allocVector(REALSXP, n));          /* z mean */

  rmu    = REAL(VECTOR_ELT(out, 0));
  rtau2  = REAL(VECTOR_ELT(out, 1));
  ralpha = REAL(VECTOR_ELT(out, 2));
  rsigma = REAL(VECTOR_ELT(out, 3));
  rC1    = REAL(VECTOR_ELT(out, 4));
  rC2    = REAL(VECTOR_ELT(out, 5));
  rzbar  = REAL(VECTOR_ELT(out, 6));

  for(t = 0; t < n; t++)
    rzbar[t] = 0.0;

  /* --- Initial values ---------------------------------------------------- */
  mu[0] = REAL(init)[0];
  mu[1] = REAL(init)[1];
  tau2[0] = REAL(init)[2];
  tau2[1] = REAL(init)[3];

  for(t = 0; t < n; t++)
    alpha[t] = 0.5;
  for(t = 0; t < n; t++)
    z[t] = 0;

  GetRNGstate();

  /* The allocation variables start from the weights implied by the initial
   * component parameters, as in the reference implementation. */
  delta = mu[1] - mu[0];
  if(R_FINITE(delta) && delta != 0.0){
    for(t = 0; t < npad; t++)
      w[t] = (ry[idx[t]] - mu[0])/delta;
    WaveDecP(w, npad, 0, rwfilter, N, wc, work);
    WBBayesThresh(wc, npad, jlo, rhyper[0], rhyper[1], rhyper[2], rhyper[3],
                  rhyper[4], work, &btfit);
    WaveRecP(wc, npad, 0, rwfilter, N, fit, work);
    for(t = 0; t < n; t++)
      alpha[t] = fmin2(1.0, fmax2(0.0, fit[t]));
  }

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
    WBDrawComponent(ry, z, n, 0, rprior[0], rprior[1], rprior[2], rprior[3],
                    &mu[0], &tau2[0]);
    WBDrawComponent(ry, z, n, 1, rprior[4], rprior[5], rprior[6], rprior[7],
                    &mu[1], &tau2[1]);

    if(mu[1] < mu[0]){
      tmp = mu[0];    mu[0] = mu[1];       mu[1] = tmp;
      tmp = tau2[0];  tau2[0] = tau2[1];   tau2[1] = tmp;
    }

    /* --- Mixture weights, by Bayesian wavelet regularization ------------- */
    delta = mu[1] - mu[0];

    if(R_FINITE(delta) && delta != 0.0){
      for(t = 0; t < npad; t++)
        w[t] = (ry[idx[t]] - mu[0])/delta;

      WaveDecP(w, npad, 0, rwfilter, N, wc, work);
      WBBayesThresh(wc, npad, jlo, rhyper[0], rhyper[1], rhyper[2], rhyper[3],
                    rhyper[4], work, &btfit);
      WaveRecP(wc, npad, 0, rwfilter, N, fit, work);

      /* The weights are probabilities, so the reconstruction is truncated. */
      for(t = 0; t < n; t++)
        alpha[t] = fmin2(1.0, fmax2(0.0, fit[t]));
    }
    else{
      /* Coincident component means leave the transformation undefined. The
       * previous weights are kept and the event is reported at the end. */
      ndeg++;
      btfit.sigma = NA_REAL;
      btfit.C1 = NA_REAL;
      btfit.C2 = NA_REAL;
    }

    /* --- Allocation variables -------------------------------------------- */
    sd[0] = 1.0/sqrt(tau2[0]);
    sd[1] = 1.0/sqrt(tau2[1]);
    WBDrawAlloc(ry, alpha, n, mu, sd, z);

    /* --- Storage ---------------------------------------------------------- */
    if(i > burn && (i - burn) % lag == 0){
      rmu[keep] = mu[0];
      rmu[keep + nchain] = mu[1];
      rtau2[keep] = tau2[0];
      rtau2[keep + nchain] = tau2[1];
      rsigma[keep] = btfit.sigma;
      rC1[keep] = btfit.C1;
      rC2[keep] = btfit.C2;

      for(t = 0; t < n; t++){
        ralpha[keep + nchain*t] = alpha[t];
        rzbar[t] += z[t];
      }

      keep++;
    }

    if(nverb > 0 && i % nverb == 0)
      Rprintf("  %d iterations of %d.\n", i, total);
  }

  PutRNGstate();

  for(t = 0; t < n; t++)
    rzbar[t] /= nchain;

  if(ndeg > 0)
    warning("The component means coincided in %d of the %d sweeps, whose mixture weights were kept unchanged.",
            ndeg, total);

  PROTECT(nms = allocVector(STRSXP, 7));
  SET_STRING_ELT(nms, 0, mkChar("mu"));
  SET_STRING_ELT(nms, 1, mkChar("tau2"));
  SET_STRING_ELT(nms, 2, mkChar("alpha"));
  SET_STRING_ELT(nms, 3, mkChar("sigma"));
  SET_STRING_ELT(nms, 4, mkChar("C1"));
  SET_STRING_ELT(nms, 5, mkChar("C2"));
  SET_STRING_ELT(nms, 6, mkChar("z"));
  setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(3);
  return out;
}
