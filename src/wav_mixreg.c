/**
 * @file wav_mixreg.c
 * @brief Wavelet-based estimators for mixture regression.
 * @details Implements the raw scaling-coefficient estimators (8) and (9) of
 *          Montoril, Pinheiro and Vidakovic (2019), the synthesis of a
 *          periodized scaling-function expansion, and the local linear
 *          smoother used for the regularization in the time domain.
 *
 *          Two design decisions drive the performance of this file. First,
 *          every routine accumulates directly into the periodized coefficient
 *          vector of length 2^J, so the n x 2^J design matrix built by the
 *          original MATLAB code never exists. Second, the integrals of
 *          estimator (9) are evaluated through the primitive of the scaling
 *          function rather than by quadrature at each observation, replacing
 *          O((b - a)/delta) evaluations of phi per observation by two table
 *          lookups per non-null coefficient.
 * @author Michel H. Montoril
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_mixreg.h"
#include "wav_utilities.h"
#include "phi_psi_vec.h"
#include "phi_psi_interp.h"
#include "utils.h"

//==============================================================================
// PRIMITIVE OF THE SCALING FUNCTION
//==============================================================================

/**
 * @brief Dense sample of phi and of its primitive over the compact support.
 *
 * @details The field @c f holds phi at the points q/G, q = 0, ..., nsamp - 1,
 *          covering the whole support [0, N - 1], and @c F holds
 *          Phi(q/G) = int_0^{q/G} phi. Both are indexed by the same q, so an
 *          evaluation is a division, an integer truncation and a fused
 *          correction inside the grid cell.
 */
typedef struct {
  double *f;    /**< phi sampled at q/G, length nsamp.                        */
  double *F;    /**< primitive of phi at the same points, length nsamp.       */
  int     G;    /**< number of sub-intervals per unit interval.               */
  int     nsamp;/**< number of samples, (N - 1)*G + 1.                        */
} PhiPrim;

/**
 * @brief Computes the scaling function at the integers of its support.
 *
 * @details At integer arguments the two-scale relation
 *          \f$\phi(t) = \sqrt 2 \sum_l h_l \phi(2t - l)\f$ closes on itself, so
 *          the interior values solve the eigenvalue problem
 *          \f$u_m = \sqrt 2 \sum_j h_{2m-j} u_j\f$ associated with the
 *          eigenvalue one. The dominant eigenvector is obtained by power
 *          iteration and fixed by the partition of unity, which normalizes the
 *          values to sum to one. The end points of the compact support vanish.
 *
 * @param[out] u      Values of phi at 0, 1, ..., N - 1 (length N).
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter size.
 */
static void MixPhiIntegers(double *u, const double *filter, int N)
{
  double acc, s;
  double *v, *vnew;
  int d, idx, it, j, m;

  d = N - 2;

  /* The Haar scaling function is the indicator of the unit interval. */
  if(d <= 0){
    u[0] = 1.0;
    u[N - 1] = 0.0;
    return;
  }

  v = (double *) R_alloc(d, sizeof(double));
  vnew = (double *) R_alloc(d, sizeof(double));

  for(j = 0; j < d; j++)
    v[j] = 1.0/d;

  for(it = 0; it < 1000; it++){

    for(m = 1; m <= d; m++){
      acc = 0.0;
      for(j = 1; j <= d; j++){
        idx = 2*m - j;
        if(idx >= 0 && idx < N)
          acc += filter[idx]*v[j - 1];
      }
      vnew[m - 1] = M_SQRT2*acc;
    }

    s = 0.0;
    for(j = 0; j < d; j++)
      s += vnew[j];

    if(s == 0.0)
      break;

    for(j = 0; j < d; j++)
      v[j] = vnew[j]/s;
  }

  u[0] = 0.0;
  for(j = 1; j <= d; j++)
    u[j] = v[j - 1];
  u[N - 1] = 0.0;
}

/**
 * @brief Samples the scaling function on a dyadic grid by the cascade algorithm.
 *
 * @details Starting from the values at the integers, each pass halves the grid
 *          spacing: the points already available are spread over the even
 *          positions and the new odd positions are filled by the two-scale
 *          relation, whose right-hand side only involves points of the previous
 *          grid. The whole sample therefore costs O(N^2 2^m), against the
 *          O(N^3 2^m prec) of one Daubechies-Lagarias run per grid point, and
 *          the values it produces are exact at every dyadic argument.
 *
 * @param[out] f      Values of phi at q/2^m, q = 0, ..., (N - 1)*2^m
 *                    (length (N - 1)*2^m + 1).
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter size.
 * @param[in]  m      Number of refinement passes.
 */
static void MixPhiCascade(double *f, const double *filter, int N, int m)
{
  double acc;
  int idx, l, lev, nc, nf, q, r, step;

  MixPhiIntegers(f, filter, N);

  for(lev = 0; lev < m; lev++){

    nc = (N - 1) << lev;
    nf = nc << 1;
    step = 1 << (lev + 1);

    /* Spread the coarse values downwards so that nothing is overwritten
     * before it has been read. */
    for(q = nc; q >= 0; q--)
      f[2*q] = f[q];

    /* The new points sit at odd indices, and every argument of the two-scale
     * relation lands on an even index of the same array. */
    for(r = 1; r <= nf; r += 2){
      acc = 0.0;
      for(l = 0; l < N; l++){
        idx = 2*r - l*step;
        if(idx >= 0 && idx <= nf)
          acc += filter[l]*f[idx];
      }
      f[r] = M_SQRT2*acc;
    }
  }
}

/**
 * @brief Builds the primitive of phi by cumulative Simpson quadrature.
 *
 * @details Odd indices are chained from a three-point Newton-Cotes start and
 *          even indices from Simpson's rule, so both sequences are accurate to
 *          O(h^4) while every grid point receives a value. Plain cumulative
 *          trapezoids would only be O(h^2) and would dominate the error budget
 *          of the estimator.
 *
 * @param[out] F     Primitive at the sample points, length nsamp.
 * @param[in]  f     Sample of phi, length nsamp (nsamp >= 3).
 * @param[in]  nsamp Number of samples.
 * @param[in]  h     Distance between consecutive sample points.
 */
static void MixPrimBuild(double *F, const double *f, int nsamp, double h)
{
  int q;

  F[0] = 0.0;
  F[1] = h*(5.0*f[0] + 8.0*f[1] - f[2])/12.0;

  for(q = 2; q < nsamp; q++)
    F[q] = F[q - 2] + h*(f[q - 2] + 4.0*f[q - 1] + f[q])/3.0;
}

/**
 * @brief Evaluates the primitive of phi at an arbitrary point.
 *
 * @details Outside the support the primitive is constant, equal to 0 on the
 *          left and to the total mass on the right. Inside a grid cell the
 *          remaining mass is obtained by integrating the linear interpolant of
 *          phi, which keeps the local error at O(h^3) without any extra
 *          storage.
 *
 * @param[in] pp Primitive structure.
 * @param[in] u  Evaluation point.
 * @return The value of the primitive at @p u.
 */
static double MixPrimEval(const PhiPrim *pp, double u)
{
  double s, v;
  int q;

  if(u <= 0.0)
    return 0.0;

  v = u * pp->G;
  q = (int) v;

  if(q >= pp->nsamp - 1)
    return pp->F[pp->nsamp - 1];

  s = v - q;

  return pp->F[q] + s*(pp->f[q] + 0.5*s*(pp->f[q + 1] - pp->f[q]))/pp->G;
}

//==============================================================================
// ACCUMULATORS FOR THE RAW ESTIMATORS
//==============================================================================

/**
 * @brief Accumulates one observation into the coefficients of estimator (8).
 *
 * @details Adds coefw * phi_{Jk}(x) to every periodized position k for which
 *          the scaling function is non-null at x. Only N - 1 positions are
 *          touched, whatever the resolution level.
 *
 * @param[in,out] coef  Coefficient vector of length p, indexed by k mod p.
 * @param[in]     p     2^J.
 * @param[in]     sqp   Square root of p.
 * @param[in]     phi   Values returned by PhiVec()/PhiVecInterp() at p*x.
 * @param[in]     px    The scaled point p*x.
 * @param[in]     N     Filter size.
 * @param[in]     coefw Multiplier of the observation.
 */
static void MixAccumPhi(double *coef, int p, double sqp, const double *phi,
                        double px, int N, double coefw)
{
  int j, k, lkmin;

  lkmin = (int) floor(px) - N + 2;

  for(j = 0; j < (N - 1); j++){
    k = lkmin + j;
    coef[mod(k, p)] += coefw * sqp * phi[j];
  }
}

/**
 * @brief Accumulates one observation into the coefficients of estimator (9).
 *
 * @details Uses the change of variables u = 2^J t - k to write the integral of
 *          phi_{Jk} over [a/p, b/p] as a difference of the primitive, so each
 *          non-null coefficient costs two lookups. The loop covers the k for
 *          which the support of phi_{Jk} meets the interval; the few boundary
 *          values of k that contribute exactly zero are harmless.
 *
 * @param[in,out] coef  Coefficient vector of length p, indexed by k mod p.
 * @param[in]     p     2^J.
 * @param[in]     sqp   Square root of p.
 * @param[in]     pp    Primitive structure of phi.
 * @param[in]     a     Scaled lower end point, p*lower.
 * @param[in]     b     Scaled upper end point, p*upper.
 * @param[in]     N     Filter size.
 * @param[in]     coefw Multiplier of the observation.
 */
static void MixAccumInt(double *coef, int p, double sqp, const PhiPrim *pp,
                        double a, double b, int N, double coefw)
{
  double val;
  int k, kA, kB;

  kA = (int) floor(a) - N + 2;
  kB = (int) floor(b);

  for(k = kA; k <= kB; k++){
    val = MixPrimEval(pp, b - k) - MixPrimEval(pp, a - k);
    coef[mod(k, p)] += coefw * val / sqp;
  }
}

//==============================================================================
// R INTERFACE
//==============================================================================

SEXP C_MixCoef(SEXP w, SEXP x, SEXP lower, SEXP upper, SEXP s, SEXP J,
               SEXP family, SEXP fs, SEXP prec, SEXP waveletfilter,
               SEXP what, SEXP intmethod, SEXP intpars, SEXP evaltab,
               SEXP ngrid)
{
  double a, b, delta, dl, hgrid, px, rs, sqp, wi;
  double *cint, *cphi, *phi, *prod, *rlower, *rtmp, *rupper, *rw, *rwfilter, *rx;
  int doint, dophi, G, i, imeth, leng, m, mgrid, minpts, n, N, nsamp, p, rJ,
      rprec, rwhat;
  const double *revaltab;
  PhiPrim pp;
  SEXP coefint, coefphi, out;

  n = LENGTH(x);
  rJ = INTEGER(J)[0];
  rprec = INTEGER(prec)[0];
  rwhat = INTEGER(what)[0];
  imeth = INTEGER(intmethod)[0];
  rs = REAL(s)[0];
  delta = REAL(intpars)[0];
  minpts = (int) REAL(intpars)[1];

  rw = REAL(w);
  rx = REAL(x);
  rlower = REAL(lower);
  rupper = REAL(upper);

  dophi = (rwhat == 1 || rwhat == 3);
  doint = (rwhat == 2 || rwhat == 3);

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  p = 1 << rJ;
  sqp = sqrt((double) p);

  revaltab = NULL;
  G = 0;
  if(evaltab != R_NilValue){
    revaltab = REAL(evaltab);
    G = ncols(evaltab) - 1;
    if(nrows(evaltab) != (N - 1))
      error("'wavelet.table' does not match the requested filter.");
  }

  /* The primitive is built by the cascade algorithm, so it is independent of
   * how phi is evaluated at the observation times and costs a negligible
   * fraction of one Daubechies-Lagarias run per grid point. */
  pp.f = NULL;
  pp.F = NULL;
  pp.G = 0;
  pp.nsamp = 0;
  if(doint && imeth == 0){
    mgrid = INTEGER(ngrid)[0];
    if(mgrid < 1 || mgrid > 24)
      error("'ngrid' must be a power of two between 2 and 2^24.");
    nsamp = ((N - 1) << mgrid) + 1;
    hgrid = 1.0 / (1 << mgrid);
    pp.f = (double *) R_alloc(nsamp, sizeof(double));
    pp.F = (double *) R_alloc(nsamp, sizeof(double));
    pp.G = 1 << mgrid;
    pp.nsamp = nsamp;
    MixPhiCascade(pp.f, rwfilter, N, mgrid);
    MixPrimBuild(pp.F, pp.f, nsamp, hgrid);
  }

  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  rtmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  PROTECT(coefphi = dophi ? allocVector(REALSXP, p) : R_NilValue);
  PROTECT(coefint = doint ? allocVector(REALSXP, p) : R_NilValue);

  cphi = dophi ? REAL(coefphi) : NULL;
  cint = doint ? REAL(coefint) : NULL;

  for(i = 0; i < p; i++){
    if(dophi)
      cphi[i] = 0.0;
    if(doint)
      cint[i] = 0.0;
  }

  for(i = 0; i < n; i++){

    if(!R_FINITE(rx[i]) || !R_FINITE(rw[i]))
      continue;

    wi = rw[i] / rs;

    /* --- Estimator (8): the rectangle rule of the paper --- */
    if(dophi){
      px = p * rx[i];

      if(revaltab)
        PhiVecInterp(phi, px, revaltab, N, G);
      else
        PhiVec(phi, px, rwfilter, N, rprec, prod, rtmp);

      MixAccumPhi(cphi, p, sqp, phi, px, N, wi*(rupper[i] - rlower[i]));
    }

    /* --- Estimator (9): the integrated scaling functions --- */
    if(doint){
      a = p * rlower[i];
      b = p * rupper[i];

      if(imeth == 0){
        MixAccumInt(cint, p, sqp, &pp, a, b, N, wi);
      }
      else{
        /* Literal trapezoidal rule of the reference implementation: 'leng'
         * equally spaced points spanning the interval, with half weights at
         * the end points. */
        leng = (int) ceil((rupper[i] - rlower[i]) / delta);
        if(leng < minpts)
          leng = minpts;
        if(leng < 2)
          leng = 2;

        dl = (rupper[i] - rlower[i]) / (leng - 1);

        for(m = 0; m < leng; m++){
          px = p * (rlower[i] + m*dl);

          if(revaltab)
            PhiVecInterp(phi, px, revaltab, N, G);
          else
            PhiVec(phi, px, rwfilter, N, rprec, prod, rtmp);

          MixAccumPhi(cint, p, sqp, phi, px, N,
                      wi * dl * ((m == 0 || m == leng - 1) ? 0.5 : 1.0));
        }
      }
    }

    R_CheckUserInterrupt();
  }

  PROTECT(out = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, coefphi);
  SET_VECTOR_ELT(out, 1, coefint);

  UNPROTECT(4);
  return out;
}

SEXP C_MixEval(SEXP coef, SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec,
               SEXP waveletfilter, SEXP phitab)
{
  double acc, px, sqp;
  double *phi, *prod, *rcoef, *rout, *rtmp, *rwfilter, *rx;
  int G, i, j, lkmin, M, N, p, rJ, rprec;
  const double *rphitab;
  SEXP out;

  M = LENGTH(x);
  rJ = INTEGER(J)[0];
  rprec = INTEGER(prec)[0];
  rcoef = REAL(coef);
  rx = REAL(x);

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  p = 1 << rJ;
  sqp = sqrt((double) p);

  if(LENGTH(coef) != p)
    error("'coef' must have length 2^J.");

  rphitab = NULL;
  G = 0;
  if(phitab != R_NilValue){
    rphitab = REAL(phitab);
    G = ncols(phitab) - 1;
    if(nrows(phitab) != (N - 1))
      error("'wavelet.table' does not match the requested filter.");
  }

  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  rtmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  PROTECT(out = allocVector(REALSXP, M));
  rout = REAL(out);

  for(i = 0; i < M; i++){

    if(!R_FINITE(rx[i])){
      rout[i] = NA_REAL;
      continue;
    }

    px = p * rx[i];

    if(rphitab)
      PhiVecInterp(phi, px, rphitab, N, G);
    else
      PhiVec(phi, px, rwfilter, N, rprec, prod, rtmp);

    lkmin = (int) floor(px) - N + 2;

    acc = 0.0;
    for(j = 0; j < (N - 1); j++)
      acc += rcoef[mod(lkmin + j, p)] * phi[j];

    rout[i] = sqp * acc;
  }

  UNPROTECT(2);
  return out;
}

SEXP C_LocLin(SEXP x, SEXP y, SEXP xout, SEXP bw)
{
  double d, den, h, num, s0, s1, s2, sy1, sy2, z;
  double *rout, *rx, *rxout, *ry;
  int i, m, M, n;
  SEXP out;

  n = LENGTH(x);
  M = LENGTH(xout);
  h = REAL(bw)[0];
  rx = REAL(x);
  ry = REAL(y);
  rxout = REAL(xout);

  PROTECT(out = allocVector(REALSXP, M));
  rout = REAL(out);

  for(m = 0; m < M; m++){

    s0 = s1 = s2 = sy1 = sy2 = 0.0;

    /* One pass accumulates the local moments of the kernel weights and the
     * two weighted sums of the responses needed by the closed form below. */
    for(i = 0; i < n; i++){
      d = rxout[m] - rx[i];
      z = exp(-0.5*(d/h)*(d/h));

      s0 += z;
      s1 += d*z;
      s2 += d*d*z;
      sy1 += z*ry[i];
      sy2 += d*z*ry[i];
    }

    num = s2*sy1 - s1*sy2;
    den = s2*s0 - s1*s1;

    rout[m] = (den == 0.0) ? NA_REAL : num/den;

    R_CheckUserInterrupt();
  }

  UNPROTECT(1);
  return out;
}
