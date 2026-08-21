/**
 * @file wav_polyagamma.c
 * @brief Exact sampler for the Polya-Gamma distribution PG(1, z).
 * @author Michel H. Montoril
 *
 * @details Devroye's method, as described in Section 4 of Polson, Scott and
 *          Windle (2013). The variable is written as
 *          \f$\omega = J^*(1, z/2)/4\f$, and the Jacobi variable is drawn by
 *          rejection: the density is split at the point 0.64, below which a
 *          truncated inverse Gaussian proposal dominates and above which a
 *          truncated exponential one does, and the acceptance is decided by the
 *          alternating series of the density, whose partial sums bracket it.
 *
 *          The series is the reason the draw is exact: a partial sum of an
 *          alternating series with decreasing terms lies on one side of the
 *          limit, so the comparison against a uniform variate is decided after
 *          finitely many terms, without the tail ever being evaluated.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "wav_polyagamma.h"

/** @name Constants of the split of the density. */
/**@{*/
#define WB_PG_TRUNC       0.64   /**< Point where the two proposals meet.   */
#define WB_PG_TRUNC_RECIP 1.5625 /**< Its reciprocal, 1/0.64.               */
#define WB_PG_MAXIT       1000   /**< Cap on the rejection loop.            */
/**@}*/

/**
 * @brief Term n of the alternating series of the Jacobi density.
 *
 * @details The density of \f$J^*(1, 0)\f$ has two series representations, one
 *          converging quickly in the right tail and the other in the left one,
 *          and the split at WB_PG_TRUNC is where each of them is used.
 *
 * @param[in] n Index of the term.
 * @param[in] x Point at which the density is evaluated.
 * @return The term, a positive number.
 */
static double WBPGCoef(int n, double x){

  double k = (n + 0.5)*M_PI;

  if(x > WB_PG_TRUNC)
    return k*exp(-0.5*k*k*x);

  if(x > 0.0)
    return exp(-1.5*(log(0.5*M_PI) + log(x)) + log(k)
               - 2.0*(n + 0.5)*(n + 0.5)/x);

  return 0.0;
}

/**
 * @brief Probability that the proposal is the truncated exponential one.
 *
 * @details The ratio of the two masses is evaluated on the logarithmic scale,
 *          so that a large tilting parameter does not overflow.
 *
 * @param[in] z Half the tilting parameter, non-negative.
 * @return The mixing probability.
 */
static double WBPGMassTexpon(double z){

  double t = WB_PG_TRUNC;
  double fz = 0.125*M_PI*M_PI + 0.5*z*z;
  double b = sqrt(1.0/t)*(t*z - 1.0);
  double a = -sqrt(1.0/t)*(t*z + 1.0);
  double x0 = log(fz) + fz*t;
  double xb = x0 - z + pnorm(b, 0.0, 1.0, 1, 1);
  double xa = x0 + z + pnorm(a, 0.0, 1.0, 1, 1);

  return 1.0/(1.0 + 4.0/M_PI*(exp(xb) + exp(xa)));
}

/**
 * @brief Draws an inverse Gaussian variable truncated to (0, WB_PG_TRUNC).
 *
 * @details The mean of the untruncated variable is 1/z. When it lies to the
 *          right of the truncation point the draw comes from Devroye's
 *          rejection scheme built on the exponential distribution, and
 *          otherwise from the Michael, Schucany and Haas transformation,
 *          repeated until it lands inside the interval.
 *
 * @param[in] z Half the tilting parameter, non-negative.
 * @return The draw, in (0, WB_PG_TRUNC).
 */
static double WBPGRTIGauss(double z){

  double t = WB_PG_TRUNC, x = t + 1.0, alpha, e1, e2, mu, y, hmu, muy;

  if(z < WB_PG_TRUNC_RECIP){

    alpha = 0.0;

    while(unif_rand() > alpha){

      e1 = exp_rand();
      e2 = exp_rand();

      while(e1*e1 > 2.0*e2/t){
        e1 = exp_rand();
        e2 = exp_rand();
      }

      x = 1.0 + e1*t;
      x = t/(x*x);
      alpha = exp(-0.5*z*z*x);
    }
  }
  else{

    mu = 1.0/z;

    while(x > t){

      y = norm_rand();
      y = y*y;
      hmu = 0.5*mu;
      muy = mu*y;
      x = mu + hmu*muy - hmu*sqrt(4.0*muy + muy*muy);

      if(unif_rand() > mu/(mu + x))
        x = mu*mu/x;
    }
  }

  return x;
}

double WBDrawPG1(double z){

  double half = 0.5*fabs(z);
  double fz = 0.125*M_PI*M_PI + 0.5*half*half;
  double x, s, u;
  int n, go, it = 0;

  while(it++ < WB_PG_MAXIT){

    if(unif_rand() < WBPGMassTexpon(half))
      x = WB_PG_TRUNC + exp_rand()/fz;
    else
      x = WBPGRTIGauss(half);

    /* The partial sums of the alternating series bracket the density, so the
     * odd ones decide acceptance and the even ones decide rejection. */
    s = WBPGCoef(0, x);
    u = unif_rand()*s;
    n = 0;
    go = 1;

    while(go){

      n++;

      if(n % 2){
        s -= WBPGCoef(n, x);
        if(u <= s)
          return 0.25*x;
      }
      else{
        s += WBPGCoef(n, x);
        if(u > s)
          go = 0;
      }
    }
  }

  error("The Polya-Gamma sampler did not accept a draw.");
  return 0.0;
}

SEXP C_rpg1(SEXP z){

  int i, n = length(z);
  double *rz = REAL(z), *out;
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, n));
  out = REAL(ans);

  GetRNGstate();

  for(i = 0; i < n; i++)
    out[i] = WBDrawPG1(rz[i]);

  PutRNGstate();
  UNPROTECT(1);

  return ans;
}
