#ifndef WAVEBASED_WAV_BAYESTHRESH_H
#define WAVEBASED_WAV_BAYESTHRESH_H

/**
 * @file wav_bayesthresh.h
 * @brief Bayesian wavelet thresholding (BayesThresh) for the R interface.
 * @author Michel H. Montoril
 *
 * @details Implements the Bayesian wavelet regularization of Abramovich,
 *          Sapatinas and Silverman (1998), which assigns to each detail
 *          coefficient the spike and slab prior
 *          \f[ \pi_j N(0, v_j^2) + (1 - \pi_j)\delta_0(\theta_{jk}), \qquad
 *              v_j^2 = 2^{-\alpha j} C_1, \quad
 *              \pi_j = \min(1, 2^{-\beta j} C_2), \f]
 *          and shrinks the coefficients towards the posterior median. The
 *          hyperparameters \f$C_1\f$ and \f$C_2\f$ are estimated from the data.
 *
 * @note The routines of this file were implemented based on the "BayesThresh"
 *       policy of the function threshold.wd() of the R package 'wavethresh'
 *       (G. P. Nason), which is the implementation used by Motta and Montoril
 *       (2026). The estimation of C1 relies on a one-dimensional minimizer that
 *       was, in turn, implemented based on Brent_fmin() of R (the engine of
 *       optimize()), itself a translation of the FMIN algorithm of Brent
 *       (1973). Both are used under the terms of the GNU General Public
 *       License.
 *
 * @see Abramovich, F., Sapatinas, T. and Silverman, B. W. (1998). Wavelet
 *      thresholding via a Bayesian approach. Journal of the Royal Statistical
 *      Society Series B, 60(4), 725-749.
 * @see Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
 *      Prentice-Hall, Englewood Cliffs.
 */

#include <R.h>
#include <Rinternals.h>

/** @brief Quantities selected by BayesThresh from the data. */
typedef struct {
  double sigma; /**< Noise level, the mad of the finest detail level. */
  double C1;    /**< Hyperparameter of the prior variances. */
  double C2;    /**< Hyperparameter of the prior mixing weights. */
} WBBayesThreshFit;

/**
 * @brief Minimizes a univariate function over an interval (internal).
 *
 * @details Implemented based on Brent_fmin() of R, the engine of optimize(),
 *          which is a C translation of the FMIN algorithm of Brent (1973). It
 *          combines golden section search with successive parabolic
 *          interpolation, and is reproduced here so that the estimation of the
 *          hyperparameter C1 inside the Gibbs sampler needs no callback into R.
 *
 * @param[in] ax   Lower end of the search interval.
 * @param[in] bx   Upper end of the search interval.
 * @param[in] f    Objective function, evaluated as f(x, info).
 * @param[in] info Opaque pointer forwarded to the objective function.
 * @param[in] tol  Desired accuracy. optimize() uses .Machine$double.eps^0.25.
 * @return The abscissa of the minimum found.
 */
double WBBrentFmin(double ax, double bx, double (*f)(double, void *),
                   void *info, double tol);

/**
 * @brief Applies BayesThresh to a vector of wavelet coefficients, in place.
 *
 * @details The coefficients must be laid out as
 *          \f$(c_{00}, d_{00}, d_{1,\cdot}, \ldots, d_{J-1,\cdot})\f$, that is,
 *          decomposed down to the level zero, so that the detail level j
 *          occupies the positions \f$2^j, \ldots, 2^{j+1} - 1\f$. The levels
 *          jlo up to J - 1 are thresholded, and the coarser ones are left
 *          untouched, as in the default of threshold.wd().
 *
 *          The noise level is the mad of the finest detail level. When C1 is
 *          not supplied, it is estimated by the marginal maximum likelihood
 *          described in Abramovich, Sapatinas and Silverman (1998), restricted
 *          to the coefficients that survive a universal hard threshold, and the
 *          minimization is carried out by WBBrentFmin() over
 *          \f$[0, 50\sqrt{C_{1,start}}]\f$. When C2 is not supplied, it follows
 *          in closed form.
 *
 * @param[in,out] wc      Wavelet coefficients of length n, thresholded in place.
 * @param[in]     n       Length of wc, a power of two.
 * @param[in]     jlo     Coarsest detail level to threshold, 0 <= jlo < log2(n).
 * @param[in]     alpha   Prior parameter alpha, non-negative.
 * @param[in]     beta    Prior parameter beta, non-negative.
 * @param[in]     C1      Hyperparameter C1, or NA_REAL to estimate it.
 * @param[in]     C2      Hyperparameter C2, or NA_REAL to estimate it.
 * @param[in]     C1start Starting scale of the search interval of C1.
 * @param[out]    work    Scratch buffer with room for n doubles.
 * @param[out]    fit     Quantities selected from the data. May be NULL.
 */
void WBBayesThresh(double *wc, int n, int jlo, double alpha, double beta,
                   double C1, double C2, double C1start, double *work,
                   WBBayesThreshFit *fit);

/**
 * @brief Bayesian wavelet regularization of a signal (R interface).
 *
 * @param[in] x             Real SEXP: signal of length n, a power of two.
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets,
 *                          2 = Symmlets, 3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] j0            Integer SEXP: coarsest detail level to threshold.
 * @param[in] hyper         Real SEXP of length 5 holding alpha, beta, C1, C2
 *                          and C1.start. C1 and C2 may be NA.
 * @return SEXP list with the regularized signal, the thresholded wavelet
 *         coefficients, and the values of sigma, C1 and C2.
 */
SEXP C_BayesThresh(SEXP x, SEXP family, SEXP fs, SEXP waveletfilter, SEXP j0,
                   SEXP hyper);

#endif /* WAVEBASED_WAV_BAYESTHRESH_H */
