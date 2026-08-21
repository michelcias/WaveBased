#ifndef WAVEBASED_WAV_POLYAGAMMA_H
#define WAVEBASED_WAV_POLYAGAMMA_H

/**
 * @file wav_polyagamma.h
 * @brief Exact sampler for the Polya-Gamma distribution PG(1, z).
 * @author Michel H. Montoril
 *
 * @details The Polya-Gamma distribution of Polson, Scott and Windle (2013) is
 *          what turns a logistic likelihood into a Gaussian one: for
 *          \f$\omega \sim PG(1, 0)\f$,
 *          \f[ \frac{e^{\eta\kappa}}{1 + e^{\eta}} \propto
 *              E\left[\exp\left(\kappa\eta -
 *              \frac{\omega\eta^2}{2}\right)\right], \f]
 *          so that, given \f$\omega\f$, the linear predictor enters through a
 *          normal kernel of mean \f$\kappa/\omega\f$ and variance
 *          \f$1/\omega\f$.
 *
 *          Only the unit shape is needed, which is the case Devroye's method
 *          covers exactly: PG(1, z) is a quarter of the Jacobi variable
 *          \f$J^*(1, z/2)\f$, and the latter is drawn by rejection from a
 *          mixture of a truncated exponential and a truncated inverse Gaussian
 *          proposal, accepted through the alternating series of its density.
 *          The draw is exact, and not a truncation of the infinite sum
 *          representation.
 *
 * @see Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference
 *      for logistic models using Polya-Gamma latent variables. Journal of the
 *      American Statistical Association, 108(504), 1339-1349.
 * @see Devroye, L. (2009). On exact simulation algorithms for some
 *      distributions related to Jacobi theta functions. Statistics and
 *      Probability Letters, 79(21), 2251-2259.
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Draws a single Polya-Gamma variable of unit shape.
 *
 * @details The caller is responsible for GetRNGstate() and PutRNGstate().
 *
 * @param[in] z Tilting parameter, the linear predictor. Its sign is immaterial.
 * @return The draw, a positive number.
 */
double WBDrawPG1(double z);

/**
 * @brief R entry point of the sampler, used by the unit tests.
 *
 * @param[in] z Real SEXP with one tilting parameter per draw.
 * @return Real SEXP of the same length with the draws.
 */
SEXP C_rpg1(SEXP z);

#endif /* WAVEBASED_WAV_POLYAGAMMA_H */
