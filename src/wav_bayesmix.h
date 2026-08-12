#ifndef WAVEBASED_WAV_BAYESMIX_H
#define WAVEBASED_WAV_BAYESMIX_H

/**
 * @file wav_bayesmix.h
 * @brief Gibbs sampler for the wavelet-based Bayesian mixture regression.
 * @author Michel H. Montoril
 *
 * @details Implements Algorithm 1 of Motta and Montoril (2026) for the dynamic
 *          two-component Gaussian mixture
 *          \f[ y_t = (1 - z_t)x_{1t} + z_t x_{2t}, \qquad
 *              x_{kt} \sim N(\mu_k, \tau_k^{-2}), \qquad
 *              z_t \sim \mathrm{Bern}(\alpha_t), \f]
 *          whose mixture weights \f$\alpha_t\f$ vary with the data index and
 *          are estimated by Bayesian wavelet regularization of
 *          \f$m_t = (y_t - \mu_1)/(\mu_2 - \mu_1)\f$.
 *
 *          The whole sweep runs in C: the transform uses the pyramid algorithm
 *          of WaveDecP()/WaveRecP(), the regularization uses WBBayesThresh(),
 *          and the component parameters are drawn from the conjugate full
 *          conditionals (3) and (4) of the article. No memory is allocated
 *          inside the loop.
 *
 * @see Motta, F. C. and Montoril, M. H. (2026). A Bayesian estimation approach
 *      for the wavelet-based mixture regression. Communications in Statistics -
 *      Simulation and Computation, 55(6), 2426-2434.
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Runs the Gibbs sampler of the dynamic Gaussian mixture model.
 *
 * @details Each sweep draws \f$\mu_k\f$ and \f$\tau_k^2\f$ from their conjugate
 *          full conditionals, permutes the labels whenever
 *          \f$\mu_2 < \mu_1\f$, rebuilds the mixture weights by Bayesian
 *          wavelet thresholding of the transformed responses, and draws the
 *          allocation variables. Only the draws after the burn-in period, taken
 *          every 'lag' sweeps, are stored.
 *
 * @param[in] y             Real SEXP: observed mixture of length n.
 * @param[in] padidx        Integer SEXP of length npad (a power of two) with
 *                          zero-based positions of y, describing how the sample
 *                          is extended to a dyadic length. The first n entries
 *                          must be 0, 1, ..., n - 1.
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets,
 *                          2 = Symmlets, 3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] j0            Integer SEXP: coarsest detail level to threshold.
 * @param[in] hyper         Real SEXP of length 5: alpha, beta, C1, C2 and
 *                          C1.start of the BayesThresh prior. C1 and C2 may
 *                          be NA, in which case they are estimated.
 * @param[in] prior         Real SEXP of length 8 holding, for each component,
 *                          the mean and the variance of the normal prior of
 *                          mu, and the shape and the rate of the gamma prior of
 *                          tau^2, in the order b01, B01, v01, V01, b02, B02,
 *                          v02, V02.
 * @param[in] init          Real SEXP of length 4: initial values of mu1, mu2,
 *                          tau1^2 and tau2^2.
 * @param[in] mcmc          Integer SEXP of length 4: number of draws to keep,
 *                          burn-in, lag, and the reporting interval (zero for
 *                          a silent run).
 *
 * @return SEXP list with the retained draws of mu (nchain x 2), tau2
 *         (nchain x 2) and alpha (nchain x n), the BayesThresh quantities
 *         sigma, C1 and C2 of each retained sweep, and the posterior mean of
 *         the allocation variables.
 */
SEXP C_BayesMixReg(SEXP y, SEXP padidx, SEXP family, SEXP fs,
                   SEXP waveletfilter, SEXP j0, SEXP hyper, SEXP prior,
                   SEXP init, SEXP mcmc);

#endif /* WAVEBASED_WAV_BAYESMIX_H */
