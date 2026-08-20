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
 * @brief Draws the mean and the precision of one mixture component.
 *
 * @details Both full conditionals are the conjugate ones of expression (3) of
 *          the article: the mean is normal given the previous precision, and
 *          the precision is gamma given the mean just drawn. Components with no
 *          allocated observation fall back to their priors, which the formulas
 *          already produce.
 *
 *          The component parameters of the dynamic mixture do not depend on how
 *          the mixture weights are built, so this routine is shared with the
 *          regime sampler of C_BayesRegime(). It also serves the half-t prior
 *          of the component scales, which is a gamma full conditional of the
 *          same shape: it is enough to call this function with the shape
 *          nu/2 and the rate nu*b, where b is the auxiliary variable of the
 *          scale mixture drawn by WBDrawHalfTAux(). A prior of the component
 *          scales that is not of that form belongs beside these two functions,
 *          and not inside the samplers that call them.
 *
 * @param[in]     y     Observed values of length n.
 * @param[in]     z     Allocation variables of length n.
 * @param[in]     n     Sample size.
 * @param[in]     grp   Group whose parameters are drawn, 0 or 1.
 * @param[in]     b0    Prior mean of mu.
 * @param[in]     B0    Prior variance of mu.
 * @param[in]     v0    Prior shape of tau^2.
 * @param[in]     V0    Prior rate of tau^2.
 * @param[in,out] mu    On entry unused; on exit the drawn mean.
 * @param[in,out] tau2  On entry the current precision, on exit the drawn one.
 */
void WBDrawComponent(const double *y, const int *z, int n, int grp,
                     double b0, double B0, double v0, double V0,
                     double *mu, double *tau2);

/**
 * @brief Draws the auxiliary variable of a half-t prior of a component scale.
 *
 * @details A half-t prior of the standard deviation of a component,
 *          \f$\tau_k^{-1} \sim \mathrm{half}\mathrm{-}t(\nu, A)\f$, is
 *          written as the inverse gamma scale mixture of Wand et al. (2011),
 *          \f[ \tau_k^{-2} \mid a_k \sim \mathrm{IG}(\nu/2, \nu/a_k),
 *              \qquad a_k \sim \mathrm{IG}(1/2, 1/A^2), \f]
 *          which leaves every full conditional of the sampler in closed form.
 *          Only the reciprocal \f$b_k = 1/a_k\f$ is ever needed, and its full
 *          conditional is
 *          \f[ b_k \mid \tau_k^2 \sim
 *              \Gamma\left(\frac{\nu + 1}{2},\ \nu\tau_k^2 +
 *              \frac{1}{A^2}\right), \f]
 *          an exponential distribution when \f$\nu = 1\f$, which is the
 *          half-Cauchy prior of Gelman (2006).
 *
 *          The precision itself is drawn by WBDrawComponent(), called with the
 *          shape \f$\nu/2\f$ and the rate \f$\nu b_k\f$; together with this
 *          function that is the whole of the half-t step.
 *
 * @param[in] tau2  Current precision of the component, that is, the reciprocal
 *                  of its variance.
 * @param[in] df    Degrees of freedom of the half-t prior; one is the
 *                  half-Cauchy.
 * @param[in] scale Scale of the half-t prior.
 * @return The drawn auxiliary variable, floored at the machine epsilon.
 *
 * @see Gelman, A. (2006). Prior distributions for variance parameters in
 *      hierarchical models. Bayesian Analysis, 1(3), 515-534.
 * @see Wand, M. P., Ormerod, J. T., Padoan, S. A. and Fruhwirth, R. (2011).
 *      Mean field variational Bayes for elaborate distributions. Bayesian
 *      Analysis, 6(4), 847-900.
 */
double WBDrawHalfTAux(double tau2, double df, double scale);

/**
 * @brief Draws the allocation variables from their full conditionals.
 *
 * @details This is expression (4) of the article, evaluated on the logarithmic
 *          scale so that the ratio stays defined when a weight reaches zero or
 *          one, or when both normal densities underflow. One uniform variate is
 *          consumed per observation regardless of the path taken, which keeps
 *          the stream of the sampler independent of the data. It is shared with
 *          the regime sampler of C_BayesRegime().
 *
 * @param[in]  y     Observed values of length n.
 * @param[in]  alpha Mixture weights of length n, in [0, 1].
 * @param[in]  n     Sample size.
 * @param[in]  mu    Means of the two components.
 * @param[in]  sd    Standard deviations of the two components.
 * @param[out] z     Allocation variables of length n.
 */
void WBDrawAlloc(const double *y, const double *alpha, int n,
                 const double *mu, const double *sd, int *z);

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
