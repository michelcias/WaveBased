#ifndef WAVEBASED_WAV_MIXREG_H
#define WAVEBASED_WAV_MIXREG_H

/**
 * @file wav_mixreg.h
 * @brief Wavelet-based estimators for mixture regression (R interface).
 * @details Implements the two raw scaling-coefficient estimators of Montoril,
 *          Pinheiro and Vidakovic (2019), together with the basis synthesis and
 *          the local linear smoother used to regularize them in the time
 *          domain. All routines accumulate directly into the periodized
 *          coefficient vector, so no n x 2^J design matrix is ever formed.
 * @author Michel H. Montoril
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Computes the raw scaling-coefficient estimators of the mixture
 *        regression model.
 *
 * @details Let \f$W_i\f$ be the transformed responses and \f$T_i\f$ the
 *          observation times. The routine evaluates
 *          \f[ \hat c_{Jk} = \frac{1}{s} \sum_i (u_i - l_i)\,\phi_{Jk}(T_i)\,W_i \f]
 *          (estimator (8) of the paper) and
 *          \f[ \tilde c_{Jk} = \frac{1}{s} \sum_i W_i \int_{l_i}^{u_i} \phi_{Jk}(t)\,dt \f]
 *          (estimator (9)), where \f$[l_i, u_i] = [T_{i-s/2}, T_{i+s/2}]\f$.
 *
 *          The integrals of estimator (9) are obtained from the primitive
 *          \f$\Phi(u) = \int_{-\infty}^{u} \phi\f$ through the identity
 *          \f[ \int_a^b \phi_{Jk} = 2^{-J/2}\,[\Phi(2^J b - k) - \Phi(2^J a - k)], \f]
 *          which costs O(n N) instead of the O(n N (b - a)/delta) of the
 *          trapezoidal rule used in the original MATLAB implementation. The
 *          trapezoidal rule is still available for reference and validation.
 *
 * @param[in] w             Real SEXP: transformed responses W of length n.
 * @param[in] x             Real SEXP: observation times T of length n.
 * @param[in] lower         Real SEXP: lower interval end points l of length n.
 * @param[in] upper         Real SEXP: upper interval end points u of length n.
 * @param[in] s             Real SEXP: the constant s of the paper (s = 2 gives
 *                          the interval [T_{i-1}, T_{i+1}]).
 * @param[in] J             Integer SEXP: resolution level.
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets,
 *                          2 = Symmlets, 3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: number of Daubechies-Lagarias iterations.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] what          Integer SEXP: 1 for estimator (8) only, 2 for
 *                          estimator (9) only, 3 for both.
 * @param[in] intmethod     Integer SEXP: 0 for the primitive of phi, 1 for the
 *                          trapezoidal rule.
 * @param[in] intpars       Real SEXP of length 2: the sub-interval length delta
 *                          and the minimum number of points of the trapezoidal
 *                          rule (used only when intmethod == 1).
 * @param[in] evaltab       Real SEXP: the (N - 1) x (G + 1) phi table built by
 *                          C_WavTable() to be used for the pointwise evaluation
 *                          of phi, or R_NilValue for the exact
 *                          Daubechies-Lagarias algorithm.
 * @param[in] ngrid         Integer SEXP: base-two logarithm of the number of
 *                          sub-intervals per unit interval of the dyadic grid
 *                          on which the primitive is built. Used only when
 *                          intmethod == 0.
 * @return SEXP list of length 2 holding the coefficient vectors of estimators
 *         (8) and (9), each of length 2^J, with R_NilValue in the slot of an
 *         estimator that was not requested.
 */
SEXP C_MixCoef(SEXP w, SEXP x, SEXP lower, SEXP upper, SEXP s, SEXP J,
               SEXP family, SEXP fs, SEXP prec, SEXP waveletfilter,
               SEXP what, SEXP intmethod, SEXP intpars, SEXP evaltab,
               SEXP ngrid);

/**
 * @brief Evaluates a periodized scaling-function expansion at arbitrary points.
 *
 * @details Computes \f$f(x_i) = \sum_{k=0}^{2^J - 1} c_k \phi^{(p)}_{Jk}(x_i)\f$
 *          by accumulating the N - 1 non-null terms of each point, which costs
 *          O(M N) and O(1) extra memory instead of the O(M 2^J) of an explicit
 *          basis matrix product.
 *
 * @param[in] coef          Real SEXP: coefficient vector of length 2^J.
 * @param[in] x             Real SEXP: evaluation points of length M.
 * @param[in] J             Integer SEXP: resolution level.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: number of Daubechies-Lagarias iterations.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] phitab        Real SEXP: phi table for interpolated evaluation, or
 *                          R_NilValue for the exact algorithm.
 * @return Real SEXP of length M with the fitted values.
 */
SEXP C_MixEval(SEXP coef, SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec,
               SEXP waveletfilter, SEXP phitab);

/**
 * @brief Local linear kernel smoother with a Gaussian kernel.
 *
 * @details Solves, at every output point, the weighted least squares problem
 *          \f$\min \sum_i (y_i - \alpha - \beta (t - x_i))^2 K((t - x_i)/h)\f$
 *          and returns the intercept, which is the local linear estimate. This
 *          is the regularization in the time domain of Section 3.3.3 of the
 *          paper.
 *
 * @param[in] x    Real SEXP: design points of length n.
 * @param[in] y    Real SEXP: responses of length n.
 * @param[in] xout Real SEXP: points at which the smoother is evaluated, length M.
 * @param[in] bw   Real SEXP: the bandwidth h.
 * @return Real SEXP of length M with the smoothed values.
 */
SEXP C_LocLin(SEXP x, SEXP y, SEXP xout, SEXP bw);

#endif /* WAVEBASED_WAV_MIXREG_H */
