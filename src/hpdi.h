#ifndef WAVEBASED_HPDI_H
#define WAVEBASED_HPDI_H

/**
 * @file hpdi.h
 * @brief Header file for the highest posterior density interval (HPDI) routine.
 * @author Michel H. Montoril
 *
 * @details Declares the C entry point that computes the shortest contiguous
 *          interval containing a given probability mass of a posterior sample.
 *          Accepts either a numeric vector (a single chain) or a numeric matrix
 *          whose columns are independent chains (e.g. one column per time point
 *          of a latent trajectory).
 *
 * @note This file was brought over from the author's 'pdm' package, where the
 *       same routine supports the summaries of the polynomial dynamic models.
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Locates the shortest window of fixed span in a sorted vector.
 *
 * @param[in]  sorted Ascending-sorted values of length n.
 * @param[in]  n      Number of values.
 * @param[in]  m      Index span between the lower and upper ends of each
 *                    window, m = floor(prob * n), with 1 <= m <= n - 1.
 * @param[out] lower  Lower end of the shortest window.
 * @param[out] upper  Upper end of the shortest window.
 */
void WBShortestWindow(const double *sorted, int n, int m, double *lower,
                      double *upper);

/**
 * @brief Computes the highest posterior density interval (HPDI) of MCMC samples.
 *
 * @details For a vector of n samples, the HPDI at level prob is the shortest
 *          interval [L, U] containing a proportion prob of the draws. It is
 *          obtained by sorting the sample and, over all windows of fixed span
 *          m = floor(prob * n), selecting the one of minimum width
 *          sorted[i + m] - sorted[i]. For a matrix input the same computation
 *          is applied independently to each column.
 *
 *          Unlike an equal-tailed (quantile) interval, the HPDI is the shortest
 *          contiguous interval and may be asymmetric for skewed posteriors.
 *
 * @param[in] data A REALSXP vector (length n) or matrix (n rows by p columns).
 *                 Rows are MCMC samples; columns are independent chains. Values
 *                 are assumed finite (no NA / NaN).
 * @param[in] prob A scalar probability in the open interval (0, 1).
 *
 * @return For a vector input, a named numeric vector of length 2 (lower,
 *         upper). For a matrix input, a p by 2 numeric matrix with column names
 *         lower and upper, one row per input column.
 *
 * @note The input is not modified; each column is copied into a scratch buffer
 *       before sorting.
 */
SEXP C_hpdi(SEXP data, SEXP prob);

#endif /* WAVEBASED_HPDI_H */
