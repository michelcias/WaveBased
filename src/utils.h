#ifndef WAVEBASED_UTILS_H
#define WAVEBASED_UTILS_H

/**
 * @file utils.h
 * @brief General-purpose mathematical utility functions.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <R.h>

/** @brief Modular arithmetic macro that handles negative values correctly. */
#define mod(a, b) ((((a) % (b)) + (b)) % (b))

/**
 * @brief Computes the minimum and maximum of finite values in an array.
 *
 * @param[out] xmin Pointer to store the minimum finite value.
 * @param[out] xmax Pointer to store the maximum finite value.
 * @param[in]  x    Array of doubles to scan.
 * @param[in]  n    Length of the array.
 */
void Range(double *xmin, double *xmax, double *x, int n);

/**
 * @brief Periodizes an auxiliary matrix into a target matrix.
 *
 * @param[out] pmat  Periodized output matrix (nrows x p).
 * @param[in]  amat  Input auxiliary matrix (nrows x (k2-k1+1)).
 * @param[in]  nrows Number of rows.
 * @param[in]  k1    Lower column index bound.
 * @param[in]  k2    Upper column index bound.
 * @param[in]  p     Period length.
 */
void Periodize(double *pmat, double *amat, int nrows, int k1, int k2, int p);

/** @brief Consistency constant of the median absolute deviation, as in mad(). */
#define WB_MAD_CONSTANT 1.4826

/**
 * @brief Median of an array, computed by selection and destroying the input.
 *
 * @details Uses Hoare's selection algorithm, which reorders the array around
 *          the requested order statistic in O(n) expected time instead of the
 *          O(n log n) of a full sort. Even sample sizes average the two central
 *          order statistics, so the result matches median() exactly.
 *
 * @param[in,out] x Array of length n. It is permuted by the call.
 * @param[in]     n Length of the array.
 * @return The median, or NA_REAL when n is not positive.
 */
double WBMedianInPlace(double *x, int n);

/**
 * @brief Median absolute deviation of an array, as computed by mad().
 *
 * @details Returns WB_MAD_CONSTANT times the median of the absolute deviations
 *          from the median. The input is left untouched, all the work being
 *          done in the buffer supplied by the caller, so that the routine can
 *          be called repeatedly (as in an MCMC loop) without allocating.
 *
 * @param[in]  x    Array of length n, not modified.
 * @param[in]  n    Length of the array.
 * @param[out] work Scratch buffer with room for n doubles.
 * @return The median absolute deviation, or NA_REAL when n is not positive.
 */
double WBMad(const double *x, int n, double *work);

#endif /* WAVEBASED_UTILS_H */
