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

#endif /* WAVEBASED_UTILS_H */
