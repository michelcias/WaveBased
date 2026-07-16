#ifndef WAVEBASED_WAV_BASIS_H
#define WAVEBASED_WAV_BASIS_H

/**
 * @file wav_basis.h
 * @brief PHI/PSI matrix computation and full wavelet basis for the R interface.
 * @author Michel H. Montoril
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Computes the matrix of PHI[Jk](x[i]) for all non-null k values (R interface).
 *
 * @param[in] x             Real SEXP: evaluation points.
 * @param[in] J             Integer SEXP: resolution level.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: dyadic refinement precision.
 * @param[in] periodic      Integer SEXP: 0 for raw, 1 for periodic, 2 for interval (CDV).
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] wtab          List SEXP with the phi and psi interpolation tables
 *                          built by C_WavTable(), or R_NilValue for the exact
 *                          Daubechies-Lagarias evaluation. Ignored when
 *                          periodic == 2 (the boundary evaluation is always
 *                          exact).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when periodic == 2).
 * @return SEXP matrix of PHI values.
 */
SEXP C_PHImat(SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks);

/**
 * @brief Computes the matrix of PSI[Jk](x[i]) for all non-null k values (R interface).
 *
 * @param[in] x             Real SEXP: evaluation points.
 * @param[in] J             Integer SEXP: resolution level.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: dyadic refinement precision.
 * @param[in] periodic      Integer SEXP: 0 for raw, 1 for periodic, 2 for interval (CDV).
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] wtab          List SEXP with the phi and psi interpolation tables
 *                          built by C_WavTable(), or R_NilValue for the exact
 *                          Daubechies-Lagarias evaluation. Ignored when
 *                          periodic == 2 (the boundary evaluation is always
 *                          exact).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when periodic == 2).
 * @return SEXP matrix of PSI values.
 */
SEXP C_PSImat(SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks);

/**
 * @brief Computes the matrix of PHI[Jk](x[i]) for all non-null k values (internal).
 *
 * Used internally by C_WavBasis.
 *
 * @param[in]  x        Array of evaluation points of length n.
 * @param[in]  n        Number of evaluation points.
 * @param[in]  p        2^J (number of columns for periodic case).
 * @param[in]  filter   Wavelet filter of length N.
 * @param[in]  N        Filter length.
 * @param[in]  prec     Dyadic refinement precision.
 * @param[in]  kmin     Minimum k index.
 * @param[in]  kmax     Maximum k index.
 * @param[in]  phisl    Left support of phi.
 * @param[in]  phisr    Right support of phi.
 * @param[in]  periodic 1 for periodic, 0 otherwise.
 * @param[in]  phitab   Interpolation table for phi ((N - 1) x (G + 1),
 *                      column-major) or NULL for the exact evaluation.
 * @param[in]  G        Number of grid intervals of the table (ignored when
 *                      phitab is NULL).
 * @return SEXP matrix of PHI values.
 */
SEXP PHImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int phisl, int phisr, int periodic, const double *phitab, int G);

/**
 * @brief Computes the matrix of PSI[Jk](x[i]) for all non-null k values (internal).
 *
 * Used internally by C_WavBasis.
 *
 * @param[in]  x        Array of evaluation points of length n.
 * @param[in]  n        Number of evaluation points.
 * @param[in]  p        2^J.
 * @param[in]  filter   Wavelet filter of length N.
 * @param[in]  N        Filter length.
 * @param[in]  prec     Dyadic refinement precision.
 * @param[in]  kmin     Minimum k index.
 * @param[in]  kmax     Maximum k index.
 * @param[in]  psilh    Left support of psi.
 * @param[in]  psirh    Right support of psi.
 * @param[in]  periodic 1 for periodic, 0 otherwise.
 * @param[in]  psitab   Interpolation table for psi ((N - 1) x (G + 1),
 *                      column-major) or NULL for the exact evaluation.
 * @param[in]  G        Number of grid intervals of the table (ignored when
 *                      psitab is NULL).
 * @return SEXP matrix of PSI values.
 */
SEXP PSImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int psilh, int psirh, int periodic, const double *psitab, int G);

/**
 * @brief Computes the full wavelet basis matrix combining PHI and PSI blocks.
 *
 * Returns a matrix with PHI[J0k](x[i]) and PSI[jk](x[i]) for J0 <= j <= J.
 *
 * @param[in] x             Real SEXP: evaluation points.
 * @param[in] J0            Integer SEXP: coarsest level.
 * @param[in] J             Integer SEXP: finest level.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: dyadic refinement precision.
 * @param[in] periodic      Integer SEXP: 0 for raw, 1 for periodic, 2 for interval (CDV).
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] wtab          List SEXP with the phi and psi interpolation tables
 *                          built by C_WavTable(), or R_NilValue for the exact
 *                          Daubechies-Lagarias evaluation. Ignored when
 *                          periodic == 2 (the boundary evaluation is always
 *                          exact).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when periodic == 2).
 * @return SEXP matrix with PHI and PSI columns.
 */
SEXP C_WavBasis(SEXP x, SEXP J0, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks);

#endif /* WAVEBASED_WAV_BASIS_H */
