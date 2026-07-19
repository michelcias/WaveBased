#ifndef WAVEBASED_WAV_BASIS_SPARSE_H
#define WAVEBASED_WAV_BASIS_SPARSE_H

/**
 * @file wav_basis_sparse.h
 * @brief Direct sparse (CSC) construction of the periodic decomposed basis.
 * @author Michel H. Montoril
 * @date 2026
 */

/* ------------------------------------------------------------------------
 * Pure computational core (no R API besides R_alloc-free helpers): these
 * two functions only touch caller-provided arrays, so they can be compiled
 * standalone (with a stub of R.h) for bit-level validation and benchmarks,
 * exactly like the harnesses used for the WaveDec1 and C_WavBasis work.
 *
 * Column layout (before the optional drop): [phi_{j0} | psi_{j0} | ... |
 * psi_{J-1}], with 2^J columns in total, matching C_WavBasis. When
 * dropfirst != 0, global column 0 (the first scaling function) is omitted
 * and every other column is shifted left by one, matching the removal of
 * the constant phi_{00} column performed by .wall_design().
 * ------------------------------------------------------------------------ */

/**
 * @brief Counts the structural non-zeros of each column of the sparse basis.
 *
 * @details At the level of period p, each observation activates the N - 1
 *          translations lkmin, ..., lkmax, which reduce to
 *          min(p, N - 1) distinct columns modulo p. The counts follow that
 *          pattern exactly (values are not computed here).
 *
 * @param[in]  x         Observations (already rescaled to [0, 1]), length n.
 * @param[in]  n         Number of observations.
 * @param[in]  j0        Coarsest level.
 * @param[in]  J         Finest multiresolution level (j0 <= J).
 * @param[in]  N         Filter size.
 * @param[in]  phisl     Left support bound of phi (0).
 * @param[in]  psisl     Left support bound of psi (-(N/2 - 1)).
 * @param[in]  dropfirst Non-zero to omit global column 0.
 * @param[out] cnt       Per-column counts, length 2^J - (dropfirst != 0),
 *                       zero-initialized by the caller.
 */
void WavBasisSparseCount(const double *x, int n, int j0, int J, int N,
                         int phisl, int psisl, int dropfirst, int *cnt);

/**
 * @brief Fills the CSC arrays of the sparse periodic decomposed basis.
 *
 * @details Values are computed exactly as in the dense path
 *          (PHImat/PSImat + Periodize): one PhiVec/PsiVec (or table
 *          lookup) evaluation per observation per level, scaled by
 *          sqrt(p), with translations that collide modulo p accumulated
 *          in increasing translation order -- which reproduces the
 *          Periodize() sums bit for bit. Rows arrive in increasing order
 *          within each column, as required by the CSC format.
 *
 * @param[in]  x       Observations, length n.
 * @param[in]  n       Number of observations.
 * @param[in]  j0      Coarsest level.
 * @param[in]  J       Finest multiresolution level.
 * @param[in]  filter  Wavelet filter, length N.
 * @param[in]  N       Filter size.
 * @param[in]  prec    Daubechies-Lagarias iterations (exact path).
 * @param[in]  phisl   Left support bound of phi.
 * @param[in]  psisl   Left support bound of psi.
 * @param[in]  dropfirst Non-zero to omit global column 0.
 * @param[in]  phitab  Interpolation table for phi, or NULL for exact.
 * @param[in]  psitab  Interpolation table for psi, or NULL for exact.
 * @param[in]  G       Number of grid intervals of the tables.
 * @param[in]  colptr  Column pointers (length ncol + 1), from the counts.
 * @param[out] cursor  Workspace, length ncol (overwritten).
 * @param[out] ri      Row indices (0-based), length colptr[ncol].
 * @param[out] vx      Values, length colptr[ncol].
 * @param[out] vals    Workspace, length N - 1.
 * @param[out] prod    Workspace, length (N - 1)^2.
 * @param[out] tmp     Workspace, length (N - 1)^2.
 * @param[out] acc     Workspace, length N - 1.
 */
void WavBasisSparseFill(const double *x, int n, int j0, int J,
                        double *filter, int N, int prec,
                        int phisl, int psisl, int dropfirst,
                        const double *phitab, const double *psitab, int G,
                        const int *colptr, int *cursor, int *ri, double *vx,
                        double *vals, double *prod, double *tmp, double *acc);

#ifndef WAVEBASED_STANDALONE

#include <Rinternals.h>

/**
 * @brief .Call entry point: sparse periodic decomposed basis in CSC format.
 *
 * @details Sparse counterpart of C_WavBasis for the periodic boundary,
 *          returning the design block of one covariate without ever
 *          materializing the dense n x 2^J matrix. The result is a list
 *          with components i (0-based row indices), p (column pointers)
 *          and x (values), ready for Matrix::sparseMatrix(index1 = FALSE).
 *          Explicit zeros may be present wherever the dense path would
 *          also compute a zero inside the structural pattern; the R side
 *          removes them with Matrix::drop0(), matching the behaviour of
 *          Matrix(B, sparse = TRUE) on the dense result.
 *
 * @param[in] x             Real vector of observations in [0, 1] (finite).
 * @param[in] J0            Integer, coarsest level.
 * @param[in] J             Integer, finest level (J0 <= J).
 * @param[in] family        Integer family code (1-4), as in C_WavBasis.
 * @param[in] fs            Integer filter size.
 * @param[in] prec          Integer Daubechies-Lagarias iterations.
 * @param[in] waveletfilter Real vector, custom filter (family == 4).
 * @param[in] wtab          Interpolation tables from wtable(), or NULL.
 * @param[in] dropfirst     Logical, omit the first scaling column.
 * @return List SEXP with components i, p and x.
 */
SEXP C_WavBasisSparse(SEXP x, SEXP J0, SEXP J, SEXP family, SEXP fs,
                      SEXP prec, SEXP waveletfilter, SEXP wtab,
                      SEXP dropfirst);

#endif /* WAVEBASED_STANDALONE */

#endif /* WAVEBASED_WAV_BASIS_SPARSE_H */
