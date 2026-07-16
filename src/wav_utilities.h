#ifndef WAVEBASED_WAV_UTILITIES_H
#define WAVEBASED_WAV_UTILITIES_H

/**
 * @file wav_utilities.h
 * @brief Wavelet filter and support utilities for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Provides the wavelet filter and support limit points for phi and psi.
 *
 * For computational efficiency, returns a list of 6 elements:
 * [0] left support of phi, [1] right support of phi,
 * [2] left support of psi, [3] right support of psi,
 * [4] filter coefficients, [5] filter size.
 *
 * Supported families: Daublets (1), Symmlets (2), Coiflets (3),
 * and user-supplied filter (4).
 *
 * @param[in] family        Integer SEXP: 1 = Daublets, 2 = Symmlets, 3 = Coiflets, 4 = custom.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @return SEXP list of length 6.
 */
SEXP WavUtilities(SEXP family, SEXP fs, SEXP waveletfilter);

/**
 * @brief Returns the wavelet filter coefficients of a built-in family (R interface).
 *
 * @param[in] family Integer SEXP: wavelet family (1 = Daublets, 2 = Symmlets, 3 = Coiflets).
 * @param[in] fs     Integer SEXP: filter size.
 * @return Real SEXP with the filter coefficients.
 */
SEXP C_GetFilter(SEXP family, SEXP fs);

/**
 * @brief Returns the precomputed CDV boundary blocks of a tabulated filter (R interface).
 *
 * The blocks are returned as the named list expected by CDVUnpackBlocks
 * (BL, bL, UL, uL, phi0L, BR, bR, UR, uR, phi0R, uwidth), or R_NilValue
 * when the requested filter is not tabulated.
 *
 * @param[in] family Integer SEXP: wavelet family (1 = Daublets, 2 = Symmlets).
 * @param[in] fs     Integer SEXP: filter size.
 * @return List SEXP with the boundary blocks, or R_NilValue.
 */
SEXP C_GetCDVBlocks(SEXP family, SEXP fs);

#endif /* WAVEBASED_WAV_UTILITIES_H */
