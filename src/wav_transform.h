#ifndef WAVEBASED_WAV_TRANSFORM_H
#define WAVEBASED_WAV_TRANSFORM_H

/**
 * @file wav_transform.h
 * @brief Full wavelet decomposition and reconstruction for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Computes the multi-level wavelet decomposition of a signal.
 *
 * @param[in] x             Real SEXP: signal of length n (must be a power of two).
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets, 2 = Symmlets,
 *                          3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] J0            Integer SEXP: coarsest decomposition level.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @return SEXP real vector of length n with wavelet coefficients.
 */
SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter);

/**
 * @brief Computes the multi-level wavelet reconstruction from decomposed coefficients.
 *
 * @param[in] x             Real SEXP: decomposed signal of length n.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] J0            Integer SEXP: coarsest decomposition level.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @return SEXP real vector of length n with the reconstructed signal.
 */
SEXP C_WaveRec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter);

#endif /* WAVEBASED_WAV_TRANSFORM_H */
