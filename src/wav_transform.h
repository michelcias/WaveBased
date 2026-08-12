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
 * @brief Size, in doubles, of the scratch buffer required by WaveDecP and
 *        WaveRecP for a signal of length n.
 */
#define WB_WAVE_WORK(n) (2*(n))

/**
 * @brief Multi-level periodic wavelet decomposition (internal).
 *
 * @details Cascades the single-level pyramid algorithm of WaveDec1() down to
 *          the resolution level j0, writing the coefficients in the layout
 *          \f$(c_{j_0,\cdot}, d_{j_0,\cdot}, \ldots, d_{J-1,\cdot})\f$ used
 *          throughout the package. It performs no allocation and no R memory
 *          management, so it can be called repeatedly from a tight loop (as in
 *          the Gibbs sampler of C_BayesMixReg()) with the same buffers.
 *
 * @param[in]  x      Signal of length n, not modified.
 * @param[in]  n      Length of the signal, a power of two.
 * @param[in]  j0     Coarsest resolution level, 0 <= j0 <= log2(n).
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[out] out    Coefficient vector of length n. May not alias x.
 * @param[out] work   Scratch buffer with room for WB_WAVE_WORK(n) doubles.
 */
void WaveDecP(const double *x, int n, int j0, double *filter, int N,
              double *out, double *work);

/**
 * @brief Multi-level periodic wavelet reconstruction (internal).
 *
 * @details Inverse of WaveDecP(), cascading WaveRec1() from the level j0 up to
 *          the finest one. As WaveDecP(), it neither allocates nor touches the
 *          R heap.
 *
 * @param[in]  x      Coefficient vector of length n, not modified.
 * @param[in]  n      Length of the vector, a power of two.
 * @param[in]  j0     Coarsest resolution level of the decomposition.
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[out] out    Reconstructed signal of length n. May not alias x.
 * @param[out] work   Scratch buffer with room for WB_WAVE_WORK(n) doubles.
 */
void WaveRecP(const double *x, int n, int j0, double *filter, int N,
              double *out, double *work);

/**
 * @brief Computes the multi-level wavelet decomposition of a signal.
 *
 * @param[in] x             Real SEXP: signal of length n (must be a power of two).
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets, 2 = Symmlets,
 *                          3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] J0            Integer SEXP: coarsest decomposition level.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] boundary      Integer SEXP: 1 for periodic, 2 for interval (CDV).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when boundary == 2).
 * @return SEXP real vector of length n with wavelet coefficients.
 */
SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks);

/**
 * @brief Computes the multi-level wavelet reconstruction from decomposed coefficients.
 *
 * @param[in] x             Real SEXP: decomposed signal of length n.
 * @param[in] family        Integer SEXP: wavelet family.
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] J0            Integer SEXP: coarsest decomposition level.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] boundary      Integer SEXP: 1 for periodic, 2 for interval (CDV).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when boundary == 2).
 * @return SEXP real vector of length n with the reconstructed signal.
 */
SEXP C_WaveRec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks);

#endif /* WAVEBASED_WAV_TRANSFORM_H */
