#ifndef WAVEBASED_WAV_DECOMP1_H
#define WAVEBASED_WAV_DECOMP1_H

/**
 * @file wav_decomp1.h
 * @brief Single-level wavelet decomposition and reconstruction.
 * @author Michel H. Montoril
 * @date 2026
 */

#include "utils.h"

/**
 * @brief Computes the single-level wavelet decomposition of a signal.
 *
 * @param[in]  x      Input signal of length n.
 * @param[in]  n      Length of the input signal (must be even).
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[out] sclc   Scale coefficients (length n/2).
 * @param[out] dtlc   Detail coefficients (length n/2).
 */
void WaveDec1(double *x, int n, double *filter, int N, double *sclc, double *dtlc);

/**
 * @brief Computes the single-level wavelet reconstruction from scale and detail coefficients.
 *
 * @param[in]  sclc   Scale coefficients of length n.
 * @param[in]  dtlc   Detail coefficients of length n.
 * @param[in]  n      Length of sclc and dtlc.
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[out] recvec Reconstructed signal (length 2*n).
 */
void WaveRec1(double *sclc, double *dtlc, int n, double *filter, int N, double *recvec);

/**
 * @brief One analysis step of the CDV (interval) wavelet transform.
 *
 * Applies the boundary-corrected orthogonal filter bank to a vector laid out
 * as [left edge (L/2) | interior translates m = 1..nin-L | right edge (L/2)],
 * producing scale and detail coefficient vectors with the same layout at the
 * coarser level. No periodic wrap-around is used: the edges are handled by
 * the Cohen-Daubechies-Vial edge blocks.
 *
 * @param[in]  v    Input vector of length nin (nin/2 must be >= L + 2 edges).
 * @param[in]  nin  Length of the input vector (a power of 2).
 * @param[in]  h    Low-pass wavelet filter of length L.
 * @param[in]  g    High-pass filter, g[n] = (-1)^n h[L-1-n].
 * @param[in]  L    Filter length.
 * @param[in]  uw   Number of interior columns of the edge wavelet blocks.
 * @param[in]  BL   Left edge scaling block, L/2 x L/2 (column-major).
 * @param[in]  bL   Left interior scaling block, L/2 x (L-1).
 * @param[in]  UL   Left edge wavelet block, L/2 x L/2.
 * @param[in]  uL   Left interior wavelet block, L/2 x uw.
 * @param[in]  BR   Right edge scaling block (from the reversed filter).
 * @param[in]  bR   Right interior scaling block.
 * @param[in]  UR   Right edge wavelet block.
 * @param[in]  uR   Right interior wavelet block.
 * @param[out] sclc Scale coefficients (length nin/2).
 * @param[out] dtlc Detail coefficients (length nin/2).
 */
void WaveDec1CDV(double *v, int nin, double *h, double *g, int L, int uw,
                 double *BL, double *bL, double *UL, double *uL,
                 double *BR, double *bR, double *UR, double *uR,
                 double *sclc, double *dtlc);

/**
 * @brief One synthesis step of the CDV (interval) wavelet transform.
 *
 * Inverse of WaveDec1CDV: reconstructs the finer-level vector from scale
 * and detail coefficients, all in the layout
 * [left edge (L/2) | interior translates | right edge (L/2)].
 *
 * @param[in]  sclc   Scale coefficients (length n).
 * @param[in]  dtlc   Detail coefficients (length n).
 * @param[in]  n      Length of sclc and dtlc.
 * @param[in]  h      Low-pass wavelet filter of length L.
 * @param[in]  g      High-pass filter, g[n] = (-1)^n h[L-1-n].
 * @param[in]  L      Filter length.
 * @param[in]  uw     Number of interior columns of the edge wavelet blocks.
 * @param[in]  BL     Left edge scaling block, L/2 x L/2 (column-major).
 * @param[in]  bL     Left interior scaling block, L/2 x (L-1).
 * @param[in]  UL     Left edge wavelet block, L/2 x L/2.
 * @param[in]  uL     Left interior wavelet block, L/2 x uw.
 * @param[in]  BR     Right edge scaling block (from the reversed filter).
 * @param[in]  bR     Right interior scaling block.
 * @param[in]  UR     Right edge wavelet block.
 * @param[in]  uR     Right interior wavelet block.
 * @param[out] recvec Reconstructed vector (length 2*n).
 */
void WaveRec1CDV(double *sclc, double *dtlc, int n, double *h, double *g,
                 int L, int uw,
                 double *BL, double *bL, double *UL, double *uL,
                 double *BR, double *bR, double *UR, double *uR,
                 double *recvec);

#endif /* WAVEBASED_WAV_DECOMP1_H */
