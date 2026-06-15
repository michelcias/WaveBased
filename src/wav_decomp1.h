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

#endif /* WAVEBASED_WAV_DECOMP1_H */
