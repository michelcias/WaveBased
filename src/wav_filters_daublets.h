#ifndef WAVEBASED_WAV_FILTERS_DAUBLETS_H
#define WAVEBASED_WAV_FILTERS_DAUBLETS_H

/**
 * @file wav_filters_daublets.h
 * @brief Daublets (Daubechies) wavelet filter coefficients.
 * @author Michel Cias
 * @date 2026
 */

#include <R.h>

/**
 * @brief Fills the Daublets wavelet filter for a given filter size.
 *
 * All coefficients are taken from PyWavelets with extended precision.
 * Valid filter sizes are 2, 4, 6, ..., 74 and 76.
 *
 * @param[in]  rfs      Filter size (number of coefficients).
 * @param[out] rwfilter Array of length rfs to fill with filter coefficients.
 */
void fill_filter_daublets(int rfs, double *rwfilter);

#endif /* WAVEBASED_WAV_FILTERS_DAUBLETS_H */
