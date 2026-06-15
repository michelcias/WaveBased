#ifndef WAVEBASED_WAV_FILTERS_COIFLETS_H
#define WAVEBASED_WAV_FILTERS_COIFLETS_H

/**
 * @file wav_filters_coiflets.h
 * @brief Coiflets wavelet filter coefficients.
 * @author Michel Cias
 * @date 2026
 */

#include <R.h>

/**
 * @brief Fills the Coiflets wavelet filter for a given filter size.
 *
 * All coefficients are taken from PyWavelets with extended precision
 * and scaled by M_SQRT2. Valid filter sizes are 6, 12, 18, 24, ..., 96 and 102.
 *
 * @param[in]  rfs      Filter size (number of coefficients).
 * @param[out] rwfilter Array of length rfs to fill with filter coefficients.
 */
void fill_filter_coiflets(int rfs, double *rwfilter);

#endif /* WAVEBASED_WAV_FILTERS_COIFLETS_H */
