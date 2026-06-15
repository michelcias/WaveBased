#ifndef WAVEBASED_WAV_FILTERS_SYMMLETS_H
#define WAVEBASED_WAV_FILTERS_SYMMLETS_H

/**
 * @file wav_filters_symmlets.h
 * @brief Symmlets (least asymmetric Daubechies) wavelet filter coefficients.
 * @author Michel Cias
 * @date 2026
 */

#include <R.h>

/**
 * @brief Fills the Symmlets wavelet filter for a given filter size.
 *
 * The Symmlets family represents the least asymmetric Daubechies family.
 * The k-tap least asymmetric filters (k = 8, 10, 12, 14, 20) are flipped
 * to maintain the pattern of Daubechies (1992).
 * Valid filter sizes are 4, 6, 8, ..., 38 and 40.
 *
 * @param[in]  rfs      Filter size (number of coefficients).
 * @param[out] rwfilter Array of length rfs to fill with filter coefficients.
 */
void fill_filter_symmlets(int rfs, double *rwfilter);

#endif /* WAVEBASED_WAV_FILTERS_SYMMLETS_H */
