#ifndef WAVEBASED_WAV_FILTERS_CDV_H
#define WAVEBASED_WAV_FILTERS_CDV_H

/**
 * @file wav_filters_cdv.h
 * @brief Registry of the precomputed Cohen-Daubechies-Vial boundary blocks.
 * @author Michel H. Montoril
 * @date 2026
 */

/**
 * @brief One tabulated filter in the CDV boundary block registry.
 *
 * The data pointer holds the left and right edge blocks concatenated as
 * B (Nv x Nv), b (Nv x (L-1)), U (Nv x Nv), u (Nv x uwidth) and phi0 (Nv),
 * all column-major, with Nv = fs/2, first for the left edge and then for
 * the right edge (derived from the reversed filter).
 */
typedef struct {
  int family;          /**< 1 = Daublets, 2 = Symmlets.                    */
  int fs;              /**< Filter size (length of the wavelet filter).    */
  int uwidth;          /**< Interior columns of the edge wavelet blocks.   */
  const double *data;  /**< Coefficients, laid out as described above.     */
} cdv_table_entry;

/** @brief Precomputed boundary blocks (see wav_filters_cdv.c). */
extern const cdv_table_entry cdv_table[];

/** @brief Number of entries in cdv_table. */
extern const int cdv_table_size;

#endif /* WAVEBASED_WAV_FILTERS_CDV_H */
