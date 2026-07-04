#ifndef WAVEBASED_WAV_TABLE_H
#define WAVEBASED_WAV_TABLE_H

/**
 * @file wav_table.h
 * @brief Generation of interpolation tables for the scaling and wavelet functions.
 * @author Michel H. Montoril
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Tabulates the scaling function and the mother wavelet on a fine grid.
 *
 * @details Both PhiVec() and PsiVec() depend on the evaluation point only
 *          through its fractional part w = px - floor(px). This routine calls
 *          them on the regular grid w = g/ngrid, g = 0, ..., ngrid - 1, storing
 *          each returned vector as one column of a (N - 1) x (ngrid + 1)
 *          matrix. The last column (w = 1) is obtained exactly by shifting the
 *          w = 0 column, since phi and psi vanish at the end points of their
 *          compact supports. The resulting tables are independent of the
 *          resolution level J and can be reused by every basis routine.
 *
 * @param[in] family        Integer SEXP: wavelet family (1 = Daublets,
 *                          2 = Symmlets, 3 = Coiflets, 4 = custom).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: number of Daubechies-Lagarias iterations.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] ngrid         Integer SEXP: number of grid intervals in [0, 1].
 * @return SEXP list of length 2: the phi table and the psi table, both
 *         (N - 1) x (ngrid + 1) real matrices whose g-th column holds the
 *         output of PhiVec()/PsiVec() at w = g/ngrid.
 */
SEXP C_WavTable(SEXP family, SEXP fs, SEXP prec, SEXP waveletfilter, SEXP ngrid);

#endif /* WAVEBASED_WAV_TABLE_H */
