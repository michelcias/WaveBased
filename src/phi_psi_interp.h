#ifndef WAVEBASED_PHI_PSI_INTERP_H
#define WAVEBASED_PHI_PSI_INTERP_H

/**
 * @file phi_psi_interp.h
 * @brief Fast evaluation of scaling function and wavelet vectors by table lookup.
 * @author Michel H. Montoril
 */

/**
 * @brief Approximates the output of PhiVec() by linear interpolation of a
 *        precomputed table.
 *
 * @details The exact routine PhiVec() depends on the evaluation point only
 *          through its fractional part w = px - floor(px). Given a table whose
 *          g-th column holds the exact PhiVec() output at w = g/G (built by
 *          C_WavTable()), this routine locates the grid cell containing w and
 *          interpolates linearly between the two neighbouring columns. It is a
 *          drop-in replacement for PhiVec() in the basis-matrix routines.
 *
 * @param[out] phi Output vector of length N - 1, same convention as PhiVec().
 * @param[in]  px  Evaluation point (only its fractional part is used).
 * @param[in]  tab Table of exact values, (N - 1) x (G + 1), column-major.
 * @param[in]  N   Filter size.
 * @param[in]  G   Number of grid intervals in [0, 1].
 */
void PhiVecInterp(double *phi, double px, const double *tab, int N, int G);

/**
 * @brief Approximates the output of PsiVec() by linear interpolation of a
 *        precomputed table.
 *
 * @details Same lookup-and-interpolate scheme as PhiVecInterp(), applied to
 *          the mother-wavelet table built by C_WavTable(). Drop-in replacement
 *          for PsiVec() in the basis-matrix routines.
 *
 * @param[out] psi Output vector of length N - 1, same convention as PsiVec().
 * @param[in]  px  Evaluation point (only its fractional part is used).
 * @param[in]  tab Table of exact values, (N - 1) x (G + 1), column-major.
 * @param[in]  N   Filter size.
 * @param[in]  G   Number of grid intervals in [0, 1].
 */
void PsiVecInterp(double *psi, double px, const double *tab, int N, int G);

#endif /* WAVEBASED_PHI_PSI_INTERP_H */
