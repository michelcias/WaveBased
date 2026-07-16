#ifndef WAVEBASED_CDV_EDGE_H
#define WAVEBASED_CDV_EDGE_H

/**
 * @file cdv_edge.h
 * @brief Pointwise evaluation of the Cohen-Daubechies-Vial boundary functions.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <R.h>

/**
 * @brief Number of doubles of workspace required by CDVEdgePhiVec.
 *
 * @param[in] Nv Number of edge functions (filter length / 2).
 * @param[in] L  Filter length.
 */
#define CDV_EDGE_PHI_WORK(Nv, L) \
  (2 * (Nv) * (Nv) + 2 * (Nv) + 2 * ((L) - 1) + 2 * ((L) - 1) * ((L) - 1))

/**
 * @brief Number of doubles of workspace required by CDVEdgePsiVec.
 *
 * @param[in] Nv Number of edge functions (filter length / 2).
 * @param[in] L  Filter length.
 * @param[in] uw Number of interior columns of the edge wavelet block.
 */
#define CDV_EDGE_PSI_WORK(Nv, L, uw) \
  ((Nv) + (uw) + 2 * ((L) - 1) * ((L) - 1) + ((L) - 1) + \
   CDV_EDGE_PHI_WORK(Nv, L))

/**
 * @brief Evaluates the CDV boundary scaling functions at one point.
 *
 * Computes out[k] = phiL_k(y), k = 0, ..., Nv - 1, for the left-edge scaling
 * functions of the interval multiresolution analysis at the reference scale
 * (supports contained in [0, L-1]). The evaluation unrolls the two-scale
 * relation
 *   PhiL(y) = sqrt(2) * B * PhiL(2y) + sqrt(2) * b * (phi(2y - m))_{m=1..L-1}
 * until 2^t y leaves the support, so the recursion terminates exactly; the
 * interior values are obtained by the Daubechies-Lagarias algorithm (PhiVec).
 * The right edge is handled by calling this routine with the blocks derived
 * from the reversed filter and y measured from the right endpoint.
 *
 * @param[out] out    Array of length Nv receiving the values.
 * @param[in]  y      Evaluation point (y >= 0).
 * @param[in]  B      Nv x Nv edge two-scale block (column-major).
 * @param[in]  bmat   Nv x (L-1) interior two-scale block (column-major).
 * @param[in]  phi0   Values of the edge functions at y = 0 (length Nv).
 * @param[in]  Nv     Number of edge functions.
 * @param[in]  L      Filter length.
 * @param[in]  filter Wavelet filter of length L.
 * @param[in]  prec   Daubechies-Lagarias iterations for interior values.
 * @param[in]  work   Workspace of CDV_EDGE_PHI_WORK(Nv, L) doubles.
 */
void CDVEdgePhiVec(double *out, double y, const double *B, const double *bmat,
                   const double *phi0, int Nv, int L, const double *filter,
                   int prec, double *work);

/**
 * @brief Evaluates the CDV boundary wavelets at one point.
 *
 * Computes out[k] = psiL_k(y), k = 0, ..., Nv - 1, for the left-edge wavelets
 * at the reference scale, through the single two-scale step
 *   psiL_k(y) = sqrt(2) [ sum_l U[k,l] phiL_l(2y)
 *                         + sum_{m=1}^{uw} u[k,m] phi(2y - m) ].
 *
 * @param[out] out    Array of length Nv receiving the values.
 * @param[in]  y      Evaluation point (y >= 0).
 * @param[in]  U      Nv x Nv edge wavelet block (column-major).
 * @param[in]  umat   Nv x uw interior wavelet block (column-major).
 * @param[in]  uw     Number of interior columns of the edge wavelet block.
 * @param[in]  B      Nv x Nv edge two-scale block of the scaling functions.
 * @param[in]  bmat   Nv x (L-1) interior two-scale block.
 * @param[in]  phi0   Values of the edge scaling functions at 0 (length Nv).
 * @param[in]  Nv     Number of edge functions.
 * @param[in]  L      Filter length.
 * @param[in]  filter Wavelet filter of length L.
 * @param[in]  prec   Daubechies-Lagarias iterations for interior values.
 * @param[in]  work   Workspace of CDV_EDGE_PSI_WORK(Nv, L, uw) doubles.
 */
void CDVEdgePsiVec(double *out, double y, const double *U, const double *umat,
                   int uw, const double *B, const double *bmat,
                   const double *phi0, int Nv, int L, const double *filter,
                   int prec, double *work);

#endif /* WAVEBASED_CDV_EDGE_H */
