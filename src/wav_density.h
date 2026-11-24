#ifndef WAVEBASED_WAV_DENSITY_H
#define WAVEBASED_WAV_DENSITY_H

/**
 * @file wav_density.h
 * @brief Bona fide correction of the wavelet density estimates (R interface).
 * @author Michel H. Montoril
 * @date 2026
 */

#include <R.h>
#include <Rinternals.h>

/**
 * @brief Evaluates a scaling-function expansion on a fine dyadic grid of [0, 1].
 *
 * @details Returns g(u) = sum_k coefs[k] phi_{Jk}(u) at the points
 *          u = l / 2^(J + r), l = 0, ..., 2^(J + r), with the basis of PHI()
 *          (periodized, or the boundary-corrected basis of Cohen, Daubechies
 *          and Vial). No basis matrix is formed. The point 2^J u has one of
 *          the 2^r fractional parts s / 2^r, and the Daubechies-Lagarias
 *          vector of the scaling function depends on the fractional part only,
 *          so 2^r evaluations of PhiVec give every value needed: g is then a
 *          circular convolution of the coefficients with these vectors. With
 *          the interval basis, the edge functions are evaluated by
 *          CDVEdgePhiVecDyadic from the same vectors (and those of the
 *          reversed filter, for the right edge). The values agree with
 *          PHI(u, J) %*% coefs up to rounding.
 *
 * @param[in] coefs         Real SEXP: coefficients, of length 2^J.
 * @param[in] J             Integer SEXP: resolution level.
 * @param[in] r             Integer SEXP: number of dyadic refinements of the
 *                          grid (2^r points per interval of level J).
 * @param[in] family        Integer SEXP: wavelet family (1 to 4).
 * @param[in] fs            Integer SEXP: filter size.
 * @param[in] prec          Integer SEXP: Daubechies-Lagarias iterations.
 * @param[in] waveletfilter Real SEXP: custom filter (used only when family == 4).
 * @param[in] periodic      Integer SEXP: 1 for periodic, 2 for interval (CDV).
 * @param[in] cdvblocks     List SEXP: CDV boundary blocks (used only when
 *                          periodic == 2).
 * @return Real SEXP of length 2^(J + r) + 1 with the values of g.
 */
SEXP C_WavFineGrid(SEXP coefs, SEXP J, SEXP r, SEXP family, SEXP fs, SEXP prec,
                   SEXP waveletfilter, SEXP periodic, SEXP cdvblocks);

/**
 * @brief Constant of the correction of Glad, Hjort and Ushakov (2003).
 *
 * @details Returns the xi that solves sum_i w[i] * max(0, f[i] - xi) = 1, for
 *          the values f of an estimate at the nodes of a quadrature with
 *          weights w. The left side is continuous, piecewise linear and
 *          decreasing in xi wherever it is positive, so the solution exists,
 *          is unique and is found exactly after sorting f: with the k largest
 *          values active, xi = (sum w f - 1) / (sum w) over them, and k is the
 *          largest for which the k-th largest value exceeds that xi. Then
 *          max(0, f - xi) is the projection of f, in the inner product of the
 *          quadrature, on the set of densities; xi is negative when the
 *          positive part of f has mass below one. Nodes with zero weight are
 *          ignored.
 *
 * @param[in] f Real SEXP: values of the estimate at the nodes.
 * @param[in] w Real SEXP: quadrature weights (non-negative, not all zero).
 * @return Real SEXP of length 1 with xi.
 */
SEXP C_BonaFideXi(SEXP f, SEXP w);

#endif /* WAVEBASED_WAV_DENSITY_H */
