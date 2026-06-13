#ifndef WAVEBASED_PHI_PSI_VEC_H
#define WAVEBASED_PHI_PSI_VEC_H

/**
 * @file phi_psi_vec.h
 * @brief Computation of scaling function and wavelet vectors at a given point.
 * @author Michel Cias
 * @date 2026
 */

#include "utils.h"

/**
 * @brief Computes the vector phi[0k](px) for all k values where phi is non-null.
 *
 * Based on the 'phi' function from 'WAVDE.c' of the wavethresh package.
 *
 * @param[out] phi    Output scaling function values (length N-1).
 * @param[in]  px     Evaluation point.
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[in]  prec   Number of dyadic refinement iterations.
 * @param[in,out] prod Workspace matrix of size (N-1)^2.
 * @param[in,out] tmp  Workspace matrix of size (N-1)^2.
 */
void PhiVec(double *phi, double px, double *filter, int N, int prec, double *prod, double *tmp);

/**
 * @brief Computes the vector psi[0k](px) for all k values where psi is non-null.
 *
 * Uses the result in Theorem 3.5.5 of Vidakovic (2002, p. 90).
 *
 * @param[out] psi    Output wavelet values (length N-1).
 * @param[in]  px     Evaluation point.
 * @param[in]  filter Wavelet filter of length N.
 * @param[in]  N      Filter length.
 * @param[in]  prec   Number of dyadic refinement iterations.
 * @param[in]  kmin   Minimum k index for the support.
 * @param[in,out] prod Workspace matrix of size (N-1)^2.
 * @param[in,out] tmp  Workspace matrix of size (N-1)^2.
 */
void PsiVec(double *psi, double px, double *filter, int N, int prec, int kmin, double *prod, double *tmp);

#endif /* WAVEBASED_PHI_PSI_VEC_H */
