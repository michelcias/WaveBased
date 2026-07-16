/**
 * @file cdv_edge.c
 * @brief Pointwise evaluation of the Cohen-Daubechies-Vial boundary functions.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include "cdv_edge.h"
#include "phi_psi_vec.h"

/* Fills wvec[m-1] = phi(z - m), m = 1, ..., L-1, using the
 * Daubechies-Lagarias algorithm. phi, prod and tmp are workspaces of sizes
 * (L-1), (L-1)^2 and (L-1)^2. */
static void interior_phi_window(double *wvec, double z, const double *filter,
                                int L, int prec, double *phi, double *prod,
                                double *tmp) {
  int m, lkmin, lkmax, mlo, mhi;

  for (m = 0; m < L - 1; m++)
    wvec[m] = 0.0;

  lkmax = (int) floor(z);
  lkmin = lkmax - L + 2;
  mlo = lkmin < 1 ? 1 : lkmin;
  mhi = lkmax > L - 1 ? L - 1 : lkmax;
  if (mlo > mhi)
    return;

  PhiVec(phi, z, (double *) filter, L, prec, prod, tmp);
  for (m = mlo; m <= mhi; m++)
    wvec[m - 1] = phi[m - lkmin];
}

void CDVEdgePhiVec(double *out, double y, const double *B, const double *bmat,
                   const double *phi0, int Nv, int L, const double *filter,
                   int prec, double *work) {

  double *P    = work;
  double *Pn   = P    + Nv * Nv;
  double *tv   = Pn   + Nv * Nv;
  double *tv2  = tv   + Nv;
  double *wvec = tv2  + Nv;
  double *phi  = wvec + (L - 1);
  double *prod = phi  + (L - 1);
  double *tmp  = prod + (L - 1) * (L - 1);
  double z, acc;
  int i, j, k, iter;

  if (y <= 0.0) {
    /* limit values at the boundary point itself */
    for (k = 0; k < Nv; k++)
      out[k] = phi0[k];
    return;
  }

  for (k = 0; k < Nv; k++)
    out[k] = 0.0;
  for (i = 0; i < Nv * Nv; i++)
    P[i] = 0.0;
  for (i = 0; i < Nv; i++)
    P[i + Nv * i] = 1.0;

  for (iter = 0; iter < 1200; iter++) {

    if (y >= (double)(L - 1))
      return;                       /* remaining term vanishes exactly */

    z = 2.0 * y;
    interior_phi_window(wvec, z, filter, L, prec, phi, prod, tmp);

    /* out += sqrt(2) * P %*% (bmat %*% wvec) */
    for (k = 0; k < Nv; k++) {
      acc = 0.0;
      for (j = 0; j < L - 1; j++)
        acc += bmat[k + Nv * j] * wvec[j];
      tv[k] = acc;
    }
    for (k = 0; k < Nv; k++) {
      acc = 0.0;
      for (j = 0; j < Nv; j++)
        acc += P[k + Nv * j] * tv[j];
      out[k] += M_SQRT2 * acc;
    }

    /* P = sqrt(2) * P %*% B */
    for (j = 0; j < Nv; j++) {
      for (k = 0; k < Nv; k++) {
        acc = 0.0;
        for (i = 0; i < Nv; i++)
          acc += P[k + Nv * i] * B[i + Nv * j];
        Pn[k + Nv * j] = M_SQRT2 * acc;
      }
    }
    for (i = 0; i < Nv * Nv; i++)
      P[i] = Pn[i];

    y = z;
  }

  /* Pathologically tiny y (below 2^-1200 after doubling): close the
   * recursion with the limit values at 0. */
  for (k = 0; k < Nv; k++) {
    acc = 0.0;
    for (j = 0; j < Nv; j++)
      acc += P[k + Nv * j] * phi0[j];
    out[k] += acc;
  }
}

void CDVEdgePsiVec(double *out, double y, const double *U, const double *umat,
                   int uw, const double *B, const double *bmat,
                   const double *phi0, int Nv, int L, const double *filter,
                   int prec, double *work) {

  double *evec  = work;
  double *wvec  = evec  + Nv;
  double *phi   = wvec  + uw;
  double *prod  = phi   + (L - 1);
  double *tmp   = prod  + (L - 1) * (L - 1);
  double *ework = tmp   + (L - 1) * (L - 1);
  double z, acc;
  int j, k, m, lkmin, lkmax, mlo, mhi;

  z = 2.0 * y;

  /* edge scaling functions at 2y */
  if (z < (double)(L - 1))
    CDVEdgePhiVec(evec, z, B, bmat, phi0, Nv, L, filter, prec, ework);
  else
    for (k = 0; k < Nv; k++)
      evec[k] = 0.0;

  /* interior scaling functions at 2y: wvec[m-1] = phi(2y - m), m = 1..uw */
  for (m = 0; m < uw; m++)
    wvec[m] = 0.0;
  lkmax = (int) floor(z);
  lkmin = lkmax - L + 2;
  mlo = lkmin < 1 ? 1 : lkmin;
  mhi = lkmax > uw ? uw : lkmax;
  if (mlo <= mhi && z > 0.0) {
    PhiVec(phi, z, (double *) filter, L, prec, prod, tmp);
    for (m = mlo; m <= mhi; m++)
      wvec[m - 1] = phi[m - lkmin];
  }

  for (k = 0; k < Nv; k++) {
    acc = 0.0;
    for (j = 0; j < Nv; j++)
      acc += U[k + Nv * j] * evec[j];
    for (j = 0; j < uw; j++)
      acc += umat[k + Nv * j] * wvec[j];
    out[k] = M_SQRT2 * acc;
  }
}
