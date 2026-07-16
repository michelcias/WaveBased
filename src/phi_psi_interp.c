/**
 * @file phi_psi_interp.c
 * @brief Fast evaluation of scaling function and wavelet vectors by table lookup.
 * @details Replaces the per-point Daubechies-Lagarias iteration (prec matrix
 *          products of dimension (N - 1) x (N - 1)) with a linear interpolation
 *          between two precomputed columns, at a cost of O(N) operations per
 *          evaluation point. The tables are built by C_WavTable().
 * @author Michel H. Montoril
 */

#include <math.h>
#include "phi_psi_interp.h"

void PhiVecInterp(double *phi, double px, const double *tab, int N, int G){

  int g, m;
  double t, u, w;
  const double *c0, *c1;

  w = px - floor(px);

  u = w * G;
  g = (int) floor(u);
  if(g < 0)
    g = 0;
  if(g > G - 1)
    g = G - 1;
  t = u - g;

  c0 = tab + (N - 1)*g;
  c1 = c0 + (N - 1);

  for(m = 0; m < N - 1; m++)
    phi[m] = (1.0 - t)*c0[m] + t*c1[m];
}

void PsiVecInterp(double *psi, double px, const double *tab, int N, int G){

  int g, m;
  double t, u, w;
  const double *c0, *c1;

  w = px - floor(px);

  u = w * G;
  g = (int) floor(u);
  if(g < 0)
    g = 0;
  if(g > G - 1)
    g = G - 1;
  t = u - g;

  c0 = tab + (N - 1)*g;
  c1 = c0 + (N - 1);

  for(m = 0; m < N - 1; m++)
    psi[m] = (1.0 - t)*c0[m] + t*c1[m];
}
