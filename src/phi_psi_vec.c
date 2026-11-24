/**
 * @file phi_psi_vec.c
 * @brief Computation of scaling function and wavelet vectors at a given point.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include "phi_psi_vec.h"

/* The binary digits of a double in (0, 1) beyond the 1074th are all zero,
 * since it is a multiple of 2^-1074. */
#define PHIVEC_MAXDIG 1100

void PhiVec(double *phi, double px, double *filter, int N, int prec, double *prod, double *tmp){

  double acc, w, *v = prod;
  int d, i, j, jhi, jlo, l, n1 = N - 1;
  unsigned char dig[PHIVEC_MAXDIG];

  w = px - floor(px);

  if(w == 0)
    w = w + 1e-9;
  else if(w == 1)
    w = w - 1e-9;

  // The first prec digits of the dyadic expansion of w {dig[l] --> d[l]}
  for(l = 0; l < prec && l < PHIVEC_MAXDIG; l++){
    w = 2*w;
    d = floor(w);
    dig[l] = (unsigned char) d;
    w = w - d;
  }

  /* --- Calculating phi[0k](px) --- */
  // The values of phi, in reverse order, are the means of the columns of
  // T(d[0]) T(d[1]) ... T(d[prec - 1]), where T(d)[i, j] = sqrt(2)
  // filter[2i - j + d]: the product applied to the vector with all entries
  // 1/(N - 1). The matrices are applied to that vector from the right, one
  // matrix-vector product each, over the band where T(d) is not null.
  for(i = 0; i < n1; i++)
    v[i] = 1.0/n1;

  for(l = prec - 1; l >= 0; l--){
    d = (l < PHIVEC_MAXDIG) ? dig[l] : 0;

    for(i = 0; i < n1; i++){
      jlo = 2*i + d - n1;
      if(jlo < 0)
        jlo = 0;
      jhi = 2*i + d;
      if(jhi > n1 - 1)
        jhi = n1 - 1;

      acc = 0.0;
      for(j = jlo; j <= jhi; j++)
        acc += filter[2*i - j + d]*v[j];
      tmp[i] = M_SQRT2*acc;
    }

    for(i = 0; i < n1; i++)
      v[i] = tmp[i];
  }

  for(i = 0; i < n1; i++)
    phi[N - 2 - i] = v[i];
  /* --- phi[0k](px) calculated --- */
}

void PsiVec(double *psi, double px, double *filter, int N, int prec, int kmin, double *prod, double *tmp){

  double acc, uval, *phi = prod;
  int i, j, n1 = N - 1, off, u;

  // phi[0k](2px), computed in psi and kept in the first N - 1 entries of
  // prod, which PhiVec no longer needs once it returns
  PhiVec(psi, 2*px, filter, N, prec, prod, tmp);
  for(i = 0; i < n1; i++)
    phi[i] = psi[i];

  off = 1 + 2*kmin - (int) floor(2*px);

  /* --- Calculating psi[0k](px) based on {psi <-- Umat %*% phi} --- */
  // Umat[i, j] = (-1)^(u - 1) filter[u], u = 2i + j + off, null unless
  // 0 <= u <= N - 1
  for(i = 0; i < n1; i++){
    acc = 0.0;
    for(j = 0; j < n1; j++){
      u = 2*i + j + off;
      if(u < 0 || u > n1)
        continue;
      uval = (u % 2) ? filter[u] : -filter[u];
      acc += M_SQRT2 * uval * phi[n1 - 1 - j];
    }
    psi[i] = acc;
  }
  /* --- psi[0k](px) calculated --- */
}
