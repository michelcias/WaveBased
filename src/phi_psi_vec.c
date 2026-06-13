/**
 * @file phi_psi_vec.c
 * @brief Computation of scaling function and wavelet vectors at a given point.
 * @author Michel Cias
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include "phi_psi_vec.h"

void PhiVec(double *phi, double px, double *filter, int N, int prec, double *prod, double *tmp){
  
  double tval, w;
  int d, i, ind, j, k, l;
  
  /* --- Calculating phi[0k](px) --- */
  // Initializing prod as an identity matrix
  for(i = 0; i < (N - 1); i++)
    for(j = 0; j < (N - 1); j++)
      prod[i + (N - 1)*j] = (i == j) ? 1.0 : 0.0;
  
  w = px - floor(px);
  
  if(w == 0)
    w = w + 1e-9;
  else if(w == 1)
    w = w - 1e-9;
  
  for(l = 0; l < prec; l++){
    
    // Calculating the l-th coefficient of the diadic expansion of w
    // {d[l] --> d}
    w = 2*w;
    d = floor(w);
    
    for(i = 0; i < (N - 1); i++){
      for(j = 0; j < (N - 1); j++){
        tmp[i + (N - 1)*j] = 0;
        
        for(k = 0; k < (N - 1); k++){
          
          // Calculating {T[k,j] --> tval}
          ind = 2*(k + 1) - (j + 1) - 1 + d;
          if(ind < 0 || ind > (N - 1))
            tval = 0;
          else
            tval = M_SQRT2*filter[ind];
          
          tmp[i + (N - 1)*j] += prod[i + (N - 1)*k]*tval;
        }
      }
    }
    
    for(i = 0; i < (N - 1)*(N - 1); i++)
      prod[i] = tmp[i];
    
    w = w - d; // updating w to the next diadic coefficient
  }
  
  // Calculating phi[i] as the average of the columns from the prod matrix
  for(i = 0; i < (N - 1); i++){
    phi[N - 2 - i] = 0;
    for(j = 0; j < (N - 1); j++){
      phi[N - 2 - i] += prod[i + (N - 1)*j]/(N - 1);
    }
  }
  /* --- phi[0k](px) calculated --- */
}

void PsiVec(double *psi, double px, double *filter, int N, int prec, int kmin, double *prod, double *tmp){
  
  double px2 = 2*px, uval, *phi;
  int i, j, fpx2;
  int *idxu;
  
  fpx2 = floor(px2);
  idxu = (int *) R_alloc((3*(N - 2) + 1), sizeof(int));
  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  
  PhiVec(phi, px2, filter, N, prec, prod, tmp);
  
  for(i = 0; i < (3*(N - 2) + 1); i++)
    idxu[i] = i + 1 + 2 * kmin - fpx2;
  
  /* --- Calculating psi[0k](px) based on {psi <-- Umat %*% phi} --- */
  for(i = 0; i < (N - 1); i++){
    psi[i] = 0;
    for(j = 0; j < (N - 1); j++){
      if(idxu[2*i + j] < 0 || idxu[2*i + j] > (N - 1))
        uval = 0.0;
      else
        uval = pow((-1), idxu[2*i + j] - 1) * filter[idxu[2*i + j]];
      
      psi[i] += M_SQRT2 * uval * phi[(N - 1) - j - 1];
    }
  }
  /* --- psi[0k](px) calculated --- */
}
