/**
 * @file utils.c
 * @brief General-purpose mathematical utility functions.
 * @author Michel H. Montoril
 * @date 2026
 */

#include "utils.h"

void Range(double *xmin, double *xmax, double *x, int n){

   int i, k = 0;
   // Initializing xmin as the first finite observation
   for(i = 0; i < n; i++){
     if(R_FINITE(x[i])){
       xmin[0] = x[i];
       k = i;
       break;
     }
   }
   
   // (Just being careful!) If there is no finite observation
   // in the data x, finish the code
   if(ISNA(xmin[0]))
     return;
   
   // Initializing xmax as the first finite observation
   xmax[0] = xmin[0];
   
   // Selecting the minimum (xmin) and the maximum (xmax)
   // among all the finite observations
   for(i = k; i < n; i++){
     if(xmin[0] > x[i] && R_FINITE(x[i]))
       xmin[0] = x[i];
     
     if(xmax[0] < x[i] && R_FINITE(x[i]))
       xmax[0] = x[i];
   }
   
 }

void Periodize(double *pmat, double *amat, int nrows, int k1, int k2, int p){

  int i, j, k, l;
  
  for(i = 0; i < nrows; i++){
    if(!R_FINITE(amat[i])){
      for(k = 0; k < p; k++)
        pmat[i + nrows*k] = NA_REAL;
      continue;
    }
    
    for(j = k1; j < (k1 + p); j++){
      k = mod(j, p);
      pmat[i + nrows*k] = 0.0;
      for(l = j; l < (k2 + 1); l += p)
        pmat[i + nrows*k] += amat[i + nrows*(l - k1)];
    }
  }
  
}
