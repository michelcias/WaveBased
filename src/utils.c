/**
 * @file utils.c
 * @brief General-purpose mathematical utility functions.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <string.h>
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

/**
 * @brief Places the k-th order statistic of an array at position k.
 *
 * @details Hoare's selection algorithm: the array is partitioned around the
 *          element currently sitting at the requested position, and only the
 *          side that contains the order statistic is partitioned again. On
 *          return, x[k] holds the k-th smallest value and every element to its
 *          left is not greater than it.
 *
 * @param[in,out] x Array of length n, permuted by the call.
 * @param[in]     n Length of the array.
 * @param[in]     k Index of the requested order statistic, 0 <= k < n.
 */
static void WBSelect(double *x, int n, int k){

  int l = 0, r = n - 1, i, j;
  double pivot, tmp;

  while(l < r){
    pivot = x[k];
    i = l;
    j = r;
    do{
      while(x[i] < pivot) i++;
      while(pivot < x[j]) j--;
      if(i <= j){
        tmp = x[i]; x[i] = x[j]; x[j] = tmp;
        i++;
        j--;
      }
    } while(i <= j);
    if(j < k) l = i;
    if(k < i) r = j;
  }
}

double WBMedianInPlace(double *x, int n){

  int i, half = n/2;
  double lo;

  if(n <= 0)
    return NA_REAL;

  WBSelect(x, n, half);

  if(n % 2)
    return x[half];

  /* Everything to the left of the selected position is not greater than it,
   * so the preceding order statistic is the maximum of that side. */
  lo = x[0];
  for(i = 1; i < half; i++)
    if(x[i] > lo)
      lo = x[i];

  return 0.5*(lo + x[half]);
}

double WBMad(const double *x, int n, double *work){

  int i;
  double med;

  if(n <= 0)
    return NA_REAL;

  memcpy(work, x, ((size_t) n)*sizeof(double));
  med = WBMedianInPlace(work, n);

  for(i = 0; i < n; i++)
    work[i] = fabs(x[i] - med);

  return WB_MAD_CONSTANT*WBMedianInPlace(work, n);
}
