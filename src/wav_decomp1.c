/**
 * @file wav_decomp1.c
 * @brief Single-level wavelet decomposition and reconstruction.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include "wav_decomp1.h"

void WaveDec1(double *x, int n, double *filter, int N, double *sclc, double *dtlc){
   
  int i, j;
  
  /* --- Calculating the scale and detail coefficients --- */
  for(i = 0; i < n/2; i++){
    sclc[i] = 0.0;
    dtlc[i] = 0.0;
    for(j = 0; j < N; j++){
      sclc[i] += x[mod(j + 2*i, n)] * filter[j];
      dtlc[i] += x[mod(1 - j + 2*i, n)] * filter[j] * pow(-1, j + 1);
    }
  }
  /* --- Scale and detail coefficients calculated --- */
}

void WaveRec1(double *sclc, double *dtlc, int n, double *filter, int N, double *recvec){
    
  int i, j, k;
  
  /* --- Calculating the reconstructed vector --- */
  for(i = 0; i < 2*n; i++){
    recvec[i] = 0.0;
    for(j = 0; j < N/2; j ++){
      k = i%2;
      recvec[i] += sclc[mod((i - 2*j - k)/2, n)]*filter[2*j + k] + dtlc[mod((i + 2*j - k)/2, n)]*filter[2*j + 1 - k]*pow(-1, 2*j - k);
    }
  }
  /* --- Reconstructed vector calculated --- */
}
