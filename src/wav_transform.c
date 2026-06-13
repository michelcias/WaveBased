/**
 * @file wav_transform.c
 * @brief Full wavelet decomposition and reconstruction for the R interface.
 * @author Michel Cias
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_transform.h"
#include "wav_utilities.h"
#include "wav_decomp1.h"

SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter){
  
  int i, j0, j, J, n, N, tmpscl, tmpn;
  double *dtlc, *rwdx, *rwfilter, *rx, *sclc, *tmp;
  SEXP wdx;
  
  rx = REAL(x);
  n = length(x);
  J = log2(n);
  j0 = INTEGER(J0)[0];
  
  if(J != trunc(J))
    error("This sample size must be a power of two.");
  
  if(j0 > J){
    error("2^j0 must be smaller than the sample size of the data.");
  }
  else if(j0 == J){
    warning("2^j0 is equal to the sample size. Therefore, there will be no decomposition.");
    return x;
  }
  else{
    
    SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
    SEXP wfilter = VECTOR_ELT(wutils, 4);
    SEXP fsize = VECTOR_ELT(wutils, 5);
    rwfilter = REAL(wfilter);
    N = INTEGER(fsize)[0];//LENGTH(wfilter);//INTEGER(fs)[0];
    
    sclc = (double *) R_alloc(n/2, sizeof(double));
    dtlc = (double *) R_alloc(n/2, sizeof(double));
    tmp  = (double *) R_alloc(n/2, sizeof(double));
    
    PROTECT(wdx = allocVector(REALSXP, n));
    rwdx = REAL(wdx);
    
    tmpscl = n/2;
    WaveDec1(rx, n, rwfilter, N, sclc, dtlc);
    for(i = 0; i < tmpscl; i++){
      rwdx[i + tmpscl] = dtlc[i];
      tmp[i] = sclc[i];
    }
    
    if(j0 == J - 1){
      for(i = 0; i < tmpscl; i++){
        rwdx[i] = sclc[i];
      }
      
      UNPROTECT(2);
      return wdx;
    }
    
    for(j = J - 2; j > j0; j--){
      tmpn = tmpscl;
      tmpscl /= 2;
      WaveDec1(tmp, tmpn, rwfilter, N, sclc, dtlc);
      for(i = 0; i < tmpscl; i++){
        rwdx[i + tmpscl] = dtlc[i];
        tmp[i] = sclc[i];
      }
    }
    
    tmpn = tmpscl;
    tmpscl /= 2;
    WaveDec1(tmp, tmpn, rwfilter, N, sclc, dtlc);
    for(i = 0; i < tmpscl; i++){
      rwdx[i + tmpscl] = dtlc[i];
      rwdx[i] = sclc[i];
    }
    
    UNPROTECT(2);
    return wdx;
  }
}

SEXP C_WaveRec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter){
  
  int i, j, j0, J, n, N, tmpscl;
  double *dtlc, *sclc, *rwfilter, *rwrx, *rx;
  SEXP wrx;
  
  rx = REAL(x);
  n = length(x);
  J = log2(n);
  j0 = INTEGER(J0)[0];
  
  if(j0 > J)
    error("2^j0 must be smaller than the sample size of the data.");
  else if(j0 == J){
    warning("2^j0 is equal to the sample size. Therefore, there will be no decomposition.");
    return x;
  }
  else{
    
    SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
    SEXP wfilter = VECTOR_ELT(wutils, 4);
    SEXP fsize = VECTOR_ELT(wutils, 5);
    rwfilter = REAL(wfilter);
    N = INTEGER(fsize)[0];//LENGTH(wfilter);//INTEGER(fs)[0];
    
    sclc = (double *) R_alloc(n/2, sizeof(double));
    dtlc = (double *) R_alloc(n/2, sizeof(double));

    PROTECT(wrx = allocVector(REALSXP, n));
    rwrx = REAL(wrx);
    
    tmpscl = pow(2, j0);
    for(i = 0; i < tmpscl; i++){
      sclc[i] = rx[i];
      dtlc[i] = rx[i + tmpscl];
    }
    WaveRec1(sclc, dtlc, tmpscl, rwfilter, N, rwrx);
    
    for(j = j0 + 1; j < J; j++){
      tmpscl *= 2;
      for(i = 0; i < tmpscl; i++){
        sclc[i] = rwrx[i];
        dtlc[i] = rx[i + tmpscl];
      }
      WaveRec1(sclc, dtlc, tmpscl, rwfilter, N, rwrx);
    }
    
    UNPROTECT(2);
    return wrx;
  }
}
