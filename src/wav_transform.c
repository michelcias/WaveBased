/**
 * @file wav_transform.c
 * @brief Full wavelet decomposition and reconstruction for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_transform.h"
#include "wav_utilities.h"
#include "wav_decomp1.h"
#include "cdv_edge.h"

SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks){
  
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

    if(INTEGER(boundary)[0] == 2){
      /* boundary-corrected (CDV) analysis: cascade the interval filter
       * bank instead of the periodic one */
      CDVBlocks blk;
      double *g, *tmpv;
      int len;

      CDVUnpackBlocks(cdvblocks, N, &blk);
      g = (double *) R_alloc(N, sizeof(double));
      for(i = 0; i < N; i++)
        g[i] = (i % 2 ? -1.0 : 1.0) * rwfilter[N - 1 - i];

      sclc = (double *) R_alloc(n/2, sizeof(double));
      dtlc = (double *) R_alloc(n/2, sizeof(double));
      tmpv = (double *) R_alloc(n, sizeof(double));

      PROTECT(wdx = allocVector(REALSXP, n));
      rwdx = REAL(wdx);

      for(i = 0; i < n; i++)
        tmpv[i] = rx[i];

      len = n;
      while(len > (1 << j0)){
        WaveDec1CDV(tmpv, len, rwfilter, g, N, blk.uw,
                    blk.BL, blk.bL, blk.UL, blk.uL,
                    blk.BR, blk.bR, blk.UR, blk.uR,
                    sclc, dtlc);
        for(i = 0; i < len/2; i++){
          rwdx[i + len/2] = dtlc[i];
          tmpv[i] = sclc[i];
        }
        len /= 2;
      }
      for(i = 0; i < len; i++)
        rwdx[i] = tmpv[i];

      UNPROTECT(2);
      return wdx;
    }

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

SEXP C_WaveRec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks){
  
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

    if(INTEGER(boundary)[0] == 2){
      /* boundary-corrected (CDV) synthesis */
      CDVBlocks blk;
      double *g;

      CDVUnpackBlocks(cdvblocks, N, &blk);
      g = (double *) R_alloc(N, sizeof(double));
      for(i = 0; i < N; i++)
        g[i] = (i % 2 ? -1.0 : 1.0) * rwfilter[N - 1 - i];

      sclc = (double *) R_alloc(n/2, sizeof(double));
      dtlc = (double *) R_alloc(n/2, sizeof(double));

      PROTECT(wrx = allocVector(REALSXP, n));
      rwrx = REAL(wrx);

      tmpscl = pow(2, j0);
      for(i = 0; i < tmpscl; i++){
        sclc[i] = rx[i];
        dtlc[i] = rx[i + tmpscl];
      }
      WaveRec1CDV(sclc, dtlc, tmpscl, rwfilter, g, N, blk.uw,
                  blk.BL, blk.bL, blk.UL, blk.uL,
                  blk.BR, blk.bR, blk.UR, blk.uR, rwrx);

      for(j = j0 + 1; j < J; j++){
        tmpscl *= 2;
        for(i = 0; i < tmpscl; i++){
          sclc[i] = rwrx[i];
          dtlc[i] = rx[i + tmpscl];
        }
        WaveRec1CDV(sclc, dtlc, tmpscl, rwfilter, g, N, blk.uw,
                    blk.BL, blk.bL, blk.UL, blk.uL,
                    blk.BR, blk.bR, blk.UR, blk.uR, rwrx);
      }

      UNPROTECT(2);
      return wrx;
    }

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
