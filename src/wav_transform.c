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

void WaveDecP(const double *x, int n, int j0, double *filter, int N,
              double *out, double *work){

  int i, half, len = n, coarse = 1 << j0;
  double *tmp = work, *sclc = work + n, *dtlc = work + n + n/2;

  for(i = 0; i < n; i++)
    tmp[i] = x[i];

  while(len > coarse){
    WaveDec1(tmp, len, filter, N, sclc, dtlc);
    half = len/2;
    for(i = 0; i < half; i++){
      out[half + i] = dtlc[i];
      tmp[i] = sclc[i];
    }
    len = half;
  }

  for(i = 0; i < len; i++)
    out[i] = tmp[i];
}

void WaveRecP(const double *x, int n, int j0, double *filter, int N,
              double *out, double *work){

  int i, len = 1 << j0;
  double *sclc = work, *dtlc = work + n/2;

  if(len >= n){
    for(i = 0; i < n; i++)
      out[i] = x[i];
    return;
  }

  for(i = 0; i < len; i++){
    sclc[i] = x[i];
    dtlc[i] = x[i + len];
  }
  WaveRec1(sclc, dtlc, len, filter, N, out);

  for(len *= 2; len < n; len *= 2){
    for(i = 0; i < len; i++){
      sclc[i] = out[i];
      dtlc[i] = x[i + len];
    }
    WaveRec1(sclc, dtlc, len, filter, N, out);
  }
}

SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks){

  int i, j0, J, n, N;
  double *dtlc, *rwdx, *rwfilter, *rx, *sclc;
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

    PROTECT(wdx = allocVector(REALSXP, n));

    WaveDecP(rx, n, j0, rwfilter, N, REAL(wdx),
             (double *) R_alloc(WB_WAVE_WORK(n), sizeof(double)));

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

    PROTECT(wrx = allocVector(REALSXP, n));

    WaveRecP(rx, n, j0, rwfilter, N, REAL(wrx),
             (double *) R_alloc(WB_WAVE_WORK(n), sizeof(double)));

    UNPROTECT(2);
    return wrx;
  }
}
