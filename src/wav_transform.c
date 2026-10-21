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

/* Interval (CDV) counterparts of WaveDecP() and WaveRecP(), with the same
 * layout of the coefficients and the same scratch size, WB_WAVE_WORK(n) */
static void WaveDecCDV(const double *x, int n, int j0, double *filter,
                       double *g, int N, CDVBlocks *blk, double *out,
                       double *work){

  int i, len = n;
  double *tmpv = work, *sclc = work + n, *dtlc = work + n + n/2;

  for(i = 0; i < n; i++)
    tmpv[i] = x[i];

  while(len > (1 << j0)){
    WaveDec1CDV(tmpv, len, filter, g, N, blk->uw,
                blk->BL, blk->bL, blk->UL, blk->uL,
                blk->BR, blk->bR, blk->UR, blk->uR,
                sclc, dtlc);
    for(i = 0; i < len/2; i++){
      out[i + len/2] = dtlc[i];
      tmpv[i] = sclc[i];
    }
    len /= 2;
  }

  for(i = 0; i < len; i++)
    out[i] = tmpv[i];
}

static void WaveRecCDV(const double *x, int n, int j0, double *filter,
                       double *g, int N, CDVBlocks *blk, double *out,
                       double *work){

  int i, len = 1 << j0;
  double *sclc = work, *dtlc = work + n/2;

  for(i = 0; i < len; i++){
    sclc[i] = x[i];
    dtlc[i] = x[i + len];
  }
  WaveRec1CDV(sclc, dtlc, len, filter, g, N, blk->uw,
              blk->BL, blk->bL, blk->UL, blk->uL,
              blk->BR, blk->bR, blk->UR, blk->uR, out);

  for(len *= 2; len < n; len *= 2){
    for(i = 0; i < len; i++){
      sclc[i] = out[i];
      dtlc[i] = x[i + len];
    }
    WaveRec1CDV(sclc, dtlc, len, filter, g, N, blk->uw,
                blk->BL, blk->bL, blk->UL, blk->uL,
                blk->BR, blk->bR, blk->UR, blk->uR, out);
  }
}

/* Decomposition (inverse = 0) or reconstruction (inverse = 1) of a vector, or
 * of each row of a matrix. The filter and the CDV blocks are set up once, and
 * every row goes through the same buffers. */
static SEXP WaveRows(SEXP x, SEXP family, SEXP fs, SEXP J0,
                     SEXP waveletfilter, SEXP boundary, SEXP cdvblocks,
                     int inverse){

  int i, j0, J, n, N, nr;
  R_xlen_t k, r;
  double *g = NULL, *orow, *out, *rwfilter, *rx, *work, *xrow;
  CDVBlocks blk;
  SEXP res;

  rx = REAL(x);
  nr = isMatrix(x) ? nrows(x) : 1;
  n = isMatrix(x) ? ncols(x) : length(x);

  if(n < 2 || (n & (n - 1)) != 0)
    error("This sample size must be a power of two, and at least two.");
  J = (int) round(log2((double) n));
  j0 = INTEGER(J0)[0];

  if(j0 > J)
    error("2^j0 must be smaller than the sample size of the data.");
  if(j0 == J){
    warning("2^j0 is equal to the sample size. Therefore, there will be no decomposition.");
    return x;
  }

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rwfilter = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  if(INTEGER(boundary)[0] == 2){
    /* boundary-corrected (CDV) transform: cascade the interval filter bank
     * instead of the periodic one */
    CDVUnpackBlocks(cdvblocks, N, &blk);
    g = (double *) R_alloc(N, sizeof(double));
    for(i = 0; i < N; i++)
      g[i] = (i % 2 ? -1.0 : 1.0) * rwfilter[N - 1 - i];
  }

  PROTECT(res = isMatrix(x) ? allocMatrix(REALSXP, nr, n)
                            : allocVector(REALSXP, n));
  out = REAL(res);
  xrow = (double *) R_alloc(n, sizeof(double));
  orow = (double *) R_alloc(n, sizeof(double));
  work = (double *) R_alloc(WB_WAVE_WORK(n), sizeof(double));

  for(r = 0; r < nr; r++){
    for(k = 0; k < n; k++)
      xrow[k] = rx[r + k*nr];

    if(g == NULL && !inverse)
      WaveDecP(xrow, n, j0, rwfilter, N, orow, work);
    else if(g == NULL)
      WaveRecP(xrow, n, j0, rwfilter, N, orow, work);
    else if(!inverse)
      WaveDecCDV(xrow, n, j0, rwfilter, g, N, &blk, orow, work);
    else
      WaveRecCDV(xrow, n, j0, rwfilter, g, N, &blk, orow, work);

    for(k = 0; k < n; k++)
      out[r + k*nr] = orow[k];

    if(r % 1024 == 1023)
      R_CheckUserInterrupt();
  }

  UNPROTECT(2);
  return res;
}

SEXP C_WaveDec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks){
  return WaveRows(x, family, fs, J0, waveletfilter, boundary, cdvblocks, 0);
}

SEXP C_WaveRec(SEXP x, SEXP family, SEXP fs, SEXP J0, SEXP waveletfilter, SEXP boundary, SEXP cdvblocks){
  return WaveRows(x, family, fs, J0, waveletfilter, boundary, cdvblocks, 1);
}
