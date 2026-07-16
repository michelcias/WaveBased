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

void WaveDec1CDV(double *v, int nin, double *h, double *g, int L, int uw,
                 double *BL, double *bL, double *UL, double *uL,
                 double *BR, double *bR, double *UR, double *uR,
                 double *sclc, double *dtlc){

  int Nv = L / 2;
  int nint_in = nin - L;          /* interior translates m = 1..nint_in    */
  int nout = nin / 2;
  int nint_out = nout - L;        /* interior translates kept at the output */
  int j, k, m, m1;
  double acc;

  /* Input layout: v = [edgeL (Nv) | interior m = 1..nint_in | edgeR (Nv)].
   * Interior index m lives at v[Nv + m - 1]; the mirrored interior index
   * (m-th translate counted from the right end) lives at v[Nv + nint_in - m]. */

  /* --- left edge --- */
  for(k = 0; k < Nv; k++){
    acc = 0.0;
    for(j = 0; j < Nv; j++)
      acc += BL[k + Nv*j] * v[j];
    for(m = 1; m <= L - 1; m++)
      acc += bL[k + Nv*(m - 1)] * v[Nv + m - 1];
    sclc[k] = acc;

    acc = 0.0;
    for(j = 0; j < Nv; j++)
      acc += UL[k + Nv*j] * v[j];
    for(m = 1; m <= uw; m++)
      acc += uL[k + Nv*(m - 1)] * v[Nv + m - 1];
    dtlc[k] = acc;
  }

  /* --- interior --- */
  for(m1 = 1; m1 <= nint_out; m1++){
    acc = 0.0;
    for(j = 0; j < L; j++)
      acc += h[j] * v[Nv + 2*m1 + j - 1];
    sclc[Nv + m1 - 1] = acc;

    acc = 0.0;
    for(j = 0; j < L; j++)
      acc += g[j] * v[Nv + 2*m1 + j - 1];
    dtlc[Nv + m1 - 1] = acc;
  }

  /* --- right edge --- */
  for(k = 0; k < Nv; k++){
    acc = 0.0;
    for(j = 0; j < Nv; j++)
      acc += BR[k + Nv*j] * v[Nv + nint_in + j];
    for(m = 1; m <= L - 1; m++)
      acc += bR[k + Nv*(m - 1)] * v[Nv + nint_in - m];
    sclc[nout - Nv + k] = acc;

    acc = 0.0;
    for(j = 0; j < Nv; j++)
      acc += UR[k + Nv*j] * v[Nv + nint_in + j];
    for(m = 1; m <= uw; m++)
      acc += uR[k + Nv*(m - 1)] * v[Nv + nint_in - m];
    dtlc[nout - Nv + k] = acc;
  }
}

void WaveRec1CDV(double *sclc, double *dtlc, int n, double *h, double *g,
                 int L, int uw,
                 double *BL, double *bL, double *UL, double *uL,
                 double *BR, double *bR, double *UR, double *uR,
                 double *recvec){

  int Nv = L / 2;
  int nin = 2*n;                  /* length of the reconstructed vector    */
  int nint_in = nin - L;          /* interior translates at the fine level */
  int nint_out = n - L;
  int j, k, m, m1;
  double dk, sk;

  /* The synthesis step is the transpose of WaveDec1CDV: each coarse
   * coefficient scatters its edge/interior filter row into the fine
   * vector, laid out as [edgeL (Nv) | interior m = 1..nint_in | edgeR]. */

  for(j = 0; j < nin; j++)
    recvec[j] = 0.0;

  /* --- left edge --- */
  for(k = 0; k < Nv; k++){
    sk = sclc[k];
    dk = dtlc[k];
    for(j = 0; j < Nv; j++)
      recvec[j] += BL[k + Nv*j]*sk + UL[k + Nv*j]*dk;
    for(m = 1; m <= L - 1; m++)
      recvec[Nv + m - 1] += bL[k + Nv*(m - 1)]*sk;
    for(m = 1; m <= uw; m++)
      recvec[Nv + m - 1] += uL[k + Nv*(m - 1)]*dk;
  }

  /* --- interior --- */
  for(m1 = 1; m1 <= nint_out; m1++){
    sk = sclc[Nv + m1 - 1];
    dk = dtlc[Nv + m1 - 1];
    for(j = 0; j < L; j++)
      recvec[Nv + 2*m1 + j - 1] += h[j]*sk + g[j]*dk;
  }

  /* --- right edge --- */
  for(k = 0; k < Nv; k++){
    sk = sclc[n - Nv + k];
    dk = dtlc[n - Nv + k];
    for(j = 0; j < Nv; j++)
      recvec[Nv + nint_in + j] += BR[k + Nv*j]*sk + UR[k + Nv*j]*dk;
    for(m = 1; m <= L - 1; m++)
      recvec[Nv + nint_in - m] += bR[k + Nv*(m - 1)]*sk;
    for(m = 1; m <= uw; m++)
      recvec[Nv + nint_in - m] += uR[k + Nv*(m - 1)]*dk;
  }
}
