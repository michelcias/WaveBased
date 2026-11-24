/**
 * @file wav_density.c
 * @brief Bona fide correction of the wavelet density estimates (R interface).
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_density.h"
#include "wav_utilities.h"
#include "phi_psi_vec.h"
#include "cdv_edge.h"
#include "utils.h"

SEXP C_WavFineGrid(SEXP coefs, SEXP J, SEXP r, SEXP family, SEXP fs, SEXP prec,
                   SEXP waveletfilter, SEXP periodic, SEXP cdvblocks){

  int k, L, l, lkmin, lkmax, m, mhi, mlo, Nv, p, q, R, rJ, rr, rprec, rper, s, t;
  double acc, sqp, *c, *g, *h, *hrev, *prod, *tmp, *taps, *tapsrev, *evec, *work;
  SEXP out;

  rJ = INTEGER(J)[0];
  rr = INTEGER(r)[0];
  rprec = INTEGER(prec)[0];
  rper = INTEGER(periodic)[0];

  if(rJ < 0 || rr < 0 || rJ + rr > 30)
    error("The dyadic grid must satisfy J >= 0, r >= 0 and J + r <= 30.");
  if(rper != 1 && rper != 2)
    error("The fine grid is available for the periodic and interval bases only.");

  p = 1 << rJ;
  R = 1 << rr;

  if(LENGTH(coefs) != p)
    error("'coefs' must have length 2^J.");
  c = REAL(coefs);

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  h = REAL(VECTOR_ELT(wutils, 4));
  L = INTEGER(VECTOR_ELT(wutils, 5))[0];

  /* taps[(L - 1)*s + t] = phi(s/R + L - 2 - t), the vector of PhiVec at the
   * fractional part s/R, laid out as PHImat reads it */
  prod = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));
  tmp  = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));
  taps = (double *) R_alloc((size_t) (L - 1)*R, sizeof(double));
  for(s = 0; s < R; s++)
    PhiVec(taps + (size_t) (L - 1)*s, (double) s / R, h, L, rprec, prod, tmp);

  PROTECT(out = allocVector(REALSXP, (R_xlen_t) p*R + 1));
  g = REAL(out);
  sqp = sqrt((double) p);

  if(rper == 1){
    /* periodized basis: phi_{Jk}(u) = sqp * phi(2^J u - k), k taken mod p */
    for(l = 0; l <= p*R; l++){
      q = l / R;
      s = l % R;
      lkmin = q - L + 2;
      acc = 0.0;
      for(t = 0; t < L - 1; t++)
        acc += taps[(size_t) (L - 1)*s + t] * c[mod(lkmin + t, p)];
      g[l] = sqp * acc;
    }
  }
  else{
    /* interval basis: Nv left edge functions, the interior translates
     * m = 1, ..., p - L, and Nv right edge functions (as in PHImatCDV) */
    CDVBlocks blk;
    CDVUnpackBlocks(cdvblocks, L, &blk);
    Nv = blk.Nv;

    hrev = (double *) R_alloc(L, sizeof(double));
    for(k = 0; k < L; k++)
      hrev[k] = h[L - 1 - k];
    tapsrev = (double *) R_alloc((size_t) (L - 1)*R, sizeof(double));
    for(s = 0; s < R; s++)
      PhiVec(tapsrev + (size_t) (L - 1)*s, (double) s / R, hrev, L, rprec,
             prod, tmp);

    evec = (double *) R_alloc(Nv, sizeof(double));
    work = (double *) R_alloc(CDV_EDGE_PHI_DYADIC_WORK(Nv, L), sizeof(double));

    for(l = 0; l <= p*R; l++){
      acc = 0.0;

      /* left edge functions, supported where 2^J u < L - 1 */
      if(l < (L - 1)*R){
        CDVEdgePhiVecDyadic(evec, l, R, blk.BL, blk.bL, blk.phi0L, Nv, L,
                            taps, work);
        for(k = 0; k < Nv; k++)
          acc += c[k] * evec[k];
      }

      q = l / R;
      s = l % R;
      lkmax = q;
      lkmin = lkmax - L + 2;
      mlo = lkmin < 1 ? 1 : lkmin;
      mhi = lkmax > p - L ? p - L : lkmax;
      for(m = mlo; m <= mhi; m++)
        acc += c[Nv + m - 1] * taps[(size_t) (L - 1)*s + (m - lkmin)];

      /* right edge functions, at the distance p - 2^J u to the right end */
      if(p*R - l < (L - 1)*R){
        CDVEdgePhiVecDyadic(evec, p*R - l, R, blk.BR, blk.bR, blk.phi0R, Nv,
                            L, tapsrev, work);
        for(k = 0; k < Nv; k++)
          acc += c[p - Nv + k] * evec[k];
      }

      g[l] = sqp * acc;
    }
  }

  UNPROTECT(2);
  return out;
}

SEXP C_BonaFideXi(SEXP f, SEXP w){

  int i, k, m, n, *idx;
  double cand, C, W, xi, *fv, *wv, *rf, *rw;

  n = LENGTH(f);
  if(LENGTH(w) != n)
    error("'f' and 'w' must have the same length.");
  rf = REAL(f);
  rw = REAL(w);

  fv  = (double *) R_alloc(n, sizeof(double));
  wv  = (double *) R_alloc(n, sizeof(double));
  idx = (int *) R_alloc(n, sizeof(int));

  m = 0;
  for(i = 0; i < n; i++){
    if(!R_FINITE(rf[i]) || !R_FINITE(rw[i]))
      error("The estimate is not finite on the domain of the correction.");
    if(rw[i] < 0)
      error("The quadrature weights must be non-negative.");
    if(rw[i] > 0){
      fv[m] = rf[i];
      wv[m] = rw[i];
      idx[m] = m;
      m++;
    }
  }
  if(m == 0)
    error("The quadrature has no positive weight.");

  /* fv ascending, idx carrying the weights along */
  rsort_with_index(fv, idx, m);

  /* Activate the values from the largest down. The k-th largest value stays
   * above the candidate xi exactly while sum_{j < k} W_j (f_(j) - f_(j+1))
   * < 1, W_j being the weight of the j largest values. The sum grows with k,
   * so the first failure ends the search. The largest value always passes,
   * so xi is set. */
  C = 0.0;
  W = 0.0;
  xi = 0.0;
  for(k = m - 1; k >= 0; k--){
    C += wv[idx[k]] * fv[k];
    W += wv[idx[k]];
    cand = (C - 1.0) / W;
    if(fv[k] > cand)
      xi = cand;
    else
      break;
  }

  return ScalarReal(xi);
}
