/**
 * @file wav_basis_sparse.c
 * @brief Direct sparse (CSC) construction of the periodic decomposed basis.
 * @details Builds the design block of one covariate -- the same matrix
 *          produced by C_WavBasis with the periodic boundary -- directly in
 *          compressed sparse column format, without the dense n x 2^J
 *          intermediate. The values reproduce the dense path bit for bit:
 *          the same PhiVec/PsiVec (or table lookup) evaluations, the same
 *          sqrt(p) scaling, and sums of colliding translations performed in
 *          the same order used by Periodize().
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include "utils.h"
#include "phi_psi_vec.h"
#include "phi_psi_interp.h"
#include "wav_basis_sparse.h"

/* ------------------------------------------------------------------------
 * Pure computational core (standalone-friendly).
 * ------------------------------------------------------------------------ */

/**
 * @brief Counts the entries contributed by one level block.
 *
 * @details At the block of period p, observation i activates the columns
 *          mod(lkmin, p), ..., mod(lkmax, p): a cyclic run of
 *          min(p, N - 1) distinct columns starting at mod(lkmin, p).
 */
static void CountBlock(const double *x, int n, int p, int sl, int colbase,
                       int N, int dropfirst, int *cnt){

  int c, c0, g, i, len, lkmax, lkmin, u;
  double px;

  len = (N - 1 < p) ? N - 1 : p;

  for(i = 0; i < n; i++){
    px = p*x[i];
    lkmax = floor(px - sl);
    lkmin = lkmax - N + 2;
    c0 = mod(lkmin, p);
    for(u = 0; u < len; u++){
      c = c0 + u;
      if(c >= p)
        c -= p;
      g = colbase + c;
      if(dropfirst){
        if(g == 0)
          continue;
        g--;
      }
      cnt[g]++;
    }
  }
}

void WavBasisSparseCount(const double *x, int n, int j0, int J, int N,
                         int phisl, int psisl, int dropfirst, int *cnt){

  int colbase, j, p;

  p = 1 << j0;
  /* With dropfirst and j0 = 0 the phi block has a single column, the
   * dropped one: skip it entirely. */
  if(!(dropfirst && j0 == 0))
    CountBlock(x, n, p, phisl, 0, N, dropfirst, cnt);

  colbase = p;
  for(j = j0; j < J; j++){
    p = 1 << j;
    CountBlock(x, n, p, psisl, colbase, N, dropfirst, cnt);
    colbase += p;
  }
}

/**
 * @brief Fills the entries contributed by one level block.
 *
 * @details Reproduces PHImat/PSImat + Periodize bit for bit: the value of
 *          column mod(l, p) is the sum of sqp * vals[l - lkmin] over the
 *          active translations l in the same residue class, accumulated in
 *          increasing l order from 0.0 -- the order used by Periodize().
 *          When p >= N - 1 there are no collisions and the entry is the
 *          plain sqp * vals[t] term (adding it to Periodize's initial 0.0
 *          does not change its bits).
 */
static void FillBlock(const double *x, int n, int p, int sl, int colbase,
                      double *filter, int N, int prec, int ispsi,
                      const double *tab, int G, int dropfirst,
                      int *cursor, int *ri, double *vx,
                      double *vals, double *prod, double *tmp, double *acc){

  int c, c0, g, i, k, lkmax, lkmin, t, u;
  double px, sqp;
  const void *vmax;

  sqp = sqrt((double) p);

  for(i = 0; i < n; i++){
    px = p*x[i];
    lkmax = floor(px - sl);
    lkmin = lkmax - N + 2;

    if(ispsi){
      if(tab)
        PsiVecInterp(vals, px, tab, N, G);
      else{
        /* PsiVec allocates internal workspace with R_alloc on every
         * call: release it per observation to keep the peak memory
         * proportional to the filter size, not to n. */
        vmax = vmaxget();
        PsiVec(vals, px, filter, N, prec, lkmin, prod, tmp);
        vmaxset((void *) vmax);
      }
    }
    else{
      if(tab)
        PhiVecInterp(vals, px, tab, N, G);
      else
        PhiVec(vals, px, filter, N, prec, prod, tmp);
    }

    c0 = mod(lkmin, p);

    if(N - 1 <= p){
      /* No collisions: one entry per translation. */
      for(t = 0; t < N - 1; t++){
        c = c0 + t;
        if(c >= p)
          c -= p;
        g = colbase + c;
        if(dropfirst){
          if(g == 0)
            continue;
          g--;
        }
        k = cursor[g]++;
        ri[k] = i;
        vx[k] = sqp * vals[t];
      }
    }
    else{
      /* p < N - 1: translations wrap and collide; accumulate per
       * residue class in increasing translation order. */
      for(u = 0; u < p; u++)
        acc[u] = 0.0;
      for(t = 0; t < N - 1; t++)
        acc[t % p] += sqp * vals[t];
      for(u = 0; u < p; u++){
        c = c0 + u;
        if(c >= p)
          c -= p;
        g = colbase + c;
        if(dropfirst){
          if(g == 0)
            continue;
          g--;
        }
        k = cursor[g]++;
        ri[k] = i;
        vx[k] = acc[u];
      }
    }
  }
}

void WavBasisSparseFill(const double *x, int n, int j0, int J,
                        double *filter, int N, int prec,
                        int phisl, int psisl, int dropfirst,
                        const double *phitab, const double *psitab, int G,
                        const int *colptr, int *cursor, int *ri, double *vx,
                        double *vals, double *prod, double *tmp, double *acc){

  int colbase, g, j, ncw, p;

  ncw = (1 << J) - (dropfirst ? 1 : 0);
  for(g = 0; g < ncw; g++)
    cursor[g] = colptr[g];

  p = 1 << j0;
  if(!(dropfirst && j0 == 0))
    FillBlock(x, n, p, phisl, 0, filter, N, prec, 0, phitab, G,
              dropfirst, cursor, ri, vx, vals, prod, tmp, acc);

  colbase = p;
  for(j = j0; j < J; j++){
    p = 1 << j;
    FillBlock(x, n, p, psisl, colbase, filter, N, prec, 1, psitab, G,
              dropfirst, cursor, ri, vx, vals, prod, tmp, acc);
    colbase += p;
  }
}

/* ------------------------------------------------------------------------
 * R interface.
 * ------------------------------------------------------------------------ */

#ifndef WAVEBASED_STANDALONE

#include <limits.h>
#include "wav_utilities.h"

/** @brief Local copy of the table extractor of wav_basis.c. */
static const double *TableGetSparse(SEXP wtab, int slot, int N, int *G){

  int *dims;
  SEXP tab;

  if(isNull(wtab))
    return NULL;

  tab = VECTOR_ELT(wtab, slot);
  dims = INTEGER(getAttrib(tab, R_DimSymbol));

  if(dims[0] != N - 1)
    error("'wavelet.table' does not match the selected wavelet filter.");

  *G = dims[1] - 1;
  return REAL(tab);
}

SEXP C_WavBasisSparse(SEXP x, SEXP J0, SEXP J, SEXP family, SEXP fs,
                      SEXP prec, SEXP waveletfilter, SEXP wtab,
                      SEXP dropfirst){

  int dropf, g, i, n, N, ncw, nnz, rJ0, rJ, rphisl, rpsisl, rprec;
  int *cnt, *cursor, *rip, *rp;
  double tot;
  double *acc, *prod, *rwfilt, *rx, *tmp, *vals;
  SEXP ans, ivec, names, pvec, xvec;

  rJ0 = INTEGER(J0)[0];
  rJ = INTEGER(J)[0];
  rprec = INTEGER(prec)[0];
  dropf = LOGICAL(dropfirst)[0];

  if(rJ0 > rJ)
    error("The coarsest level can't be greater than the finest level.");

  rx = REAL(x);
  n = LENGTH(x);

  for(i = 0; i < n; i++)
    if(!R_FINITE(rx[i]))
      error("The observations must be finite real values for the sparse basis construction.");

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  rphisl = INTEGER(VECTOR_ELT(wutils, 0))[0];
  rpsisl = INTEGER(VECTOR_ELT(wutils, 2))[0];
  rwfilt = REAL(VECTOR_ELT(wutils, 4));
  N = INTEGER(VECTOR_ELT(wutils, 5))[0];

  int G = 0;
  const double *phitab = TableGetSparse(wtab, 0, N, &G);
  const double *psitab = rJ0 < rJ ? TableGetSparse(wtab, 1, N, &G) : NULL;

  ncw = (1 << rJ) - (dropf ? 1 : 0);

  /* Pass 1: structural counts and column pointers. */
  cnt = (int *) R_alloc(ncw, sizeof(int));
  for(g = 0; g < ncw; g++)
    cnt[g] = 0;
  WavBasisSparseCount(rx, n, rJ0, rJ, N, rphisl, rpsisl, dropf, cnt);

  PROTECT(pvec = allocVector(INTSXP, ncw + 1));
  rp = INTEGER(pvec);
  rp[0] = 0;
  tot = 0.0;
  for(g = 0; g < ncw; g++){
    tot += cnt[g];
    if(tot > INT_MAX)
      error("The sparse basis has too many non-zero entries; use sparse = \"never\".");
    rp[g + 1] = rp[g] + cnt[g];
  }
  nnz = rp[ncw];

  /* Pass 2: values. */
  PROTECT(ivec = allocVector(INTSXP, nnz));
  PROTECT(xvec = allocVector(REALSXP, nnz));
  rip = INTEGER(ivec);

  cursor = (int *) R_alloc(ncw, sizeof(int));
  vals = (double *) R_alloc(N - 1, sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  acc = (double *) R_alloc(N - 1, sizeof(double));

  WavBasisSparseFill(rx, n, rJ0, rJ, rwfilt, N, rprec, rphisl, rpsisl,
                     dropf, phitab, psitab, G, rp, cursor, rip, REAL(xvec),
                     vals, prod, tmp, acc);

  PROTECT(ans = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ans, 0, ivec);
  SET_VECTOR_ELT(ans, 1, pvec);
  SET_VECTOR_ELT(ans, 2, xvec);
  PROTECT(names = allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("i"));
  SET_STRING_ELT(names, 1, mkChar("p"));
  SET_STRING_ELT(names, 2, mkChar("x"));
  setAttrib(ans, R_NamesSymbol, names);

  UNPROTECT(6);
  return ans;
}

#endif /* WAVEBASED_STANDALONE */
