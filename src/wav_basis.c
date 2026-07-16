/**
 * @file wav_basis.c
 * @brief PHI/PSI matrix computation and full wavelet basis for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_basis.h"
#include "wav_utilities.h"
#include "wav_decomp1.h"
#include "phi_psi_vec.h"
#include "phi_psi_interp.h"
#include "cdv_edge.h"
#include "utils.h"

/**
 * @brief Extracts one interpolation table from the list built by wtable().
 *
 * @details Returns NULL when wtab is R_NilValue (exact evaluation requested).
 *          Otherwise validates that the table rows match the filter in use and
 *          returns a pointer to its values, storing the number of grid
 *          intervals in G.
 *
 * @param[in]  wtab List SEXP with the phi and psi tables, or R_NilValue.
 * @param[in]  slot Position of the desired table (0 = phi, 1 = psi).
 * @param[in]  N    Filter size expected by the caller.
 * @param[out] G    Number of grid intervals of the table.
 * @return Pointer to the table values, or NULL for exact evaluation.
 */
static const double *TableGet(SEXP wtab, int slot, int N, int *G){

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

static SEXP PHImatCDV(double *x, int n, int p, double *filter,
                      double *filtrev, int L, int prec, CDVBlocks *blk){

  int i, j, k, m, lkmin, lkmax, mlo, mhi, Nv = blk->Nv;
  double px, y, sqp;
  double *evec, *ework, *phi, *prod, *tmp, *rmat;
  SEXP mat;

  PROTECT(mat = allocMatrix(REALSXP, n, p));
  rmat = REAL(mat);
  sqp = sqrt((double) p);

  evec  = (double *) R_alloc(Nv, sizeof(double));
  ework = (double *) R_alloc(CDV_EDGE_PHI_WORK(Nv, L), sizeof(double));
  phi   = (double *) R_alloc(L - 1, sizeof(double));
  prod  = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));
  tmp   = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));

  for(i = 0; i < n; i++){

    if(!R_FINITE(x[i])){
      warning("At least one observation is not finite and its associated line in the PHI matrix will be set as NA.");
      for(j = 0; j < p; j++)
        rmat[i + n*j] = NA_REAL;
      continue;
    }

    px = p*x[i];
    for(j = 0; j < p; j++)
      rmat[i + n*j] = 0.0;

    /* left edge functions (support [0, L-1] in level coordinates) */
    if(px < (double)(L - 1)){
      CDVEdgePhiVec(evec, px, blk->BL, blk->bL, blk->phi0L, Nv, L, filter,
                    prec, ework);
      for(k = 0; k < Nv; k++)
        rmat[i + n*k] = sqp * evec[k];
    }

    /* interior translates m = 1, ..., p - L */
    lkmax = (int) floor(px);
    lkmin = lkmax - L + 2;
    mlo = lkmin < 1 ? 1 : lkmin;
    mhi = lkmax > p - L ? p - L : lkmax;
    if(mlo <= mhi){
      PhiVec(phi, px, filter, L, prec, prod, tmp);
      for(m = mlo; m <= mhi; m++)
        rmat[i + n*(Nv + m - 1)] = sqp * phi[m - lkmin];
    }

    /* right edge functions (reflected construction) */
    y = (double) p - px;
    if(y < (double)(L - 1)){
      CDVEdgePhiVec(evec, y, blk->BR, blk->bR, blk->phi0R, Nv, L, filtrev,
                    prec, ework);
      for(k = 0; k < Nv; k++)
        rmat[i + n*(p - Nv + k)] = sqp * evec[k];
    }
  }

  UNPROTECT(1);
  return mat;
}

static SEXP PSImatCDV(double *x, int n, int p, double *filter,
                      double *filtrev, double *g, int L, int prec,
                      CDVBlocks *blk){

  int i, j, k, m, r, lkmin, lkmax, mlo, mhi, Nv = blk->Nv, uw = blk->uw;
  double px, y, z, acc, sqp, psupp;
  double *evec, *ework, *phi, *prod, *tmp, *rmat;
  SEXP mat;

  PROTECT(mat = allocMatrix(REALSXP, n, p));
  rmat = REAL(mat);
  sqp = sqrt((double) p);
  psupp = 0.5 * (double)(L - 1 + uw);   /* edge wavelet support bound */

  evec  = (double *) R_alloc(Nv, sizeof(double));
  ework = (double *) R_alloc(CDV_EDGE_PSI_WORK(Nv, L, uw), sizeof(double));
  phi   = (double *) R_alloc(L - 1, sizeof(double));
  prod  = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));
  tmp   = (double *) R_alloc((L - 1)*(L - 1), sizeof(double));

  for(i = 0; i < n; i++){

    if(!R_FINITE(x[i])){
      warning("At least one observation is not finite and its associated line in the PSI matrix will be set as NA.");
      for(j = 0; j < p; j++)
        rmat[i + n*j] = NA_REAL;
      continue;
    }

    px = p*x[i];
    for(j = 0; j < p; j++)
      rmat[i + n*j] = 0.0;

    /* left edge wavelets */
    if(px < psupp){
      CDVEdgePsiVec(evec, px, blk->UL, blk->uL, uw, blk->BL, blk->bL,
                    blk->phi0L, Nv, L, filter, prec, ework);
      for(k = 0; k < Nv; k++)
        rmat[i + n*k] = sqp * evec[k];
    }

    /* interior wavelets: psi[0m](px) = sqrt(2) sum_j g[j] phi(2 px - 2m - j) */
    z = 2.0*px;
    lkmax = (int) floor(z);
    lkmin = lkmax - L + 2;
    mlo = (int) ceil((lkmin - L + 1) / 2.0);
    if(mlo < 1)
      mlo = 1;
    mhi = (int) floor(lkmax / 2.0);
    if(mhi > p - L)
      mhi = p - L;
    if(mlo <= mhi){
      PhiVec(phi, z, filter, L, prec, prod, tmp);
      for(m = mlo; m <= mhi; m++){
        acc = 0.0;
        for(j = 0; j < L; j++){
          r = 2*m + j;
          if(r >= lkmin && r <= lkmax)
            acc += g[j] * phi[r - lkmin];
        }
        rmat[i + n*(Nv + m - 1)] = sqp * M_SQRT2 * acc;
      }
    }

    /* right edge wavelets (reflected construction) */
    y = (double) p - px;
    if(y < psupp){
      CDVEdgePsiVec(evec, y, blk->UR, blk->uR, uw, blk->BR, blk->bR,
                    blk->phi0R, Nv, L, filtrev, prec, ework);
      for(k = 0; k < Nv; k++)
        rmat[i + n*(p - Nv + k)] = sqp * evec[k];
    }
  }

  UNPROTECT(1);
  return mat;
}

SEXP C_PHImat(SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks){

  int i, j, kdiff1, kmax, kmin, lkmin, lkmax, n, N, p, rJ, rper, rprec;
  double rphisl, rphisr, px, x1 = NA_REAL, xn = NA_REAL;
  double *rphimat1, *phi, *rwfilter, *rphimat2, *rx, *prod, *tmp;
  SEXP phimat;

  n = length(x);

  rJ = INTEGER(J)[0];
  rper = INTEGER(periodic)[0];
  rprec = INTEGER(prec)[0];
  rx = REAL(x);
  p = pow(2, rJ);

  /* --- Defining {min(x) --> x1} and {max(x) --> xn} --- */
  Range(&x1, &xn, rx, n);

  if(ISNA(x1))
    error("Check your data. The observations should be real values.");
  /* --- min(x) and max(x) defined --- */

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  SEXP phisl = VECTOR_ELT(wutils, 0);
  SEXP phisr = VECTOR_ELT(wutils, 1);
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  SEXP fsize = VECTOR_ELT(wutils, 5);
  rphisl = INTEGER(phisl)[0];
  rphisr = INTEGER(phisr)[0];
  rwfilter = REAL(wfilter);
  N = INTEGER(fsize)[0];//LENGTH(wfilter);//INTEGER(fs)[0];

  if(rper == 2){
    /* boundary-corrected (CDV) basis on [0, 1] */
    CDVBlocks blk;
    double *hrev;
    CDVUnpackBlocks(cdvblocks, N, &blk);
    hrev = (double *) R_alloc(N, sizeof(double));
    for(i = 0; i < N; i++)
      hrev[i] = rwfilter[N - 1 - i];
    phimat = PHImatCDV(rx, n, p, rwfilter, hrev, N, rprec, &blk);
    UNPROTECT(1);
    return phimat;
  }

  int G = 0;
  const double *phitab = TableGet(wtab, 0, N, &G);

  kmax = floor(p*xn - rphisl);
  kmin = ceil(p*x1 - rphisr + 1e-9);
  kdiff1 = kmax - kmin + 1;

  if(rper){
    rphimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
    PROTECT(phimat = allocMatrix(REALSXP, n, p));
    rphimat2 = REAL(phimat);
  }
  else{
    PROTECT(phimat = allocMatrix(REALSXP, n, kdiff1));
    rphimat1 = REAL(phimat);
  }

  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  // Let's fill the matrix PHI[Jk](x[i])!
  for(i = 0; i < n; i++){

    if(!R_FINITE(rx[i])){
      warning("At least one observation is not finite and its associated line in the PHI matrix will be set as NA.");
      for(j = 0; j < kdiff1; j++)
        rphimat1[i + n*j] = NA_REAL;
      continue;
    }

    px = p*rx[i];
    /* --- Calculating phi[0k](px) --- */
    if(phitab)
      PhiVecInterp(phi, px, phitab, N, G);
    else
      PhiVec(phi, px, rwfilter, N, rprec, prod, tmp);
    /* --- phi[0k](px) calculated --> phi --- */

    lkmax = floor(px - rphisl);
    lkmin = lkmax - N + 2;

    /* --- Putting the phi[Jk](x[i]) in the i-th row of the matrix rphimat1 --- */
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++){
        rphimat1[i + n*j] = sqrt(p) * phi[j];
      }
      for(j = (N - 1); j < kdiff1; j++){
        rphimat1[i + n*j] = 0.0;
      }
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rphimat1[i + n*j] = 0.0;

      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rphimat1[i + n*j] = sqrt(p) * phi[j + kmin - lkmin];

      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rphimat1[i + n*j] = 0.0;
      }
    }
    /* --- phi[Jk](x[i]) put in the i-th row of the matrix rphimat1 --- */
  }

  if(rper)
    Periodize(rphimat2, rphimat1, n, kmin, kmax, p);

  UNPROTECT(2);
  return phimat;
}

SEXP C_PSImat(SEXP x, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks){

  int i, j, kdiff1, kmax, kmin, lkmin, lkmax, n, N, p, rJ, rper, rprec;
  double rpsisl, rpsisr, px, x1 = NA_REAL, xn = NA_REAL;
  double *psi, *rpsimat1, *rpsimat2, *rwfilter, *rx, *prod, *tmp;
  SEXP psimat;

  n = length(x);
  //N = INTEGER(fs)[0];

  rJ = INTEGER(J)[0];
  rper = INTEGER(periodic)[0];
  rprec = INTEGER(prec)[0];
  rx = REAL(x);

  /* --- Defining {min(x) --> x1} and {max(x) --> xn} --- */
  Range(&x1, &xn, rx, n);

  if(ISNA(x1))
    error("Check your data. The observations should be real values.");
  /* --- min(x) and max(x) defined --- */

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  SEXP psisl = VECTOR_ELT(wutils, 2);
  SEXP psisr = VECTOR_ELT(wutils, 3);
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  SEXP fsize = VECTOR_ELT(wutils, 5);
  rpsisl = INTEGER(psisl)[0];
  rpsisr = INTEGER(psisr)[0];
  rwfilter = REAL(wfilter);
  N = INTEGER(fsize)[0];//LENGTH(wfilter);//INTEGER(fs)[0];

  p = pow(2, rJ);

  if(rper == 2){
    /* boundary-corrected (CDV) wavelets on [0, 1] */
    CDVBlocks blk;
    double *hrev, *g;
    CDVUnpackBlocks(cdvblocks, N, &blk);
    hrev = (double *) R_alloc(N, sizeof(double));
    g = (double *) R_alloc(N, sizeof(double));
    for(i = 0; i < N; i++){
      hrev[i] = rwfilter[N - 1 - i];
      g[i] = (i % 2 ? -1.0 : 1.0) * rwfilter[N - 1 - i];
    }
    psimat = PSImatCDV(rx, n, p, rwfilter, hrev, g, N, rprec, &blk);
    UNPROTECT(1);
    return psimat;
  }

  int G = 0;
  const double *psitab = TableGet(wtab, 1, N, &G);

  kmax = floor(p*xn - rpsisl);
  kmin = ceil(p*x1 - rpsisr + 1e-9);
  kdiff1 = kmax - kmin + 1;

  psi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  if(rper){
    rpsimat1 = (double *) R_alloc((n*kdiff1), sizeof(double));
    PROTECT(psimat = allocMatrix(REALSXP, n, p));
    rpsimat2 = REAL(psimat);
  }
  else{
    PROTECT(psimat = allocMatrix(REALSXP, n, kdiff1));
    rpsimat1 = REAL(psimat);
  }

  // Let's fill the matrix PSI[Jk](x[i])!
  for(i = 0; i < n; i++){

    if(!R_FINITE(rx[i])){
      warning("At least one observation is not finite and its associated line in the PSI matrix will be set as NA.");
      for(j = 0; j < kdiff1; j++)
        rpsimat1[i + n*j] = NA_REAL;
      continue;
    }

    px = p * rx[i];
    lkmax = floor(px - rpsisl);
    lkmin = lkmax - N + 2;

    /* --- Calculating phi[0k](px) --- */
    if(psitab)
      PsiVecInterp(psi, px, psitab, N, G);
    else
      PsiVec(psi, px, rwfilter, N, rprec, lkmin, prod, tmp);
    /* --- psi[0k](px) calculated --> psi --- */

    /* --- Putting the psi[Jk](x[i]) in the i-th row of the matrix rpsimat1 --- */
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++)
        rpsimat1[i + n*j] = sqrt(p) * psi[j];

      for(j = (N - 1); j < kdiff1; j++)
        rpsimat1[i + n*j] = 0.0;
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rpsimat1[i + n*j] = 0.0;

      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rpsimat1[i + n*j] = sqrt(p) * psi[j + kmin - lkmin];

      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rpsimat1[i + n*j] = 0.0;
      }
    }
    /* --- psi[Jk](x[i]) put in the i-th row of the matrix rpsimat1 --- */
  }

  if(rper)
    Periodize(rpsimat2, rpsimat1, n, kmin, kmax, p);

  UNPROTECT(2);
  return psimat;
}

SEXP PHImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int phisl, int phisr, int periodic, const double *phitab, int G){

  int i, j, kdiff1, lkmin, lkmax;
  double px, *rphimat1, *phi, *rphimat2, *prod, *tmp;
  SEXP phimat;

  kdiff1 = kmax - kmin + 1;

  if(periodic){
    rphimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
    PROTECT(phimat = allocMatrix(REALSXP, n, p));
    rphimat2 = REAL(phimat);
  }
  else{
    PROTECT(phimat = allocMatrix(REALSXP, n, kdiff1));
    rphimat1 = REAL(phimat);
  }

  phi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  // Let's fill the matrix PHI[Jk](x[i])!
  for(i = 0; i < n; i++){

    if(!R_FINITE(x[i])){
      for(j = 0; j < kdiff1; j++)
        rphimat1[i + n*j] = NA_REAL;
      continue;
    }

    px = p*x[i];
    /* --- Calculating phi[0k](px) --- */
    if(phitab)
      PhiVecInterp(phi, px, phitab, N, G);
    else
      PhiVec(phi, px, filter, N, prec, prod, tmp);
    /* --- phi[0k](px) calculated --> phi --- */

    lkmax = floor(px - phisl);
    lkmin = lkmax - N + 2;

    /* --- Putting the phi[Jk](x[i]) in the i-th row of the matrix rphimat1 --- */
    if(kmin == lkmin){
      for(j = 0; j < (N - 1); j++){
        rphimat1[i + n*j] = sqrt(p) * phi[j];
      }
      for(j = (N - 1); j < kdiff1; j++){
        rphimat1[i + n*j] = 0.0;
      }
    }
    else{
      for(j = 0; j < (lkmin - kmin); j++)
        rphimat1[i + n*j] = 0.0;

      for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
        rphimat1[i + n*j] = sqrt(p) * phi[j + kmin - lkmin];

      if(kmax > lkmax){
        for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
          rphimat1[i + n*j] = 0.0;
      }
    }
    /* --- phi[Jk](x[i]) put in the i-th row of the matrix rphimat1 --- */
  }

  if(periodic)
    Periodize(rphimat2, rphimat1, n, kmin, kmax, p);

  UNPROTECT(1);
  return phimat;
}

SEXP PSImat(double *x, int n, int p, double *filter, int N, int prec, int kmin, int kmax, int psilh, int psirh, int periodic, const double *psitab, int G){

    int i, j, kdiff1, lkmin, lkmax;
    double px, *psi, *rpsimat1, *rpsimat2, *prod, *tmp;
    SEXP psimat;

    kdiff1 = kmax - kmin + 1;

    psi = (double *) R_alloc((N - 1), sizeof(double));
    prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
    tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

    if(periodic){
        rpsimat1 = (double *) R_alloc(n*kdiff1, sizeof(double));
        PROTECT(psimat = allocMatrix(REALSXP, n, p));
        rpsimat2 = REAL(psimat);
    }
    else{
        PROTECT(psimat = allocMatrix(REALSXP, n, kdiff1));
        rpsimat1 = REAL(psimat);
    }

    for(i = 0; i < n; i++){

        if(!R_FINITE(x[i])){
            for(j = 0; j < kdiff1; j++)
                rpsimat1[i + n*j] = NA_REAL;
            continue;
        }

        px = p * x[i];
        lkmax = floor(px - psilh);
        lkmin = lkmax - N + 2;

        /* --- Calculating psi[0k](px) --- */
        if(psitab)
          PsiVecInterp(psi, px, psitab, N, G);
        else
          PsiVec(psi, px, filter, N, prec, lkmin, prod, tmp);
        /* --- psi[0k](px) calculated --> psi --- */

        /* --- Putting the psi[Jk](x[i]) in the i-th row of the matrix rpsimat1 --- */
        if(kmin == lkmin){
            for(j = 0; j < (N - 1); j++)
                rpsimat1[i + n*j] = sqrt(p) * psi[j];

            for(j = (N - 1); j < kdiff1; j++)
                rpsimat1[i + n*j] = 0.0;
        }
        else{
            for(j = 0; j < (lkmin - kmin); j++)
                rpsimat1[i + n*j] = 0.0;

            for(j = (lkmin - kmin); j < (lkmin - kmin + (N - 1)); j++)
                rpsimat1[i + n*j] = sqrt(p) * psi[j + kmin - lkmin];

            if(kmax > lkmax){
                for(j = (lkmin - kmin + (N - 1)); j < kdiff1; j++)
                    rpsimat1[i + n*j] = 0.0;
            }
        }
        /* --- psi[Jk](x[i]) put in the i-th row of the matrix rpsimat1 --- */
    }

    if(periodic)
        Periodize(rpsimat2, rpsimat1, n, kmin, kmax, p);

    UNPROTECT(1);
    return psimat;
}

SEXP C_WavBasis(SEXP x, SEXP J0, SEXP J, SEXP family, SEXP fs, SEXP prec, SEXP periodic, SEXP waveletfilter, SEXP wtab, SEXP cdvblocks){

  double *rx, *rwfilt, x1 = NA_REAL, xn = NA_REAL;
  int i, n, N, rJ0, rJ, kmax, kmin, p, rper, rphisl, rphisr, rpsisl, rpsisr, rprec;

  rJ0 = INTEGER(J0)[0];
  rJ = INTEGER(J)[0];
  rprec = INTEGER(prec)[0];
  rper = INTEGER(periodic)[0];

  if(rJ0 > rJ)
    error("The coarsest level can't be greater than the finest level.");

  rx = REAL(x);
  n = LENGTH(x);

  /* --- Defining {min(x) --> x1} and {max(x) --> xn} --- */
  Range(&x1, &xn, rx, n);

  if(ISNA(x1))
    error("Check your data. The observations should be real valued.");
  /* --- min(x) and max(x) defined --- */

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  SEXP phisl = VECTOR_ELT(wutils, 0);
  SEXP phisr = VECTOR_ELT(wutils, 1);
  SEXP psisl = VECTOR_ELT(wutils, 2);
  SEXP psisr = VECTOR_ELT(wutils, 3);
  SEXP wfilt = VECTOR_ELT(wutils, 4);
  SEXP fsize = VECTOR_ELT(wutils, 5);
  rphisl = INTEGER(phisl)[0];
  rphisr = INTEGER(phisr)[0];
  rpsisl = INTEGER(psisl)[0];
  rpsisr = INTEGER(psisr)[0];
  rwfilt = REAL(wfilt);
  N = INTEGER(fsize)[0];//LENGTH(wfilt);//INTEGER(fs)[0];

  SEXP wmat, pmat;

  if(rJ0 == rJ){

    p = pow(2, rJ);

    if(rper == 2){
      CDVBlocks blk;
      double *hrev;
      CDVUnpackBlocks(cdvblocks, N, &blk);
      hrev = (double *) R_alloc(N, sizeof(double));
      for(i = 0; i < N; i++)
        hrev[i] = rwfilt[N - 1 - i];
      PROTECT(pmat = PHImatCDV(rx, n, p, rwfilt, hrev, N, rprec, &blk));
      UNPROTECT(2);
      return pmat;
    }

    int G = 0;
    const double *phitab = TableGet(wtab, 0, N, &G);

    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-9);

    PROTECT(pmat = PHImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rphisl, rphisr, rper, phitab, G));
    wmat = pmat;

    UNPROTECT(2);
    return wmat;
  }
  else if(rper == 2){

    /* boundary-corrected (CDV) decomposed basis on [0, 1]: build the
     * scaling basis at the finest level and cascade each row through the
     * interval filter bank, exactly as the periodic path does with the
     * periodic filter bank. */
    CDVBlocks blk;
    double *dtlc, *hrev, *g, *rpmat, *rwmat, *sclc, *tmpv;
    int j, k, lev, len;

    CDVUnpackBlocks(cdvblocks, N, &blk);
    hrev = (double *) R_alloc(N, sizeof(double));
    g = (double *) R_alloc(N, sizeof(double));
    for(i = 0; i < N; i++){
      hrev[i] = rwfilt[N - 1 - i];
      g[i] = (i % 2 ? -1.0 : 1.0) * rwfilt[N - 1 - i];
    }

    p = pow(2, rJ);

    PROTECT(wmat = allocMatrix(REALSXP, n, p));
    rwmat = REAL(wmat);

    PROTECT(pmat = PHImatCDV(rx, n, p, rwfilt, hrev, N, rprec, &blk));
    rpmat = REAL(pmat);

    tmpv = (double *) R_alloc(p, sizeof(double));
    sclc = (double *) R_alloc(p/2, sizeof(double));
    dtlc = (double *) R_alloc(p/2, sizeof(double));

    for(i = 0; i < n; i++){
      for(j = 0; j < p; j++)
        tmpv[j] = rpmat[i + n*j];

      len = p;
      for(lev = rJ; lev > rJ0; lev--){
        WaveDec1CDV(tmpv, len, rwfilt, g, N, blk.uw,
                    blk.BL, blk.bL, blk.UL, blk.uL,
                    blk.BR, blk.bR, blk.UR, blk.uR,
                    sclc, dtlc);
        for(k = 0; k < len/2; k++){
          rwmat[i + n*(len/2 + k)] = dtlc[k];
          tmpv[k] = sclc[k];
        }
        len /= 2;
      }

      for(k = 0; k < len; k++)
        rwmat[i + n*k] = tmpv[k];
    }

    UNPROTECT(3);
    return wmat;
  }
  else if(rper){

    double *dtlc, *rpmat, *rwmat, *sclc, *tmp;
    int i, j, k, tmpn, tmpscl;

    p = pow(2, rJ);

    int G = 0;
    const double *phitab = TableGet(wtab, 0, N, &G);

    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-9);

    PROTECT(wmat = allocMatrix(REALSXP, n, p));
    rwmat = REAL(wmat);

    PROTECT(pmat = PHImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rphisl, rphisr, rper, phitab, G));
    rpmat = REAL(pmat);

    sclc = (double *) R_alloc(p/2, sizeof(double));
    dtlc = (double *) R_alloc(p/2, sizeof(double));
    tmp  = (double *) R_alloc(p, sizeof(double));

    if(rJ0 == rJ - 1){
      for(i = 0; i < n; i++){
        for(j = 0; j < p; j++)
          tmp[j] = rpmat[i + n*j];

        tmpscl = p/2;
        WaveDec1(tmp, p, rwfilt, N, sclc, dtlc);

        for(j = 0; j < tmpscl; j++){
          rwmat[i + n*(j + tmpscl)] = dtlc[j];
          rwmat[i + n*j] = sclc[j];
        }
      }

      UNPROTECT(3);
      return wmat;
    }
    else{
      for(i = 0; i < n; i++){
        for(j = 0; j < p; j++)
          tmp[j] = rpmat[i + n*j];

        tmpscl = p/2;
        WaveDec1(tmp, p, rwfilt, N, sclc, dtlc);

        for(j = 0; j < tmpscl; j++){
          rwmat[i + n*(j + tmpscl)] = dtlc[j];
          tmp[j] = sclc[j];
        }

        for(j = rJ - 2; j > rJ0; j--){
          tmpn = tmpscl;
          tmpscl /= 2;
          WaveDec1(tmp, tmpn, rwfilt, N, sclc, dtlc);
          for(k = 0; k < tmpscl; k++){
            rwmat[i + n*(k + tmpscl)] = dtlc[k];
            tmp[k] = sclc[k];
          }
        }

        tmpn = tmpscl;
        tmpscl /= 2;
        WaveDec1(tmp, tmpn, rwfilt, N, sclc, dtlc);
        for(k = 0; k < tmpscl; k++){
          rwmat[i + n*(k + tmpscl)] = dtlc[k];
          rwmat[i + n*k] = sclc[k];
        }
      }

      UNPROTECT(3);
      return wmat;
    }
  }
  else{

    double *rpmat, *rwmat;
    int j, k, nc0, nc1, ncw;

    int G = 0;
    const double *phitab = TableGet(wtab, 0, N, &G);
    const double *psitab = TableGet(wtab, 1, N, &G);

    p = pow(2, rJ0);
    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-9);

    ncw = (kmax - kmin + 1);
    for(k = rJ0; k < rJ; k++){
      p = pow(2, k);
      kmax = floor(p*xn - rpsisl);
      kmin = ceil(p*x1 - rpsisr + 1e-9);
      ncw += (kmax - kmin + 1);
    }

    PROTECT(wmat = allocMatrix(REALSXP, n, ncw));
    rwmat = REAL(wmat);

    p = pow(2, rJ0);
    kmax = floor(p*xn - rphisl);
    kmin = ceil(p*x1 - rphisr + 1e-9);

    /* --- Calculating the PHI matrix to put in rpmat --- */
    PROTECT(pmat = PHImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rphisl, rphisr, rper, phitab, G));
    rpmat = REAL(pmat);

    nc0 = 0;
    nc1 = rper ? p : (kmax - kmin + 1);

    for(i = 0; i < n; i++)
      for(j = nc0; j < nc1; j++)
        rwmat[i + n*j] = rpmat[i + n*(j - nc0)];

    UNPROTECT(1);
    /* --- PHI matrix put --- */

    /* --- Calculating the PSI matrices to put in rpmat --- */
    for(k = rJ0; k < rJ; k++){
      p = pow(2, k);
      kmax = floor(p*xn - rpsisl);
      kmin = ceil(p*x1 - rpsisr + 1e-9);
      PROTECT(pmat = PSImat(rx, n, p, rwfilt, N, rprec, kmin, kmax, rpsisl, rpsisr, rper, psitab, G));
      rpmat = REAL(pmat);

      nc0 = nc1;
      nc1 += (rper ? p : (kmax - kmin + 1));

      for(i = 0; i < n; i++)
        for(j = nc0; j < nc1; j++)
          rwmat[i + n*j] = rpmat[i + n*(j - nc0)];

      UNPROTECT(1);
      /* --- PSI matrices put --- */
    }

  }

  UNPROTECT(2);
  return wmat;

}
