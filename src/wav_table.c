/**
 * @file wav_table.c
 * @brief Generation of interpolation tables for the scaling and wavelet functions.
 * @details Builds the lookup tables used by PhiVecInterp()/PsiVecInterp() as a
 *          fast alternative to the Daubechies-Lagarias algorithm. The tables
 *          are produced by the exact routines PhiVec()/PsiVec() themselves,
 *          which guarantees that the storage conventions (ordering and
 *          translation offsets) are identical to the exact evaluation path.
 * @author Michel H. Montoril
 */

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "wav_table.h"
#include "wav_utilities.h"
#include "phi_psi_vec.h"

SEXP C_WavTable(SEXP family, SEXP fs, SEXP prec, SEXP waveletfilter, SEXP ngrid){

  int g, G, lkmax, lkmin, m, N, rprec;
  double px, rpsisl;
  double *phi, *prod, *psi, *rphitab, *rpsitab, *rwfilter, *tmp;
  const void *vmax;
  SEXP out, phitab, psitab;

  G = INTEGER(ngrid)[0];
  rprec = INTEGER(prec)[0];

  if(G < 2)
    error("'ngrid' must be an integer greater than or equal to 2.");

  SEXP wutils = PROTECT(WavUtilities(family, fs, waveletfilter));
  SEXP psisl = VECTOR_ELT(wutils, 2);
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  SEXP fsize = VECTOR_ELT(wutils, 5);
  rpsisl = INTEGER(psisl)[0];
  rwfilter = REAL(wfilter);
  N = INTEGER(fsize)[0];

  PROTECT(phitab = allocMatrix(REALSXP, N - 1, G + 1));
  PROTECT(psitab = allocMatrix(REALSXP, N - 1, G + 1));
  rphitab = REAL(phitab);
  rpsitab = REAL(psitab);

  phi = (double *) R_alloc((N - 1), sizeof(double));
  psi = (double *) R_alloc((N - 1), sizeof(double));
  prod = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));
  tmp = (double *) R_alloc((N - 1)*(N - 1), sizeof(double));

  /* Columns g = 0, ..., G - 1: exact Daubechies-Lagarias values at w = g/G.
   * PsiVec() allocates internal buffers with R_alloc on every call, so the
   * allocation stack is rewound at each iteration to keep memory bounded. */
  for(g = 0; g < G; g++){

    px = (double) g / G;

    vmax = vmaxget();

    PhiVec(phi, px, rwfilter, N, rprec, prod, tmp);

    /* Same translation offset used by C_PSImat/PSImat for this point. */
    lkmax = floor(px - rpsisl);
    lkmin = lkmax - N + 2;
    PsiVec(psi, px, rwfilter, N, rprec, lkmin, prod, tmp);

    vmaxset(vmax);

    for(m = 0; m < N - 1; m++){
      rphitab[m + (N - 1)*g] = phi[m];
      rpsitab[m + (N - 1)*g] = psi[m];
    }
  }

  /* Column G (w = 1): the m-th entry equals the (m - 1)-th entry of the w = 0
   * column, because the tabulated arguments shift by one unit; the first entry
   * sits at the right end of the compact support, where phi and psi vanish. */
  rphitab[0 + (N - 1)*G] = 0.0;
  rpsitab[0 + (N - 1)*G] = 0.0;
  for(m = 1; m < N - 1; m++){
    rphitab[m + (N - 1)*G] = rphitab[m - 1];
    rpsitab[m + (N - 1)*G] = rpsitab[m - 1];
  }

  PROTECT(out = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, phitab);
  SET_VECTOR_ELT(out, 1, psitab);

  UNPROTECT(4);
  return out;
}
