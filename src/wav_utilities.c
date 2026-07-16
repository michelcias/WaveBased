/**
 * @file wav_utilities.c
 * @brief Wavelet filter and support utilities for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include <string.h>
#include "wav_utilities.h"
#include "wav_filters_daublets.h"
#include "wav_filters_symmlets.h"
#include "wav_filters_coiflets.h"
#include "wav_filters_cdv.h"

SEXP WavUtilities(SEXP family, SEXP fs, SEXP waveletfilter){
  
  double *rwfilter;
  int rfam, rfs;
  rfam = INTEGER(family)[0];
  
  SEXP wfilter;
  
  if(rfam == 1){
    rfs = INTEGER(fs)[0];
    PROTECT(wfilter = allocVector(REALSXP, rfs));
    rwfilter = REAL(wfilter);
    fill_filter_daublets(rfs, rwfilter);
  }
  else if(rfam == 2){
    rfs = INTEGER(fs)[0];
    PROTECT(wfilter = allocVector(REALSXP, rfs));
    rwfilter = REAL(wfilter);
    fill_filter_symmlets(rfs, rwfilter);
  }
  else if(rfam == 3){
    rfs = INTEGER(fs)[0];
    PROTECT(wfilter = allocVector(REALSXP, rfs));
    rwfilter = REAL(wfilter);
    fill_filter_coiflets(rfs, rwfilter);
  }
  else if(rfam == 4){
    rfs = LENGTH(waveletfilter);
    wfilter = waveletfilter;
  }
  else
    error("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.");
  
  SEXP results = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(results, 0, ScalarInteger(0));
  SET_VECTOR_ELT(results, 1, ScalarInteger(rfs - 1));
  SET_VECTOR_ELT(results, 2, ScalarInteger(-(rfs/2 - 1)));
  SET_VECTOR_ELT(results, 3, ScalarInteger(rfs/2));
  SET_VECTOR_ELT(results, 4, wfilter);
  SET_VECTOR_ELT(results, 5, ScalarInteger(rfs));
  
  if(rfam == 4)
    UNPROTECT(1);
  else
    UNPROTECT(2);

  return results;
}

SEXP C_GetFilter(SEXP family, SEXP fs){
  SEXP wutils = PROTECT(WavUtilities(family, fs, R_NilValue));
  SEXP wfilter = VECTOR_ELT(wutils, 4);
  UNPROTECT(1);
  return wfilter;
}

/* Copies one block of the precomputed CDV table into a freshly allocated
 * matrix (or vector, when nc == 0) at position pos of the result list,
 * advancing the data cursor. */
static const double *cdv_push_block(SEXP result, int pos, const double *src,
                                    int nr, int nc){
  SEXP blk;
  int len = nc > 0 ? nr*nc : nr;

  if(nc > 0)
    PROTECT(blk = allocMatrix(REALSXP, nr, nc));
  else
    PROTECT(blk = allocVector(REALSXP, nr));

  memcpy(REAL(blk), src, len*sizeof(double));
  SET_VECTOR_ELT(result, pos, blk);
  UNPROTECT(1);
  return src + len;
}

SEXP C_GetCDVBlocks(SEXP family, SEXP fs){

  int i, k, side, L, Nv, uw, rfam, rfs;
  const cdv_table_entry *entry = NULL;
  const double *cursor;
  SEXP result, names;
  static const char *block_names[] = {"BL", "bL", "UL", "uL", "phi0L",
                                      "BR", "bR", "UR", "uR", "phi0R",
                                      "uwidth"};

  rfam = INTEGER(family)[0];
  rfs = INTEGER(fs)[0];

  for(i = 0; i < cdv_table_size; i++){
    if(cdv_table[i].family == rfam && cdv_table[i].fs == rfs){
      entry = &cdv_table[i];
      break;
    }
  }

  if(entry == NULL)
    return R_NilValue;

  L = entry->fs;
  Nv = L/2;
  uw = entry->uwidth;

  PROTECT(result = allocVector(VECSXP, 11));
  cursor = entry->data;
  for(side = 0; side < 2; side++){
    k = 5*side;
    cursor = cdv_push_block(result, k + 0, cursor, Nv, Nv);      /* B    */
    cursor = cdv_push_block(result, k + 1, cursor, Nv, L - 1);   /* b    */
    cursor = cdv_push_block(result, k + 2, cursor, Nv, Nv);      /* U    */
    cursor = cdv_push_block(result, k + 3, cursor, Nv, uw);      /* u    */
    cursor = cdv_push_block(result, k + 4, cursor, Nv, 0);       /* phi0 */
  }
  SET_VECTOR_ELT(result, 10, ScalarInteger(uw));

  PROTECT(names = allocVector(STRSXP, 11));
  for(i = 0; i < 11; i++)
    SET_STRING_ELT(names, i, mkChar(block_names[i]));
  setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(2);
  return result;
}
