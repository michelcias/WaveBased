/**
 * @file wav_utilities.c
 * @brief Wavelet filter and support utilities for the R interface.
 * @author Michel H. Montoril
 * @date 2026
 */

#include "wav_utilities.h"
#include "wav_filters_daublets.h"
#include "wav_filters_symmlets.h"
#include "wav_filters_coiflets.h"

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
