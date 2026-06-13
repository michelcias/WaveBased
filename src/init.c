/**
 * @file init.c
 * @brief R entry point registration for the WaveBased package.
 * @author Michel Cias
 * @date 2026
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "wav_transform.h"
#include "wav_basis.h"
#include "wav_utilities.h"

static const R_CallMethodDef CallEntries[] = {
    {"C_WaveDec",    (DL_FUNC) &C_WaveDec,    5},
    {"C_WaveRec",    (DL_FUNC) &C_WaveRec,    5},
    {"C_PHImat",     (DL_FUNC) &C_PHImat,     7},
    {"C_PSImat",     (DL_FUNC) &C_PSImat,     7},
    {"C_WavBasis",   (DL_FUNC) &C_WavBasis,   8},
    {"WavUtilities", (DL_FUNC) &WavUtilities, 3},
    {NULL, NULL, 0}
};

void R_init_WaveBased(DllInfo *dll){
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
