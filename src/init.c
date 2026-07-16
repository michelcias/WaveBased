/**
 * @file init.c
 * @brief R package initialization and routine registration for WaveBased.
 * @details Handles dynamic loading and registration of the compiled C routines
 *          of the WaveBased package. Registers the .Call entry points used by
 *          the wavelet-based estimation functions (scaling-function and wavelet
 *          basis evaluation, wavelet decomposition and reconstruction), ensuring
 *          a proper and safe interface between R and C code. Dynamic symbol
 *          lookup is disabled to enforce encapsulation.
 * @author Michel H. Montoril
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declarations of the C entry points, now implemented in dedicated files. */
#include "wav_basis.h"
#include "wav_transform.h"
#include "wav_utilities.h"
#include "wav_table.h"

//==============================================================================
// METHOD REGISTRATION TABLE
//==============================================================================

/**
 * @brief Static table defining .Call method entries for the R-C interface.
 *
 * @details Maps R-visible names to their corresponding C implementations with
 *          argument counts. R's dynamic loading system uses this table to route
 *          .Call() invocations to the correct C routines.
 *
 *          **Registered routines:**
 *          - C_PHImat:       scaling-function basis matrix (9 args)
 *          - C_PSImat:       mother-wavelet basis matrix (9 args)
 *          - C_WavBasis:     decomposed wavelet basis matrix (10 args)
 *          - C_WaveDec:      wavelet decomposition (7 args)
 *          - C_WaveRec:      wavelet reconstruction (7 args)
 *          - C_GetFilter:    filter coefficients of a built-in family (2 args)
 *          - C_GetCDVBlocks: precomputed CDV boundary blocks (2 args)
 *          - C_WavTable:     phi/psi interpolation tables (5 args)
 *
 * @note Names must match exactly those used in the R wrapper .Call() calls.
 * @note Argument counts are enforced by R's .Call() mechanism at runtime.
 * @note The NULL terminator is required for proper array traversal by R.
 *
 * @see R_registerRoutines
 * @see R_CallMethodDef
 */
static const R_CallMethodDef CallEntries[] = {
  {"_WaveBased_C_PHImat",       (DL_FUNC) &C_PHImat,       9},
  {"_WaveBased_C_PSImat",       (DL_FUNC) &C_PSImat,       9},
  {"_WaveBased_C_WavBasis",     (DL_FUNC) &C_WavBasis,     10},
  {"_WaveBased_C_WaveDec",      (DL_FUNC) &C_WaveDec,      7},
  {"_WaveBased_C_WaveRec",      (DL_FUNC) &C_WaveRec,      7},
  {"_WaveBased_C_GetFilter",    (DL_FUNC) &C_GetFilter,    2},
  {"_WaveBased_C_GetCDVBlocks", (DL_FUNC) &C_GetCDVBlocks, 2},
  {"_WaveBased_C_WavTable",     (DL_FUNC) &C_WavTable,     5},
  {NULL, NULL, 0}
};

//==============================================================================
// PACKAGE INITIALIZATION FUNCTION
//==============================================================================

/**
 * @brief Register compiled routines for the WaveBased package.
 *
 * @details The initialization routine registers every compiled entry point with
 *          R's dynamic loader and disables runtime symbol lookup to enforce
 *          encapsulation, so that R's .Call interface can safely dispatch to the
 *          corresponding C implementations.
 *
 * @param dll Pointer to the DllInfo structure supplied automatically by R during
 *            library loading.
 *
 * @return Nothing. Registration occurs for its side effects on R's loader state.
 *
 * @note The function name must follow the "R_init_<package>" convention so that
 *       R invokes it during library() calls.
 *
 * @see R_registerRoutines
 * @see R_useDynamicSymbols
 */
void R_init_WaveBased(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
