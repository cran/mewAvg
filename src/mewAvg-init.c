#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Visibility.h>
#include <R_ext/Rdynload.h>

SEXP addToBin (SEXP s_bin_mat, SEXP s_sample, SEXP s_bin, SEXP s_nrow);
SEXP addxSumPart (SEXP s_x_sum_part, SEXP s_bin_mat, SEXP s_i_not_new_not_old, 
                  SEXP s_n_xx, SEXP s_n_bin);
SEXP assignLongVec(SEXP s_lhs, SEXP s_rhs, SEXP s_n);
SEXP binNotFullMean (SEXP s_x_mean, SEXP s_xx, 
                     SEXP s_n_xx, SEXP s_n_bin_use);
SEXP meanCalc(SEXP s_x_mean, SEXP s_bin_mat,
              SEXP s_x_sum_part, SEXP s_fract,
              SEXP s_n_xx, SEXP s_i_old,
              SEXP s_i_new, SEXP s_n_sample,
              SEXP s_n_part);
SEXP replaceCol (SEXP s_mat, SEXP s_x, SEXP s_col, SEXP s_nrow);

static const R_CallMethodDef callMethods[] = {
    {"addToBin", (DL_FUNC)&addToBin, 4},
    {"addxSumPart", (DL_FUNC)&addxSumPart, 5},
    {"assignLongVec", (DL_FUNC)&assignLongVec, 3},
    {"binNotFullMean", (DL_FUNC)&binNotFullMean, 4},
    {"meanCalc", (DL_FUNC)&meanCalc, 9},
    {"replaceCol", (DL_FUNC)&replaceCol, 4},
    {NULL, NULL, 0}};

void attribute_visible R_init_mewAvg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
