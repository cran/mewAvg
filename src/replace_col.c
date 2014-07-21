#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP replaceCol (SEXP s_mat, SEXP s_x, SEXP s_col, SEXP s_nrow) {

  double *mat = REAL(s_mat);
  double *x = REAL(s_x);
  int *col = INTEGER(s_col);
  int *nrow = INTEGER(s_nrow);
  int i;

  for (i = 0; i < (*nrow); i++) {
    
    mat[(*col)*(*nrow) + i] = x[i];
  }

  return R_NilValue;
}
