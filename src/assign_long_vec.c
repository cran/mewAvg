#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP assignLongVec(SEXP s_lhs, SEXP s_rhs, SEXP s_n) {
  
  double *lhs = REAL(s_lhs);
  double *rhs = REAL(s_rhs);
  int *n = INTEGER(s_n);
  int i;

  for (i = 0; i < (*n); i++) {
    
    lhs[i] = rhs[i];
  }

  return R_NilValue;
}
