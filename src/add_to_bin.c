#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP addToBin (SEXP s_bin_mat, SEXP s_sample, 
	       SEXP s_bin, SEXP s_nrow) {
  
  double *bin_mat = REAL(s_bin_mat);
  double *sample = REAL(s_sample);
  int *bin = INTEGER(s_bin);
  int *nrow = INTEGER(s_nrow);
  int i;

  for (i = 0; i < (*nrow); i++) {
    
    bin_mat[(*bin)*(*nrow) + i] = 
      bin_mat[(*bin)*(*nrow) + i] + sample[i];
  }

  return R_NilValue;
}
