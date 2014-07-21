#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP binNotFullMean (SEXP s_x_mean, SEXP s_xx, 
		     SEXP s_n_xx, SEXP s_n_bin_use) {
  
  double *x_mean = REAL(s_x_mean);
  double *xx = REAL(s_xx);
  int *n_xx = INTEGER(s_n_xx);
  int *n_bin_use = INTEGER(s_n_bin_use);
  int i;

  for (i = 0; i < (*n_xx); i++) {
    
    x_mean[i] = (((*n_bin_use) - 1)*(x_mean[i]) + xx[i])/(*n_bin_use);
  }
  
  return R_NilValue;
}
