#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP addxSumPart (SEXP s_x_sum_part, SEXP s_bin_mat,
		  SEXP s_i_not_new_not_old,
		  SEXP s_n_xx, SEXP s_n_bin) {
  
  double *x_sum_part = REAL(s_x_sum_part);
  double *bin_mat = REAL(s_bin_mat);
  int *i_not_new_not_old = INTEGER(s_i_not_new_not_old);
  int *n_xx = INTEGER(s_n_xx);
  int *n_bin = INTEGER(s_n_bin);
  int i, j, index;

  for (i = 0; i < (*n_xx); i++) {
    
    x_sum_part[i] = 0.0;
    
    for (j = 0; j < ((*n_bin) - 2); j++) {

      index = (i_not_new_not_old[j] - 1)*(*n_xx) + i;
      
      x_sum_part[i] = x_sum_part[i] + bin_mat[index];
    }
  }

  return R_NilValue;
}
