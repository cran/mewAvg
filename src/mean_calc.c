#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

SEXP meanCalc (SEXP s_x_mean, SEXP s_bin_mat, 
	       SEXP s_x_sum_part, SEXP s_fract,
	       SEXP s_n_xx, SEXP s_i_old,
	       SEXP s_i_new, SEXP s_n_sample,
	       SEXP s_n_part) {
  
  double *x_mean = REAL(s_x_mean);
  double *bin_mat = REAL(s_bin_mat);
  double *x_sum_part = REAL(s_x_sum_part);
  double *fract = REAL(s_fract);
  int *n_xx = INTEGER(s_n_xx);
  int *i_old = INTEGER(s_i_old);
  int *i_new = INTEGER(s_i_new);
  int *n_sample = INTEGER(s_n_sample);
  int *n_part = INTEGER(s_n_part);
  int i, i_old_index, i_new_index;

  for (i = 0; i < (*n_xx); i++) {
    
    i_old_index = (*i_old)*(*n_xx) + i;
    i_new_index = (*i_new)*(*n_xx) + i;
    x_mean[i] = (bin_mat[i_old_index]*(*fract) + 
		 x_sum_part[i] + 
		 bin_mat[i_new_index])/
      (n_sample[(*i_old)]*(*fract) + 
       (*n_part) + 
       n_sample[(*i_new)]);
  }

  return R_NilValue;
}
