// Adam Pintar 7/11/2013

// This is the C function that will be called
// by the R wrapper mewAvg

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"mewAvgForRpkg.h"
#include<R.h>
#include<Rdefines.h>

SEXP RmewAvg (SEXP f,          // the R function that generates the 
	                       // random sequence
	      SEXP s_nBin,     // the number of bins to pass to 
	                       // mewAccum and mewMean
	      SEXP s_nXX,      // the length of the vector returned 
	                       // by f
	      SEXP s_ff,       // the fraction of the sample to keep 
	                       // at each iteration
	      SEXP s_nSave,    // the nuber of iterations at which 
	                       // to save the mean
	      SEXP s_nIter,    // the nuber of iterations to perform
	      SEXP s_iToSave,  // the iterations to save coded as 
	                       // 0 (don't) and 1 (do)
	                       // has length nIter
	      SEXP s_meanVals, // the R matrix in which to store the 
	                       // mean values that should be saved
	                       // has dimension nSave by nXX
	      SEXP rho         // the environment in which to 
	                       // evaluate f
	      ) {
  int k, l, totalSaved=0;
  mewTyp av;

  
  PROTECT(s_nBin = AS_INTEGER(s_nBin));
  PROTECT(s_nXX = AS_INTEGER(s_nXX));
  PROTECT(s_nSave = AS_INTEGER(s_nSave));
  PROTECT(s_nIter = AS_INTEGER(s_nIter));
  PROTECT(s_iToSave = AS_INTEGER(s_iToSave));
  PROTECT(s_meanVals = AS_NUMERIC(s_meanVals));
  PROTECT(s_ff = AS_NUMERIC(s_ff));
  int *nBin = INTEGER_POINTER(s_nBin);
  int *nXX = INTEGER_POINTER(s_nXX);
  int *nSave = INTEGER_POINTER(s_nSave);
  int *nIter = INTEGER_POINTER(s_nIter);
  int *iToSave = INTEGER_POINTER(s_iToSave);
  double *ff = NUMERIC_POINTER(s_ff);
  double *meanVals = NUMERIC_POINTER(s_meanVals);

  SEXP s_xx;
  double xx[(*nXX)];

  mewInit(*nBin, *nXX, *ff, &av);

  for (k=0; k<(*nIter); k++) {

    PROTECT(s_xx = eval(f, rho));
    for (l=0; l<(*nXX); l++) {
      
      xx[l] = NUMERIC_POINTER(s_xx)[l];
    }
    UNPROTECT(1);

    mewAccum(xx, &av);

    if (iToSave[k] == 1) {
      
      mewMean(&av);
      totalSaved = totalSaved + 1;
      if (totalSaved <= (*nSave)) {
	
	for (l=0; l<(*nXX); l++) {

	  meanVals[(totalSaved - 1) + (*nSave)*l] = av.xMean[l];
	}
      }
      else {
	
	error("Attempting to save more mean values than specified!");
      }
    }
  }
  
  mewFinal(&av);

  UNPROTECT(7);

  return R_NilValue;
}
