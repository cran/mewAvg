// This is a Translation of Zachary Levine's mewAvg.f90 to C
// All of Zachary's comments are included in this code
// Unless a comment is specifically noted to be Adam's, it is
// Zachary's

// Zachary Levine 23 May 2013 - 4 June 2013
// The purpose of this module is to implement a particular averaging
// scheme which allows convergence in stochastic optimization.  The
// background is described in JC Spall, "Intro. to Stochastic Search
// and Optimization" Wiley, 2003, Chap. 4.

// The idea is to obtain the average from the following sum (for N
// even):

// \bar X = \lim_{N->\infty} {2/N} \sum_{i=(N/2)+1}^N X_i

// That is, we have a moving, expanding average, where the first half
// of the samples are discarded.  The first half are discarded because
// they are obtained under conditions which are different than those
// of the converged parameters. (The "half" is parameterized in the
// implementation.)

// In order to use a fixed amount of storage as N->\infty, we will
// have a fixed number of bins (nBin) which are partial sums of the
// series.  The number of samples in each bin increases exponentially
// (by a factor of ww, rounded to an integer).  The oldest bin is
// phased out as the newest bin is filled.

// To avoid keeping track of many shapes, the X_i is taken to be a 1D
// array.

// At the begining, only one sample is stored per bin until all bins
// have at least one sample.  At the very beginning, the mean is set
// to 0.

// Usage
//   loop over independent uses
//      call mewInit
//      loop over sample acquisition and use of mean
//         call mewAccum (when new data exists)
//         call mewMean  (whenever desired)
//      call mewFinal (optional - space reuseable in any case)


#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<R.h>
#include<Rdefines.h>
#include"mewAvgForRpkg.h"

// Adam comment
// This function initializes an instance of structrue type mewTyp
// There must always be a call to mewFinal when there
// is a call to mewInit
// if mewInit is successful, zero will be returned
// if mewInit encounters an error a different integer is returned
// check the source code for the meaning of the error codes
int mewInit (int nBin, 
	     int nXX,
	     double ff,
	     mewTyp *av) {

  int i, j;
  
  (av->nBin) = nBin;
  (av->nXX) = nXX;
  (av->ff) = ff;

  // checking inputs for errors

  if ((av->nBin) < 3) {
    
    error("averageInit:  need 3<=nBin but was %d\n", (av->nBin));
  }
  else if ((av->nBin) <= 0) {
    
    error("averageInit:  need nBin>0 but was %d\n", (av->nBin));
  }

  if (((av->ff) <= 0.0) || ((av->ff) >= 1.0)) {
    
    error("averageInit:  need 0<ff<1, but was %f\n", (av->ff));
  }

  if ((av->nXX) <= 0) {
    
    error("averageInit:  need nXX>0 but was %d\n", (av->nXX));
  }
  
  // increment factor
  (av->ww) = pow((1.0 - (av->ff)), (-1.0/(av->nBin))); 

  // ideally this many samples
  (av->aSample) = 1.0;

  // no bins in use so far
  (av->nBinUse) = 0;

  // first bin to write minus 1
  (av->iNew) = -1;

  // first bin to overwrite
  (av->iOld) = 0;

  // mean known
  (av->knowMean) = 1;

  // Adam comment
  // allocate memory for xMean (of length nBin)
  // initialize xMean
  (av->xMean) = Calloc((av->nXX), double);
  
  for (i=0; i<(av->nXX); i++) {

    (av->xMean)[i] = 0.0;
  }

  //Adam comment
  // allocate memory for xSumPart (of length nXX)
  // and initialize xSumPart
  (av->xSumPart) = Calloc((av->nXX), double);

  for (i=0; i<(av->nXX); i++) {
    
    (av->xSumPart)[i] = 0.0;
  }
  
  // Adam comment
  // allocate memory for xx (of dimension nXX rows by nBins cols)
  // and initialize xx
  (av->xx) = Calloc((av->nXX), double*);

  for (i=0; i<(av->nXX); i++) {
    
    (av->xx)[i] = Calloc((av->nBin), double); 
  }

  for (i=0; i<(av->nXX); i++) {
    
    for (j=0; j<(av->nBin); j++) {

      (av->xx)[i][j] = 0.0;
    }
  }

  // Adam comment
  // Allocate memory for mSample
  // and initialize mSample
  (av->mSample) = Calloc((av->nBin), int);

  for (i=0; i<(av->nBin); i++) {
    
    (av->mSample)[i] = 1;
  }

  // Adam comment
  // Allocate memory for nSample
  // and initialze nSample
  (av->nSample) = Calloc((av->nBin), int);

  for (i=0; i<(av->nBin); i++) {
    
    (av->nSample)[i] = 0;
  }

  return(0);
}

// Adam comment
// Whenever there is a call to mewInit, there must
// be a corresponding call to mewFinal
int mewFinal(mewTyp *av) {

  int i;

  Free((av->mSample));
  Free((av->nSample));
  Free((av->xMean));
  Free((av->xSumPart));
  for(i=0; i<(av->nXX); i++) {
    
    Free((av->xx)[i]);
  }
  Free((av->xx));

  return(0);
}

int zRound (double x) {
  
  int value;

  value = (int) floor((x + 0.5));

  return(value);
}

// ///////////////////////////////////////////////////////////////////
//    Accumulate samples for average
// (1) not all bins are in use yet; assign the new sample to the next
//     bin
// (2) all bins are in use
//     (2a) accumulate in existing bin
//     (2b) eliminate an old bin, and start a new one.  
//          Also, the partial sum of all the not-old-and-not-new 
//          bins is found.
//     (2c) error
// (3) error
int mewAccum (double *xx, 
	      mewTyp *av) {
  
  int i, j; 
  int *iNotNewNotOld;
  int iNew, iOld; // unpack to simplify syntax

  iNew = (av->iNew);
  iOld = (av->iOld);

  iNotNewNotOld = Calloc(((av->nBin) - 2), int);
  
  if ((av->nBinUse) < (av->nBin)) {
    // (1)
   
    iNew = iNew + 1;
    (av->nBinUse) = (av->nBinUse) + 1;
    for (i=0; i<(av->nXX); i++) {
      
      (av->xx)[i][iNew] = xx[i];
    }
    (av->nSample)[iNew] = 1;

    if ((av->nBinUse) == 1) {
      
      for (i=0; i<(av->nXX); i++) {
	
	(av->xMean)[i] = xx[i];
      }
    }
    else {
      
      for (i=0; i<(av->nXX); i++) {
	
	(av->xMean)[i] = (((av->nBinUse) - 1)*(av->xMean)[i] + xx[i])/
	  (av->nBinUse);
      }
    }
  }
  else if ((av->nBinUse) == (av->nBin)){
    // (2)
    
    (av->knowMean) = 0;
    
    if ((av->nSample)[iNew] < (av->mSample)[iNew]) {
      // (2a)
      
      (av->nSample)[iNew] = (av->nSample)[iNew] + 1;
      for (i=0; i<(av->nXX); i++) {
	
	(av->xx)[i][iNew] = (av->xx)[i][iNew] + xx[i];
      }
    }
    else if ((av->nSample)[iNew] == (av->mSample)[iNew]) {
      // (2b)
     
      iNew = iOld;
      
      if (iOld < ((av->nBin) - 1)) {
	
	iOld = iOld + 1;
      }
      else {

	iOld = 0;
      }

      for (i=0; i<(av->nXX); i++) {

	(av->xx)[i][iNew] = xx[i];
      }
      (av->nSample)[iNew] = 1;
      (av->aSample) = (av->aSample)*(av->ww);
      (av->mSample)[iNew] = zRound((av->aSample));
      
      j = -1;
      for (i=0; i<(av->nBin); i++) {
	
	if ((i == iOld) || (i == iNew))
	  continue;
	
	j = j + 1;
	iNotNewNotOld[j] = i;
      }

      (av->nPart) = 0;
      for (i=0; i<((av->nBin) - 2); i++) {
	
	(av->nPart) = (av->nPart) + (av->nSample)[iNotNewNotOld[i]];
      }

      for (i=0; i<(av->nXX); i++) {
	
	(av->xSumPart)[i] = 0.0;

	for (j=0; j<((av->nBin) - 2); j++) {
	  
	  (av->xSumPart)[i] = (av->xSumPart)[i] + 
	    (av->xx)[i][iNotNewNotOld[j]];
	}
      }
    }
    else {
      
      Rprintf("mewAccum:  likely programming error iNew is %d\n", 
	      iNew);
      Rprintf("mewAccum:  (av->nSample)[iNew] is %d\n", 
	      (av->nSample)[iNew]);
      Rprintf("mewAccum:  need <= \n");
      error("mewAccum:  (av->mSample)[iNew] which is %d\n", 
	    (av->mSample)[iNew]);
    }
  }
  else {
    
    Rprintf("mewAccum:  likely programming error\n");
    Rprintf("mewAccum:  (av->nBinUse) is %d\n", (av->nBinUse));
    Rprintf("mewAccum:  need <= \n");
    error("mewAccum:  (av->nBin) which is %d\n", (av->nBin));
  }

  (av->iNew) = iNew;
  (av->iOld) = iOld;

  Free(iNotNewNotOld);

  return(0);
}

// ///////////////////////////////////////////////////////////////////
// Writes the mean using a windowed moving average.
// The oldest data is reduced in weight linearly as the new data comes
// in.
int mewMean (mewTyp *av) {
  
  double fract;    // fraction of oldest bin to include
  int i;
  int iNew, iOld;  // unpack to simplify syntax

  if ((av->knowMean) == 1) {
    // if the mean is known, no action is needed
    
    return (0);
  }

  iNew = (av->iNew);
  iOld = (av->iOld);

  (av->knowMean) = 1;

  fract = ((av->mSample)[iNew] + 1.0 - (av->nSample)[iNew])/
    ((av->mSample)[iNew] + 1.0);

  for (i=0; i<(av->nXX); i++) {

    (av->xMean)[i] = ((av->xx)[i][iOld]*fract + 
		      (av->xSumPart)[i] + 
		      (av->xx)[i][iNew])/
      ((av->nSample)[iOld]*fract + 
       (av->nPart) + 
       (av->nSample)[iNew]);
  }

  return(0);
}
