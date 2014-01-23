typedef struct {
  
  int iNew;         // index of bin to add
  int iOld;         // index of bin to deweight (if 0, don't use)
  int knowMean;     // flag 0: mean not known 1: mean known 
  int nBin;         // number of groups of samples allocated
  int nBinUse;      // number of groups of samples actually defined
  int nXX;          // length of each sample
  int nPart;        // number of groups summed
  int *mSample;     // max allowed samples in bin (of length nBin)
  int *nSample;     // number of samples in each bin (of length nBin)
  double *xMean;    // mean of data (of length nXX)
  double *xSumPart; // sum of the part of the data which is not 
                    // changing
  double **xx;      // data (dimension is nXX rows by nBin cols)
  double ff;        // fraction of samples to retain in average 
                    // (typ. 1/2)
  double ww;        // factor of increase in number of samples
                    // from one bin to the next
  double aSample;   // ideal number of samples in a bin 
                    // (before rounding)
} mewTyp;

int mewInit (int nBin, 
	     int nXX,
	     double ff,
	     mewTyp *av);

int mewFinal(mewTyp *av);

int zRound (double x);

int mewAccum (double *xx, 
	      mewTyp *av);

int mewMean (mewTyp *av);
