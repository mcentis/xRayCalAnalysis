/**
 * function to fit the peak of a histo
 * two fits are performed, a first to find the fit range, and a second to estimate the parameters
 * posSig and negSig steer the range
 */

#ifndef FITPEAK_HH
#define FITPEAK_HH

#include "TH1.h"
#include "TF1.h"

TF1* fitPeak(TH1* hist, double negSig = 1, double posSig = 1);

#endif // #ifndef FITPEAK_HH

