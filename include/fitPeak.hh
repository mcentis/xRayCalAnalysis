#ifndef FITPEAK_HH
#define FITPEAK_HH

#include "TH1.h"
#include "TF1.h"

class fitPeak
{
public:
  fitPeak(TH1* inHist);

private:
  TH1* hist;
  TF1* fitFunc;

};

#endif // #ifndef FITPEAK_HH

