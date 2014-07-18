/**
 * function to fit the peak of a histo
 * two fits are performed, a first to find the fit range, and a second to estimate the parameters
 */

#ifndef FITPEAK_HH
#define FITPEAK_HH

#include "TH1.h"
#include "TF1.h"

TF1* fitPeak(TH1* hist)
{
  TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetParameter(0, hist->GetMaximum());
  fitFunc->SetParameter(1, hist->GetMaximumX());
  fitFunc->SetParameter(2, hist->GetRMS());
  double start = hist->GetMaximumX() - hist->GetRMS();
  double stop = hist->GetMaximumX() + hist->GetRMS();
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  start = fitFunc->GetParameter(1) - fitFunc->GetParameter(2);
  stop = fitFunc->GetParameter(1) + fitFunc->GetParameter(2);
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  return fitFunc;
}

#endif // #ifndef FITPEAK_HH

