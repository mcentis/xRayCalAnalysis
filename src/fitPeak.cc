#include "fitPeak.hh"

#include "TH1.h"
#include "TF1.h"
#include "TAxis.h"
#include "TCanvas.h"

TF1* fitPeak(TH1* hist, double negSig, double posSig)
{
  if(hist == 0) return 0;

  int maxBin = hist->GetMaximumBin();
  double maxPos = hist->GetXaxis()->GetBinCenter(maxBin);

  TCanvas* fitCan = new TCanvas("fitCan");

  TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetParameter(0, hist->GetMaximum());
  fitFunc->SetParameter(1, maxPos);
  fitFunc->SetParameter(2, hist->GetRMS());
  double start = maxPos - hist->GetRMS();
  double stop = maxPos + hist->GetRMS();
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  start = fitFunc->GetParameter(1) - negSig * fitFunc->GetParameter(2);
  stop = fitFunc->GetParameter(1) + posSig * fitFunc->GetParameter(2);
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  delete fitCan;

  return fitFunc;
}