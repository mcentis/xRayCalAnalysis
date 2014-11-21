#include "fitPeak.hh"

#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TAxis.h"
#include "TCanvas.h"

TF1* fitPeak(TH1* hist, double negSig, double posSig)
{
  if(hist == 0) return 0;

  TCanvas* fitCan = new TCanvas("fitCan");

  TSpectrum* spec = new TSpectrum(4, 1); // find the peaks (max 4) in the histo
  int nPeaks = spec->Search(hist, 10); // sigma peak = 10
  Double_t* pos = spec->GetPositionX();
  Double_t* heigth = spec->GetPositionY();

  double maxPos = -1;
  double startHeigth = -1;
  for(int i = 0; i < nPeaks; ++i)
    if(pos[i] > maxPos)
      {
	maxPos = pos[i];
	startHeigth = heigth[i];
      }

  //  std::cout << "Fitting around " << maxPos << std::endl;

  TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetParameter(0, startHeigth);
  fitFunc->SetParameter(1, maxPos);
  // fitFunc->SetParameter(2, hist->GetRMS());
  // double start = maxPos - hist->GetRMS();
  // double stop = maxPos + hist->GetRMS();
  fitFunc->SetParameter(2, 20);
  double start = maxPos - 20;
  double stop = maxPos + 20;
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  start = fitFunc->GetParameter(1) - negSig * fitFunc->GetParameter(2);
  stop = fitFunc->GetParameter(1) + posSig * fitFunc->GetParameter(2);
  fitFunc->SetRange(start, stop);
  hist->Fit(fitFunc, "RQ");

  delete fitCan;

  return fitFunc;
}
