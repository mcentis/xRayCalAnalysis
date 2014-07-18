#include "calSingleT.hh"
#include "constants.h"
#include "fitPeak.hh"

#include "stdlib.h"

#include "TF1.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"

calSingleT::calSingleT(const char* inFilePath, float temperature)
{
  temp = temperature;
  inFile = TFile::Open(inFilePath);
  if(inFile == 0)
    {
      std::cout << "\tERROR: could not open " << inFilePath << std::endl;
      exit(1);
    }

  dir = (TDirectory*) inFile->Get("Xray");
  if(dir == 0)
    {
      std::cout << "\tERROR: TDirectory Xray not found" << std::endl;
      exit(1);
    }

  calFit = new TF1("calFit", "pol1");
  calGr = new TGraphErrors();
  calGr->SetName("calGr");
  char title[200];
  sprintf(title, "Calibration %.00f C", temp);
  calGr->SetTitle(title);

  return;
}

calSingleT::~calSingleT()
{
  inFile->Close();

  return;
}

void calSingleT::DoCal()
{
  TF1* fit;
  TH1* hist;
  int nPoint;

  for(int i = 0; i < nHist; ++i)
    {
      hist = (TH1*) dir->Get(histNames[i]);
      if(hist == 0)
      	{
      	  std::cout << "\tWARNING: histogram " << histNames[i] << " not found in file " << inFile->GetPath() << std::endl;
      	  continue;
      	}

      fit = fitPeak(hist, 1, 1.5);
      histos.push_back(hist);

      nPoint = calGr->GetN();
      calGr->SetPoint(nPoint, fit->GetParameter(1), peakIon[i]);
      calGr->SetPointError(nPoint, fit->GetParError(1), 0);
    }

  if(calGr->GetN() == 0)
    {
      std::cout << "\tWARNING: the calibration graph has no points, from file " << inFile->GetPath() << std::endl;
      return;
    }

  TCanvas* can = new TCanvas("can");
  calGr->Draw("AP");
  calGr->GetXaxis()->SetTitle("Measured ionization [Vcal]");
  calGr->GetYaxis()->SetTitle("Teoretical charge deposit [e]");

  calFit->SetRange(calGr->GetXaxis()->GetXmin(), calGr->GetXaxis()->GetXmax());
  calFit->SetParameters(0, 55); // start parameters, quite arbitrary
  calGr->Fit(calFit, "RQ");

  delete can;

  return;
}
