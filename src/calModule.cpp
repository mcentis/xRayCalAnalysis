/*
 * program to calibrate a full module
 */

#include "constants.h"
#include "fitPeak.hh"

#include "iostream"
#include "fstream"

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TStyle.h"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "\tusage: calModule listFile dataDir" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  char fileName[200];
  char title[50];
  char name[50];
  char histName[50]; // name of the histo to look for

  sprintf(fileName, "%s/calModule.root", argv[2]);
  std::cout << "\t Output file : " << fileName << std::endl;
  TFile* outFile = new TFile(fileName, "RECREATE");

  TGraphErrors* calGr[nRoc];

  for(int i = 0; i < nRoc; ++i)
    {
      sprintf(name, "calRoc_%i", i);
      sprintf(title, "Calibration ROC %i", i);
      calGr[i] = new TGraphErrors();
      calGr[i]->SetName(name);
      calGr[i]->SetTitle(title);
      calGr[i]->SetMarkerStyle(22);
    }

  TFile* inFile;
  TDirectory* dir;
  TDirectory* elementDir;
  TF1* fit;
  TH1* hist;

  int nPoint = 0;
  const double negSigma = 0.5; // fit limits
  const double posSigma = 1.5;

  double ionization;
  std::string file;
  std::string element;

  std::ifstream listStr(argv[1], std::ifstream::in);
  if(listStr.is_open() == false)
    {
      std::cout << "\tERROR: could not open " << argv[1] << std::endl;
      return 1;
    }

  std::cout << "\tFiles being read:\n";
  std::cout << "\tIMPORTANT: if a file appears many times, correct the list file for newlines" << std::endl;

  while(listStr.good()) // loop on the files
    {
      listStr >> element >> ionization >> file;
      std::cout << '\t' << file << '\t' << element << '\t' << ionization << std::endl;
      sprintf(fileName, "%s/%s", argv[2], file.c_str());
      inFile = TFile::Open(fileName);

      if(inFile == 0)
	{
	  std::cout << "\tERROR: could not open " << fileName << std::endl;
	  continue;
	}

      dir = (TDirectory*) inFile->Get("Xray");
      if(dir == 0)
	{
	  std::cout << "\tERROR: TDirectory Xray not found" << std::endl;
	  continue;
	}

      sprintf(name, "histos_%s", element.c_str());
      elementDir = outFile->mkdir(name);

      for(int iRoc = 0; iRoc < nRoc; ++iRoc)
	{

	  sprintf(histName, "q_%s_C%i_V0", element.c_str(), iRoc);

	  hist = (TH1*) dir->Get(histName);
	  if(hist == 0)
	    {
	      std::cout << "\tWARNING: histogram " << histName << " not found in file " << inFile->GetPath() << std::endl;
	      continue;
	    }

	  fit = fitPeak(hist, negSigma, posSigma);

	  elementDir->cd();
	  hist->Write();

	  nPoint = calGr[iRoc]->GetN();
	  calGr[iRoc]->SetPoint(nPoint, fit->GetParameter(1), ionization);
	  calGr[iRoc]->SetPointError(nPoint, errIon, 0);

	} // loop on the roc histos in a file
    } // loop on the files being read

  TCanvas* serv = new TCanvas("serv", "serv");
  serv->cd();

  for(int i = 0; i < nRoc; i++)
    {
      calGr[i]->Draw("AP");
      calGr[i]->GetXaxis()->SetTitle("Measured ionization [Vcal]");
      calGr[i]->GetYaxis()->SetTitle("Expected ionization [e^{-}]");

      fit = new TF1("fitCal", "pol1", 0, 500);
      fit->SetParameters(0, 55); // arbitrary start values
      calGr[i]->Fit(fit, "RQ");
    }

  delete serv;

  outFile->cd();
  for(int i = 0; i < nRoc; i++)
    calGr[i]->Write();

  outFile->Close();

  return 0;
}
