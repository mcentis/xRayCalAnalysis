/*
 * program to analyze the position of a peak as a function of the xray tube current
 */

#include "fitPeak.hh"

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "iostream"
#include "fstream"

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
      std::cout << "\tusage: peakVsCurr element listFile dataDir" << std::endl;
      return 1;
    }

  gSystem->Load("libTree"); // to avoid warnings

  TFile* inFile;
  char fileName[200];
  char title[50];
  char name[50];
  char histName[50]; // name of the histo to look for
  sprintf(histName, "q_%s_C0_V0", argv[1]);

  sprintf(fileName, "%s/peakVsCurr.root", argv[3]);
  TFile* outFile = new TFile(fileName, "RECREATE");

  TGraphErrors* peakCurr = new TGraphErrors();
  peakCurr->SetName("peakCurr");
  sprintf(title, "%s line vs tube current", argv[1]);
  peakCurr->SetTitle(title);
  int nPoint = 0;

  TDirectory* dir;
  TF1* fit;
  TH1* hist;

  float curr;
  std::string file;

  std::ifstream listStr(argv[2], std::ifstream::in);
  if(listStr.is_open() == false)
    {
      std::cout << "\tERROR: could not open " << argv[2] << std::endl;
      return 1;
    }

  std::cout << "\tFiles being read:\n";
  std::cout << "\tIMPORTANT: if a file appears many times, correct the list file for newlines, multiple entries with same current cause crashes" << std::endl;

  while(listStr.good())
    {
      listStr >> curr >> file;
      std::cout << '\t' << file << '\t' << curr << std::endl;
      sprintf(fileName, "%s/%s", argv[3], file.c_str());
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

      hist = (TH1*) dir->Get(histName);
      if(hist == 0)
      	{
      	  std::cout << "\tWARNING: histogram " << histName << " not found in file " << inFile->GetPath() << std::endl;
      	  continue;
      	}

      fit = fitPeak(hist, 1, 1.5);

      nPoint = peakCurr->GetN();
      peakCurr->SetPoint(nPoint, curr, fit->GetParameter(1));
      peakCurr->SetPointError(nPoint, 0, fit->GetParError(1));

      sprintf(name, "%s_%.02fmA", argv[1], curr);
      hist->SetName(name);

      outFile->cd();
      hist->Write();
    }

  listStr.close();

  TCanvas* peakCurrCan = new TCanvas("peakCurrCan");
  peakCurrCan->SetGridx();
  peakCurrCan->SetGridy();
  peakCurr->Draw("AP");
  peakCurr->GetXaxis()->SetTitle("Tube current [mA]");
  peakCurr->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakCurr->Write();
  peakCurrCan->Modified();
  peakCurrCan->Update();
  peakCurrCan->Write();

  outFile->Close();

  return 0;
}
