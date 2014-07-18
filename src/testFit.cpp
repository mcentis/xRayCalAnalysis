/*
 * program to test the histo fitting
 */

#include "fitPeak.hh"

#include "TSystem.h"
#include "TFile.h"
#include "iostream"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "\tusage: testFit infile outfile" << std::endl;
      return 1;
    }
  gSystem->Load("libTree"); // to avoid warnings

  TFile* inFile = TFile::Open(argv[1]);
  TFile* outFile = new TFile(argv[2], "RECREATE");

  TDirectory* dir = (TDirectory*) inFile->Get("Xray");

  const int nHi = 3;
  const char* hiNames[nHi] = {"q_Cu_C0_V0", "q_Mo_C0_V0", "q_Ag_C0_V0"};

  TF1* fit;

  for(int i = 0; i < nHi; ++i)
    {
      TH1* hist = (TH1*) dir->Get(hiNames[i]);
      if(hist == 0)
	{
	  std::cout << "\thistogram " << hiNames[i] << "not found" << std::endl;
	  continue;
	}

      fit = fitPeak(hist);
      std::cout << hiNames[i] << " peak -> " << fit->GetParameter(1) << " +- " << fit->GetParError(1) << std::endl;

      hist->Write();
    }

  inFile->Close();
  outFile->Close();

  return 0;
}
