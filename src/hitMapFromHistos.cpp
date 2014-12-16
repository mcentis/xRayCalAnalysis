/*
 * program to produce a hitmap of a module starting from the single hitmaps saved by pxar
 */

#include "rocToModuleCR.hh"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TStyle.h"

#include "iostream"
#include "fstream"
#include "math.h"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "\tusage: hitMapFromHistos element file" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  TH2* hist;
  TDirectory* dir;
  char histName[50]; // name of the histo to look for

  const char* fileName = "hitMapFromHistos.root";
  TFile* outFile = new TFile(fileName, "RECREATE");
  std::cout << "\t Output file : " << outFile->GetName() << std::endl;

  TH2D* hitMap = new TH2D("hitMap", "Module hit map (if it is a module)", 416, -0.5, 415.5, 160, -0.5, 159.5);

  TFile* inFile;
  inFile = TFile::Open(argv[2]);
  std::cout << "\t Input file : " << inFile->GetName() << std::endl;

  if(inFile->IsZombie())
    {
      std::cout << "\tERROR: could not open " << argv[1] << std::endl;
      return -1;
    }

  dir = (TDirectory*) inFile->Get("Xray");
  if(dir == 0)
    {
      std::cout << "\tERROR: TDirectory Xray not found" << std::endl;
      return -1;
    }

  int modCR[2] = {0};
  double content  = 0;

  for(int iRoc = 0; iRoc < 16; ++iRoc)
    {
      sprintf(histName, "hMap_%s_C%i_V0", argv[1], iRoc);
      hist = (TH2*) dir->Get(histName);

      if(hist == 0)
	{
	  std::cout << "\tWARNING Histogram " << histName << " not found" << std::endl;
	  continue;
	}

      std::cout << "\tINFO Found histogram " << histName << std::endl;

      for(int iColBin = 1; iColBin < nCol + 1; iColBin++)
	for(int iRowBin = 1; iRowBin < nRow + 1; iRowBin++)
	  {
	    // the + and - 1 in this scope are due to bin to pix conversion
	    moduleColRow(iRoc, iColBin - 1, iRowBin - 1, modCR);
	    content = hist->GetBinContent(iColBin, iRowBin);
	    if(iRoc < 8)
	      hitMap->SetBinContent(modCR[0] + 1, modCR[1] + 1, content);
	    else
	      hitMap->SetBinContent(modCR[0], modCR[1], content);
	  }
    }

  inFile->Close();
  outFile->cd();
  hitMap->Write();
  outFile->Close();

  return 0;
}
