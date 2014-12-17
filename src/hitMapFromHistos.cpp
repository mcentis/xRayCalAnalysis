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
  TH2D* nEntries = new TH2D("nEntries", "Pixels fired in each ROC", 8, -0.5, 7.5, 2, -0.5, 1.5); // draw whith colztext to have numbers superimposed to the bins

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

      if(iRoc < 8)
	nEntries->SetBinContent(iRoc + 1, 1, hist->GetEntries());
      else
	nEntries->SetBinContent(16 - iRoc, 2, hist->GetEntries());

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

  // projections
  TH1D* projX = hitMap->ProjectionX("projX");
  projX->SetTitle("Projection of the hitmap on X;Col;Entries");

  TH1D* projY = hitMap->ProjectionY("projY");
  projY->SetTitle("Projection of the hitmap on Y;Row;Entries");

  TCanvas* hitCan = new TCanvas("hitMapCan", "Hit map");
  hitCan->Divide(2, 2);
  hitCan->cd(1);
  projX->SetFillColor(602);
  projX->Draw();
  hitCan->cd(2);
  nEntries->Draw("TEXTCOLZ");
  hitCan->cd(3);
  hitMap->Draw("COLZ");
  hitCan->cd(4);
  projY->SetFillColor(602);
  projY->Draw("HBAR");


  outFile->cd();

  hitCan->Write();
  hitMap->Write();
  nEntries->Write(); // draw whith colztext to have numbers superimposed to the bins

  projX->SetFillColor(kWhite);
  projY->SetFillColor(kWhite);
  projX->Write();
  projY->Write();
 
  outFile->Close();

  return 0;
}
