/*
 * program to produce a hitmap from the tree with the events
 */

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TStyle.h"

#include "iostream"
#include "fstream"
#include "math.h"

void moduleColRow(UChar_t roc, UChar_t col, UChar_t row, UChar_t* modCR) // returns module col and row
{
  modCR[0] = -1;
  modCR[1] = -1;

  const UChar_t nCol = 52; // number of columns in a roc
  const UChar_t nRow = 80; // number of rows in a roc

  if(roc < 8)
    {
      modCR[0] = roc * nCol + col;
      modCR[1] = row;
    }
  else
    {
      modCR[0] = 8 * nCol - (roc - 8) * nCol - col;
      modCR[1] = 2 * nRow - row;
    }

  return;
}

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      std::cout << "\tusage: hitMapFromTree file" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  // const double evtDuration = 25e-9; // 25 ns per event [s]
  // const double sensArea = 0.6561; // single chip module area [cm^2]

  const char* fileName = "hitMapFromTree.root";
  TFile* outFile = new TFile(fileName, "RECREATE");
  std::cout << "\t Output file : " << outFile->GetName() << std::endl;

  char title[50];
  char name[50];
  // char histName[50]; // name of the histo to look for
  // sprintf(histName, "q_%s_C0_V0", argv[1]);

  TH2D* hitMap = new TH2D("hitMap", "Module hit map (if it is a module)", 416, -0.5, 415.5, 160, -0.5, 159.5);
  TH1I* pixEvt = new TH1I("pixEvt", "Fired pxels per event;Fired pix per event;Entries", 1000, -0.5, 999.5);
  TH1I* rocs = new TH1I("rocs", "Rocs where pix fired;Roc number;Entries", 16, -0.5, 15.5);

  TH2D* singleRocs[16];
  for(int i = 0; i < 16; ++i)
    {
      sprintf(name, "hitmapRoc_%i", i);
      sprintf(title, "Hit map of Roc %i;Col;Row", i);
      singleRocs[i] = new TH2D(name, title, 52, -0.5, 51.5, 80, -0.5, 79.5);
    }

  std::cout << "created histos" << std::endl;

  TDirectory* dir;
  TTree* tree;
  // long int firedPix;
  // long int hits;
  long int events;
  // double ratePix;
  // double rateXrays;
  // double ratePixErr;
  // double rateXraysErr;

  //Declaration of leaves types
  const int maxPix = 2000; // from struct TreeEvent, PixTest.hh, pXar
  UShort_t        header;
  UShort_t        trailer;
  UShort_t        npix;
  UChar_t         proc[maxPix];
  UChar_t         pcol[maxPix];
  UChar_t         prow[maxPix];
  Double_t        pval[maxPix];
  Double_t        pq[maxPix];

  // UChar_t         pcolMod[maxPix]; // module col and row
  // UChar_t         prowMod[maxPix];

  TFile* inFile;
  inFile = TFile::Open(argv[1]);
  std::cout << "\t Input file : " << inFile->GetName() << std::endl;

  if(inFile == 0)
    {
      std::cout << "\tERROR: could not open " << fileName << std::endl;
      return -1;
    }
  
  dir = (TDirectory*) inFile->Get("Xray");
  if(dir == 0)
    {
      std::cout << "\tERROR: TDirectory Xray not found" << std::endl;
      return -1;
    }

  tree = (TTree*) dir->Get("events");
  if(!tree)
    {
      std::cout << "\tERROR: Can not find TTree events" << std::endl;
      return -1;
    }

  // Set branch addresses.
  tree->SetBranchAddress("header",&header);
  tree->SetBranchAddress("trailer",&trailer);
  tree->SetBranchAddress("npix",&npix);
  tree->SetBranchAddress("proc",proc);
  tree->SetBranchAddress("pcol",pcol);
  tree->SetBranchAddress("prow",prow);
  tree->SetBranchAddress("pval",pval);
  tree->SetBranchAddress("pq",pq);

  UChar_t modCR[2] = {0}; // module column and row

  events = tree->GetEntries();

  std::cout << events << std::endl;

  for(long int i = 0; i < events; ++i)
    {
      tree->GetEntry(i);

      pixEvt->Fill(npix);

      for(int j = 0; j < npix; ++j)
	{
	  moduleColRow(proc[j], pcol[j], prow[j], modCR);
	  hitMap->Fill(modCR[0], modCR[1]);
	  // pcolMod[j] = modCR[0];
	  // prowMod[j] = modCR[1];

	  if(proc[j] < 16 && proc[j] >= 0)
	    singleRocs[proc[j]]->Fill(pcol[j], prow[j]);
	  else
	    std::cout << "\t ERROR: ROC number not good!!! Got " << (int) proc[j] << std::endl;
 
	  rocs->Fill(proc[j]);
	}
    }

  std::cout << "out loop " << std::endl;

  inFile->Close();

  std::cout << "closed in " << std::endl;

  outFile->cd();

  std::cout << "cd out " << std::endl;

  hitMap->Write();
  pixEvt->Write();
  rocs->Write();

  for(int i = 0; i < 16; ++i)
    singleRocs[i]->Write();

  std::cout << "wrote histo " << std::endl;

  outFile->Close();

  std::cout << "closed out " << std::endl;

  return 0;
}
