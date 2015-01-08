/*
 * program to produce a hitmap from the tree with the events
 */

#include "rocToModuleCR.hh"
#include "constants.h"

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

  char title[50];
  char name[50];

  TH2D* hitMap = new TH2D("hitMap", "Module hit map (if it is a module)", 416, -0.5, 415.5, 160, -0.5, 159.5);
  TH1I* pixEvt = new TH1I("pixEvt", "Fired pxels per event;Fired pix per event;Entries", 1000, -0.5, 999.5);
  TH1I* rocs = new TH1I("rocs", "Rocs where pix fired;Roc number;Entries", 16, -0.5, 15.5);

  TH2D* singleRocs[nRoc];
  for(int i = 0; i < nRoc; ++i)
    {
      sprintf(name, "hitmapRoc_%i", i);
      sprintf(title, "Hit map of Roc %i;Col;Row", i);
      singleRocs[i] = new TH2D(name, title, 52, -0.5, 51.5, 80, -0.5, 79.5);
    }

  TDirectory* dir;
  TTree* tree;
  // long int hits;
  long int events;
  long int firedPix = 0;
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

  TFile* inFile;
  inFile = TFile::Open(argv[1]);
  std::cout << "\t Input file : " << inFile->GetName() << std::endl;

  const char* fileName = "hitMapFromTree.root";
  TFile* outFile = new TFile(fileName, "RECREATE");
  std::cout << "\t Output file : " << outFile->GetName() << std::endl;

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

  int modCR[2] = {0}; // module column and row

  events = tree->GetEntries();

  for(long int i = 0; i < events; ++i)
    {
      tree->GetEntry(i);

      if(npix >= 2000)
	{
	  std::cout << "Event with " << npix << " pixels" << std::endl;
	  continue;
	}

      pixEvt->Fill(npix);

      firedPix += npix;

      for(int j = 0; j < npix; ++j)
	{
	  moduleColRow(proc[j], pcol[j], prow[j], modCR);
	  hitMap->Fill(modCR[0], modCR[1]);

	  singleRocs[proc[j]]->Fill(pcol[j], prow[j]);
 
	  rocs->Fill(proc[j]);
	}
    }

  std::cout << "Tot fired pixels " << firedPix << std::endl;
  std::cout << "Entries in hitmap " << hitMap->GetEntries() << std::endl;
  std::cout << "Rate module " << firedPix / (evtDuration * events * sensArea * 16) * 1e-3 << " kHz / cm2" << std::endl;

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
  // hitCan->cd(2);
  // nEntries->Draw("TEXTCOLZ");
  hitCan->cd(3);
  hitMap->Draw("COLZ");
  hitCan->cd(4);
  projY->SetFillColor(602);
  projY->Draw("HBAR");

  hitCan->Write();
  hitMap->Write();
  pixEvt->Write();
  rocs->Write();

  projX->SetFillColor(kWhite);
  projY->SetFillColor(kWhite);
  projX->Write();
  projY->Write();

  for(int i = 0; i < nRoc; ++i)
    singleRocs[i]->Write();

  outFile->Close();

  return 0;
}
