/*
 * program to analyze the position of a peak as a function of the xray tube current
 */

#include "fitPeak.hh"
#include "clusterDefAndFunction.hh"

#include "TH1.h"
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
  if(argc != 4)
    {
      std::cout << "\tusage: peakVsRate element listFile dataDir" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  const double evtDuration = 25e-9; // 25 ns per event [s]
  const double sensArea = 0.6561; // single chip module area [cm^2]
  const double distCut = sqrt(2); // distance cut for clustering

  TFile* inFile;
  char fileName[200];
  char title[50];
  char name[50];
  char histName[50]; // name of the histo to look for
  sprintf(histName, "q_%s_C0_V0", argv[1]);

  sprintf(fileName, "%s/peakVsRate.root", argv[3]);
  std::cout << "\t Output file : " << fileName << std::endl;
  TFile* outFile = new TFile(fileName, "RECREATE");

  TGraphErrors* peakCurr = new TGraphErrors();
  peakCurr->SetName("peakCurr");
  sprintf(title, "%s line vs tube current", argv[1]);
  peakCurr->SetTitle(title);

  TGraphErrors* sigmaCurr = new TGraphErrors();
  sigmaCurr->SetName("sigmaCurr");
  sprintf(title, "#sigma %s line vs tube current", argv[1]);
  sigmaCurr->SetTitle(title);

  TGraphErrors* entriesCurr = new TGraphErrors();
  entriesCurr->SetName("entriesCurr");
  entriesCurr->SetTitle("Number of entries vs tube current");

  TGraphErrors* peakRate_xray = new TGraphErrors();
  peakRate_xray->SetName("peakRate_xray");
  sprintf(title, "%s line vs x-ray rate", argv[1]);
  peakRate_xray->SetTitle(title);

  TGraphErrors* peakRate_pix = new TGraphErrors();
  peakRate_pix->SetName("peakRate_pix");
  sprintf(title, "%s line vs pix rate", argv[1]);
  peakRate_pix->SetTitle(title);

  TGraphErrors* rate_xrayCurr = new TGraphErrors();
  rate_xrayCurr->SetName("rate_xrayCurr");
  sprintf(title, "Rate %s line vs tube current", argv[1]);
  rate_xrayCurr->SetTitle(title);

  TGraphErrors* rate_pixCurr = new TGraphErrors();
  rate_pixCurr->SetName("rate_pixCurr");
  sprintf(title, "Rate %s line (pixels) vs tube current", argv[1]);
  rate_pixCurr->SetTitle(title);

  int nPoint = 0;

  TDirectory* dir;
  TTree* tree;
  TF1* fit;
  TH1* hist;
  TH1D* histClustered;
  TH1I* clustSize;
  long int firedPix;
  long int hits;
  long int events;
  double ratePix;
  double rateXrays;
  double ratePixErr;
  double rateXraysErr;
  std::vector<cluster> cluVec;

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

  float curr;
  std::string file;

  const double negSigma = 1; // fit limits
  const double posSigma = 1.5;

  std::ifstream listStr(argv[2], std::ifstream::in);
  if(listStr.is_open() == false)
    {
      std::cout << "\tERROR: could not open " << argv[2] << std::endl;
      return 1;
    }

  std::cout << "\tFiles being read:\n";
  std::cout << "\tIMPORTANT: if a file appears many times, correct the list file for newlines, multiple entries with same current cause crashes" << std::endl;

  while(listStr.good()) // loop on the files
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

      fit = fitPeak(hist, negSigma, posSigma);

      nPoint = peakCurr->GetN();
      peakCurr->SetPoint(nPoint, curr, fit->GetParameter(1));
      peakCurr->SetPointError(nPoint, 0.01, fit->GetParError(1));

      sigmaCurr->SetPoint(nPoint, curr, fit->GetParameter(2));
      sigmaCurr->SetPointError(nPoint, 0, fit->GetParError(2));

      entriesCurr->SetPoint(nPoint, curr, hist->GetEntries());
      entriesCurr->SetPointError(nPoint, 0.01, sqrt(hist->GetEntries()));

      sprintf(name, "%s_%.00fuA", argv[1], curr * 1000);
      hist->SetName(name);

      outFile->cd();
      hist->Write();

      tree = (TTree*) dir->Get("events");
      if(!tree)
	{
	  std::cout << "\tERROR: Can not find TTree events" << std::endl;
	  continue;
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

      firedPix = 0;
      hits = 0;

      sprintf(name, "%s_%.00fuA_clustered", argv[1], curr * 1000);
      histClustered = new TH1D(name, "Clustered charge;Charge [Vcal];Entries / bin", 301, -0.5, 300.5);

      sprintf(name, "%s_%.00fuA_clulsterSize", argv[1], curr * 1000);
      clustSize = new TH1I(name, "Cluster size (zero suppressed);Cluster size [Pixel];Entries / bin", 11, -0.5, 10.50);

      events = tree->GetEntries();
      for(long int i = 0; i < events; ++i)
	{
	  tree->GetEntry(i);
	  firedPix += npix;

	  cluVec = clusterEvent(npix, pcol, prow, pq, distCut);

	  hits += cluVec.size();

	  for(std::vector<cluster>::iterator it = cluVec.begin(); it != cluVec.end(); ++it)
	    {
	      histClustered->Fill(it->totCharge);

	      if(it->nPix != 0) clustSize->Fill(it->nPix);
	    }
	}

      fitPeak(histClustered, negSigma, posSigma);

      histClustered->Write();
      clustSize->Write();

      ratePix = firedPix / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]
      rateXrays = hits / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]

      ratePixErr = sqrt(firedPix) / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]
      rateXraysErr = sqrt(hits) / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]

      nPoint = peakRate_xray->GetN();

      peakRate_xray->SetPoint(nPoint, rateXrays, fit->GetParameter(1));
      peakRate_xray->SetPointError(nPoint, rateXraysErr, fit->GetParError(1));

      peakRate_pix->SetPoint(nPoint, ratePix, fit->GetParameter(1));
      peakRate_pix->SetPointError(nPoint, ratePixErr, fit->GetParError(1));

      rate_xrayCurr->SetPoint(nPoint, curr, rateXrays);
      rate_xrayCurr->SetPointError(nPoint, 0.01, rateXraysErr);

      rate_pixCurr->SetPoint(nPoint, curr, ratePix);
      rate_pixCurr->SetPointError(nPoint, 0.01, ratePixErr);

    }

  listStr.close();

  TCanvas* peakCurrCan = new TCanvas("peakCurrCan");
  peakCurrCan->SetGridx();
  peakCurrCan->SetGridy();
  peakCurr->SetMarkerStyle(20);
  peakCurr->SetMarkerSize(2);
  peakCurr->Draw("AP");
  peakCurr->GetXaxis()->SetTitle("Tube current [mA]");
  peakCurr->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakCurr->Write();
  peakCurrCan->Modified();
  peakCurrCan->Update();
  peakCurrCan->Write();

  TCanvas* sigmaCurrCan = new TCanvas("sigmaCurrCan");
  sigmaCurrCan->SetGridx();
  sigmaCurrCan->SetGridy();
  sigmaCurr->SetMarkerStyle(20);
  sigmaCurr->SetMarkerSize(2);
  sigmaCurr->Draw("AP");
  sigmaCurr->GetXaxis()->SetTitle("Tube current [mA]");
  sigmaCurr->GetYaxis()->SetTitle("Sigma [Vcal]");
  sigmaCurr->Write();
  sigmaCurrCan->Modified();
  sigmaCurrCan->Update();
  sigmaCurrCan->Write();

  TCanvas* entriesCurrCan = new TCanvas("entriesCurrCan");
  entriesCurrCan->SetGridx();
  entriesCurrCan->SetGridy();
  entriesCurr->SetMarkerStyle(22);
  entriesCurr->SetMarkerSize(2);
  entriesCurr->Draw("AP");
  entriesCurr->GetXaxis()->SetTitle("Tube current [mA]");
  entriesCurr->GetYaxis()->SetTitle("Entries");
  entriesCurr->Write();
  entriesCurrCan->Modified();
  entriesCurrCan->Update();
  entriesCurrCan->Write();

  TCanvas* peakRatePixCan = new TCanvas("peakRatePixCan");
  peakRatePixCan->SetGridx();
  peakRatePixCan->SetGridy();
  peakRate_pix->SetMarkerStyle(20);
  peakRate_pix->SetMarkerSize(2);
  peakRate_pix->Draw("AP");
  peakRate_pix->GetXaxis()->SetTitle("Rate fired pixels [kHz cm^{-2}]");
  peakRate_pix->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakRate_pix->Write();
  peakRatePixCan->Modified();
  peakRatePixCan->Update();
  peakRatePixCan->Write();

  TCanvas* peakRateXrayCan = new TCanvas("peakRateXrayCan");
  peakRateXrayCan->SetGridx();
  peakRateXrayCan->SetGridy();
  peakRate_xray->SetMarkerStyle(20);
  peakRate_xray->SetMarkerSize(2);
  peakRate_xray->Draw("AP");
  peakRate_xray->GetXaxis()->SetTitle("Rate x-rays [kHz cm^{-2}]");
  peakRate_xray->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakRate_xray->Write();
  peakRateXrayCan->Modified();
  peakRateXrayCan->Update();
  peakRateXrayCan->Write();

  TCanvas* xrayCurrCan = new TCanvas("xrayCurrCan");
  xrayCurrCan->SetGridx();
  xrayCurrCan->SetGridy();
  rate_xrayCurr->SetMarkerStyle(22);
  rate_xrayCurr->SetMarkerSize(2);
  rate_xrayCurr->Draw("AP");
  rate_xrayCurr->GetXaxis()->SetTitle("Tube current [mA]");
  rate_xrayCurr->GetYaxis()->SetTitle("Rate x-rays [kHz cm^{-2}]");
  rate_xrayCurr->Write();
  xrayCurrCan->Modified();
  xrayCurrCan->Update();
  xrayCurrCan->Write();

  TCanvas* pixCurrCan = new TCanvas("pixCurrCan");
  pixCurrCan->SetGridx();
  pixCurrCan->SetGridy();
  rate_pixCurr->SetMarkerStyle(22);
  rate_pixCurr->SetMarkerSize(2);
  rate_pixCurr->Draw("AP");
  rate_pixCurr->GetXaxis()->SetTitle("Tube current [mA]");
  rate_pixCurr->GetYaxis()->SetTitle("Rate fired pixels [kHz cm^{-2}]");
  rate_pixCurr->Write();
  pixCurrCan->Modified();
  pixCurrCan->Update();
  pixCurrCan->Write();

  outFile->Close();

  return 0;
}
