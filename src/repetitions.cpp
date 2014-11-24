/*
 * program to analyze the position of a peak in repeated measurements
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
      std::cout << "\tusage: repetitions element listFile dataDir" << std::endl;
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

  sprintf(fileName, "%s/repetitions.root", argv[3]);
  std::cout << "\t Output file : " << fileName << std::endl;
  TFile* outFile = new TFile(fileName, "RECREATE");

  TGraphErrors* peakRep = new TGraphErrors();
  peakRep->SetName("peakRep");
  sprintf(title, "%s line vs repetition", argv[1]);
  peakRep->SetTitle(title);

  TGraphErrors* sigmaRep = new TGraphErrors();
  sigmaRep->SetName("sigmaRep");
  sprintf(title, "#sigma %s line vs repetition", argv[1]);
  sigmaRep->SetTitle(title);

  TGraphErrors* entriesRep = new TGraphErrors();
  entriesRep->SetName("entriesRep");
  entriesRep->SetTitle("Number of entries vs repetition");

  TGraphErrors* peakRate_xray = new TGraphErrors();
  peakRate_xray->SetName("peakRate_xray");
  sprintf(title, "%s line vs x-ray rate", argv[1]);
  peakRate_xray->SetTitle(title);

  TGraphErrors* peakRate_pix = new TGraphErrors();
  peakRate_pix->SetName("peakRate_pix");
  sprintf(title, "%s line vs pix rate", argv[1]);
  peakRate_pix->SetTitle(title);

  TGraphErrors* rate_xrayRep = new TGraphErrors();
  rate_xrayRep->SetName("rate_xrayRep");
  sprintf(title, "Rate %s line vs repetition", argv[1]);
  rate_xrayRep->SetTitle(title);

  TGraphErrors* rate_pixRep = new TGraphErrors();
  rate_pixRep->SetName("rate_pixRep");
  sprintf(title, "Rate %s line (pixels) vs tube repetition", argv[1]);
  rate_pixRep->SetTitle(title);

  TH1D* peakDistr = new TH1D("peakDistr", "Distribution of the peaks;Charge [Vcal];Entries / bin", 301, -0.5, 300.5);
  TH1D* sigmaDistr = new TH1D("sigmaDistr", "Distribution of the peaks sigma;Charge [Vcal];Entries / bin", 36, -0.5, 35.5);

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

  float rep;
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
  std::cout << "\tIMPORTANT: if a file appears many times, correct the list file for newlines" << std::endl;

  while(listStr.good()) // loop on the files
    {
      listStr >> rep >> file;
      std::cout << '\t' << file << '\t' << rep << std::endl;
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

      nPoint = peakRep->GetN();
      peakRep->SetPoint(nPoint, rep, fit->GetParameter(1));
      peakRep->SetPointError(nPoint, 0.01, fit->GetParError(1));

      sigmaRep->SetPoint(nPoint, rep, fit->GetParameter(2));
      sigmaRep->SetPointError(nPoint, 0, fit->GetParError(2));

      entriesRep->SetPoint(nPoint, rep, hist->GetEntries());
      entriesRep->SetPointError(nPoint, 0.01, sqrt(hist->GetEntries()));

      peakDistr->Fill(fit->GetParameter(1));
      sigmaDistr->Fill(fit->GetParameter(2));

      sprintf(name, "%s_%.00f", argv[1], rep);
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

      sprintf(name, "%s_%.00f_clustered", argv[1], rep);
      histClustered = new TH1D(name, "Clustered charge;Charge [Vcal];Entries / bin", 301, -0.5, 300.5);

      sprintf(name, "%s_%.00f_clulsterSize", argv[1], rep);
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

      rate_xrayRep->SetPoint(nPoint, rep, rateXrays);
      rate_xrayRep->SetPointError(nPoint, 0.01, rateXraysErr);

      rate_pixRep->SetPoint(nPoint, rep, ratePix);
      rate_pixRep->SetPointError(nPoint, 0.01, ratePixErr);

    }

  listStr.close();

  TCanvas* peakRepCan = new TCanvas("peakRepCan");
  peakRepCan->SetGridx();
  peakRepCan->SetGridy();
  peakRep->SetMarkerStyle(20);
  peakRep->SetMarkerSize(2);
  peakRep->Draw("AP");
  peakRep->GetXaxis()->SetTitle("Repetition");
  peakRep->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakRep->Write();
  peakRepCan->Modified();
  peakRepCan->Update();
  peakRepCan->Write();

  TCanvas* sigmaRepCan = new TCanvas("sigmaRepCan");
  sigmaRepCan->SetGridx();
  sigmaRepCan->SetGridy();
  sigmaRep->SetMarkerStyle(20);
  sigmaRep->SetMarkerSize(2);
  sigmaRep->Draw("AP");
  sigmaRep->GetXaxis()->SetTitle("Repetition");
  sigmaRep->GetYaxis()->SetTitle("Sigma [Vcal]");
  sigmaRep->Write();
  sigmaRepCan->Modified();
  sigmaRepCan->Update();
  sigmaRepCan->Write();

  TCanvas* entriesRepCan = new TCanvas("entriesRepCan");
  entriesRepCan->SetGridx();
  entriesRepCan->SetGridy();
  entriesRep->SetMarkerStyle(22);
  entriesRep->SetMarkerSize(2);
  entriesRep->Draw("AP");
  entriesRep->GetXaxis()->SetTitle("Repetition");
  entriesRep->GetYaxis()->SetTitle("Entries");
  entriesRep->Write();
  entriesRepCan->Modified();
  entriesRepCan->Update();
  entriesRepCan->Write();

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

  TCanvas* xrayRepCan = new TCanvas("xrayRepCan");
  xrayRepCan->SetGridx();
  xrayRepCan->SetGridy();
  rate_xrayRep->SetMarkerStyle(22);
  rate_xrayRep->SetMarkerSize(2);
  rate_xrayRep->Draw("AP");
  rate_xrayRep->GetXaxis()->SetTitle("Repetition");
  rate_xrayRep->GetYaxis()->SetTitle("Rate x-rays [kHz cm^{-2}]");
  rate_xrayRep->Write();
  xrayRepCan->Modified();
  xrayRepCan->Update();
  xrayRepCan->Write();

  TCanvas* pixRepCan = new TCanvas("pixRepCan");
  pixRepCan->SetGridx();
  pixRepCan->SetGridy();
  rate_pixRep->SetMarkerStyle(22);
  rate_pixRep->SetMarkerSize(2);
  rate_pixRep->Draw("AP");
  rate_pixRep->GetXaxis()->SetTitle("Repetition");
  rate_pixRep->GetYaxis()->SetTitle("Rate fired pixels [kHz cm^{-2}]");
  rate_pixRep->Write();
  pixRepCan->Modified();
  pixRepCan->Update();
  pixRepCan->Write();

  peakDistr->Write();
  sigmaDistr->Write();

  outFile->Close();

  return 0;
}
