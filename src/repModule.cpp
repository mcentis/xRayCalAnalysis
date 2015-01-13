/*
 * program to analyze the position of a peak in repeated measurements
 */

#include "fitPeak.hh"
#include "constants.h"

#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
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
      std::cout << "\tusage: repModule element listFile dataDir" << std::endl;
      return 1;
    }

  // big axis labels
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");

  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");

  gStyle->SetTitleOffset(0.95, "x");
  gStyle->SetTitleOffset(0.95, "y");

  TFile* inFile;
  char fileName[200];
  char title[50];
  char name[50];
  char histName[50]; // name of the histo to look for
  long int totEntries;

  sprintf(fileName, "%s/repModule.root", argv[3]);
  std::cout << "\t Output file : " << fileName << std::endl;
  TFile* outFile = new TFile(fileName, "RECREATE");

  TGraphErrors* peakRepRoc[nRoc];
  TGraphErrors* sigmaRepRoc[nRoc];
  TGraphErrors* peakRate_pixRoc[nRoc];
  TGraphErrors* rate_pixRepRoc[nRoc];

  int iColor;
  int iMarker;

  for(int iRoc = 0; iRoc < nRoc; ++iRoc)
    {
      iColor = iRoc % 8 + 1;
      if(iColor >= 5) iColor++; // skip yellow

      if(iRoc < 8)
	iMarker = 24; // empty circle
      else
	iMarker = 25; // empty square

      peakRepRoc[iRoc] = new TGraphErrors();
      sprintf(name, "peakRep_roc_%i", iRoc);
      peakRepRoc[iRoc]->SetName(name);
      sprintf(title, "%s line, ROC %i", argv[1], iRoc);
      peakRepRoc[iRoc]->SetTitle(title);
      peakRepRoc[iRoc]->SetLineColor(iColor);
      peakRepRoc[iRoc]->SetFillColor(kWhite);
      peakRepRoc[iRoc]->SetMarkerColor(iColor);
      peakRepRoc[iRoc]->SetMarkerStyle(iMarker);

      sigmaRepRoc[iRoc] = new TGraphErrors();
      sprintf(name, "sigmaRep_roc_%i", iRoc);
      sigmaRepRoc[iRoc]->SetName(name);
      sprintf(title, "#sigma %s line, ROC %i", argv[1], iRoc);
      sigmaRepRoc[iRoc]->SetTitle(title);
      sigmaRepRoc[iRoc]->SetLineColor(iColor);
      sigmaRepRoc[iRoc]->SetFillColor(kWhite);
      sigmaRepRoc[iRoc]->SetMarkerColor(iColor);
      sigmaRepRoc[iRoc]->SetMarkerStyle(iMarker);

      peakRate_pixRoc[iRoc] = new TGraphErrors();
      sprintf(name, "peakRate_pix_roc_%i", iRoc);
      peakRate_pixRoc[iRoc]->SetName(name);
      sprintf(title, "%s line, ROC %i", argv[1], iRoc);
      peakRate_pixRoc[iRoc]->SetTitle(title);
      peakRate_pixRoc[iRoc]->SetLineColor(iColor);
      peakRate_pixRoc[iRoc]->SetFillColor(kWhite);
      peakRate_pixRoc[iRoc]->SetMarkerColor(iColor);
      peakRate_pixRoc[iRoc]->SetMarkerStyle(iMarker);

      rate_pixRepRoc[iRoc] = new TGraphErrors();
      sprintf(name, "rate_pixRepRoc_%i", iRoc);
      rate_pixRepRoc[iRoc]->SetName(name);
      sprintf(title, "Pix rate, ROC %i", iRoc);
      rate_pixRepRoc[iRoc]->SetTitle(title);
      rate_pixRepRoc[iRoc]->SetLineColor(iColor);
      rate_pixRepRoc[iRoc]->SetFillColor(kWhite);
      rate_pixRepRoc[iRoc]->SetMarkerColor(iColor);
      rate_pixRepRoc[iRoc]->SetMarkerStyle(iMarker);
    }

  TMultiGraph* peakRep = new TMultiGraph();
  peakRep->SetName("peakRep");
  sprintf(title, "%s line vs repetition", argv[1]);
  peakRep->SetTitle(title);

  TMultiGraph* sigmaRep = new TMultiGraph();
  sigmaRep->SetName("sigmaRep");
  sprintf(title, "#sigma %s line vs repetition", argv[1]);
  sigmaRep->SetTitle(title);

  TGraphErrors* entriesRep = new TGraphErrors();
  entriesRep->SetName("entriesRep");
  entriesRep->SetTitle("Number of entries vs repetition");

  TMultiGraph* peakRate_pix = new TMultiGraph();
  peakRate_pix->SetName("peakRate_pix");
  sprintf(title, "%s line vs pix rate", argv[1]);
  peakRate_pix->SetTitle(title);

  TGraphErrors* rate_pixRep = new TGraphErrors();
  rate_pixRep->SetName("rate_pixRep");
  sprintf(title, "Rate %s line (pixels) vs repetition, whole module", argv[1]);
  rate_pixRep->SetTitle(title);

  // TH1D* peakDistr = new TH1D("peakDistr", "Distribution of the peaks;Charge [Vcal];Entries / bin", 3010, -0.5, 300.5);
  // TH1D* sigmaDistr = new TH1D("sigmaDistr", "Distribution of the peaks sigma;Charge [Vcal];Entries / bin", 360, -0.5, 35.5);

  int nPoint = 0;

  TDirectory* dir;
  TTree* tree;
  TF1* fit;
  TH1* hist;
  long int firedPix;
  long int firedPixRoc[nRoc];
  long int events;
  double ratePix;
  double ratePixErr;
  double* ion;
  double* ionErr;
  int ionPos;

  TDirectory* rocDir[nRoc];
  for(int iRoc = 0; iRoc < nRoc; ++iRoc)
    {
      sprintf(name, "ROC%i", iRoc);
      rocDir[iRoc] = outFile->mkdir(name);
    }

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

      totEntries = 0;

      for(int iRoc = 0; iRoc < nRoc; ++iRoc)
	{
	  sprintf(histName, "q_%s_C%i_V0", argv[1], iRoc);

	  hist = (TH1*) dir->Get(histName);
	  if(hist == 0)
	    {
	      std::cout << "\tWARNING: histogram " << histName << " not found in file " << inFile->GetPath() << std::endl;
	      continue;
	    }

	  totEntries += hist->GetEntries();

	  rocDir[iRoc]->cd();

	  fit = fitPeak(hist, negSigma, posSigma);
	  sprintf(name, "%s_%.02f", argv[1], rep);
	  hist->SetName(name);
	  hist->Write();

	  nPoint = peakRepRoc[iRoc]->GetN();
	  peakRepRoc[iRoc]->SetPoint(nPoint, rep, fit->GetParameter(1));
	  peakRepRoc[iRoc]->SetPointError(nPoint, 0, fit->GetParError(1));
	  
	  sigmaRepRoc[iRoc]->SetPoint(nPoint, rep, fit->GetParameter(2));
	  sigmaRepRoc[iRoc]->SetPointError(nPoint, 0, fit->GetParError(2));
	}

      nPoint = entriesRep->GetN();
      entriesRep->SetPoint(nPoint, rep, totEntries);
      entriesRep->SetPointError(nPoint, 0, sqrt(totEntries));

      // peakDistr->Fill(fit->GetParameter(1));
      // sigmaDistr->Fill(fit->GetParameter(2));

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

      for(int iRoc = 0; iRoc < nRoc; ++iRoc)
	firedPixRoc[iRoc] = 0;

      events = tree->GetEntries();
      for(long int i = 0; i < events; ++i)
	{
	  tree->GetEntry(i);
	  firedPix += npix;

	  for(int j = 0; j < npix; ++j)
	    ++firedPixRoc[proc[j]];
	}

      std::cout << "Tot entries in histos: " << totEntries << std::endl;
      std::cout << "Tot fired pix in tree: " << firedPix << std::endl;

      for(int iRoc = 0; iRoc < nRoc; ++iRoc)
	{
	  ratePix = firedPixRoc[iRoc] / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]
	  ratePixErr = sqrt(firedPixRoc[iRoc]) / (events * evtDuration * sensArea) * 1e-3; // [kHz cm^-2]

	  nPoint = rate_pixRepRoc[iRoc]->GetN();

	  rate_pixRepRoc[iRoc]->SetPoint(nPoint, rep, ratePix);
	  rate_pixRepRoc[iRoc]->SetPointError(nPoint, 0, ratePixErr);

	  ion = peakRepRoc[iRoc]->GetY();
	  ionErr = peakRepRoc[iRoc]->GetEY();
	  ionPos = peakRepRoc[iRoc]->GetN() - 1; // take the last point (should be right)

	  peakRate_pixRoc[iRoc]->SetPoint(nPoint, ratePix, ion[ionPos]);
	  peakRate_pixRoc[iRoc]->SetPointError(nPoint, ratePixErr, ionErr[ionPos]);
	}

      ratePix = firedPix / (events * evtDuration * sensArea * 16) * 1e-3; // [kHz cm^-2]
      ratePixErr = sqrt(firedPix) / (events * evtDuration * sensArea * 16) * 1e-3; // [kHz cm^-2]

      nPoint = rate_pixRep->GetN();

      rate_pixRep->SetPoint(nPoint, rep, ratePix);
      rate_pixRep->SetPointError(nPoint, 0, ratePixErr);
    }

  listStr.close();

  TCanvas* serv = new TCanvas("serv");
  for(int iRoc = 0; iRoc < nRoc; ++iRoc)
    {
      peakRepRoc[iRoc]->Draw("AP");
      peakRepRoc[iRoc]->GetXaxis()->SetTitle("Repetition");
      peakRepRoc[iRoc]->GetYaxis()->SetTitle("Measured ionization [Vcal]");

      sigmaRepRoc[iRoc]->Draw("AP");
      sigmaRepRoc[iRoc]->GetXaxis()->SetTitle("Repetition");
      sigmaRepRoc[iRoc]->GetYaxis()->SetTitle("Sigma [Vcal]");

      rate_pixRepRoc[iRoc]->Draw("AP");
      rate_pixRepRoc[iRoc]->GetXaxis()->SetTitle("Repetition");
      rate_pixRepRoc[iRoc]->GetYaxis()->SetTitle("Fired pixel rate [kHz cm^{-2}]");

      rocDir[iRoc]->cd();
      peakRepRoc[iRoc]->Write();
      sigmaRepRoc[iRoc]->Write();
      rate_pixRepRoc[iRoc]->Write();
    }
  delete serv;

  for(int iRoc = 0; iRoc < nRoc; ++iRoc)
    {
      peakRep->Add(peakRepRoc[iRoc]);
      sigmaRep->Add(sigmaRepRoc[iRoc]);
      peakRate_pix->Add(peakRate_pixRoc[iRoc]);
    }

  outFile->cd();

  TLegend* leg;

  TCanvas* peakRepCan = new TCanvas("peakRepCan");
  peakRepCan->SetGridx();
  peakRepCan->SetGridy();
  peakRep->Draw("APL");
  peakRep->GetXaxis()->SetTitle("Repetition");
  peakRep->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakRep->Write();
  leg = peakRepCan->BuildLegend();
  leg->SetFillColor(kWhite);
  peakRepCan->Modified();
  peakRepCan->Update();
  peakRepCan->Write();

  TCanvas* sigmaRepCan = new TCanvas("sigmaRepCan");
  sigmaRepCan->SetGridx();
  sigmaRepCan->SetGridy();
  sigmaRep->Draw("APL");
  sigmaRep->GetXaxis()->SetTitle("Repetition");
  sigmaRep->GetYaxis()->SetTitle("Sigma [Vcal]");
  sigmaRep->Write();
  leg = sigmaRepCan->BuildLegend();
  leg->SetFillColor(kWhite);
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
  peakRate_pix->Draw("AP");
  peakRate_pix->GetXaxis()->SetTitle("Rate fired pixels [kHz cm^{-2}]");
  peakRate_pix->GetYaxis()->SetTitle("Measured ionization [Vcal]");
  peakRate_pix->Write();
  leg = peakRatePixCan->BuildLegend();
  leg->SetFillColor(kWhite);
  peakRatePixCan->Modified();
  peakRatePixCan->Update();
  peakRatePixCan->Write();

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

  // peakDistr->Write();
  // sigmaDistr->Write();

  outFile->Close();

  return 0;
}
