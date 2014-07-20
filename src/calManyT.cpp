/*
 * program to perform the temperature calibration at many temperatures
 */

#include "calSingleT.hh"

#include "TFile.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TSystem.h"

#include "iostream"
#include "fstream"
#include "vector"

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "\tusage: calManyT listFile dataDir" << std::endl;
      return 1;
    }

  gSystem->Load("libTree"); // to avoid warnings

  std::vector<calSingleT*> singleTs;
  calSingleT* single;
  char fileName[200];

  float temp;
  std::string file;

  std::ifstream listStr(argv[1], std::ifstream::in);
  if(listStr.is_open() == false)
    {
      std::cout << "\tERROR: could not open " << argv[1] << std::endl;
      return 1;
    }

  std::cout << "\tFiles being read:\n";

  while(listStr.good())
    {
      listStr >> temp >> file;
      std::cout << '\t' << file << '\t' << temp << std::endl;
      sprintf(fileName, "%s/%s", argv[2], file.c_str());
      single = new calSingleT(fileName, temp);
      single->DoCal();
      singleTs.push_back(single);
    }

  std::cout << "\tIMPORTANT: if a file appears many times, correct the list file for newlines, multiple entries with same temperature cause crashes" << std::endl;

  listStr.close();

  sprintf(fileName, "%s/calManyT.root", argv[2]);
  TFile* outFile = new TFile(fileName, "RECREATE");

  TDirectory* dir;
  std::vector<TH1*> hiVec;
  char name[50];
  for(unsigned int iTemp = 0; iTemp < singleTs.size(); ++iTemp) // loop on the single temperatures and write the hitos and graphs
    {
      single = singleTs.at(iTemp);
      sprintf(name, "%.00fC", single->GetTemperature());
      dir = outFile->mkdir(name);
      dir->cd();
      hiVec = single->GetHistos();

      for(unsigned int iHi = 0; iHi < hiVec.size(); ++iHi) // write the histos of this temperature
	hiVec.at(iHi)->Write();

      single->GetCalGraph()->Write();
    }

  TMultiGraph* allTcal = new TMultiGraph("allTcal", "Calibration curves taken at differen temperatures");
  TGraphErrors* gr;
  TF1* fit;
  char title[50];
  int iColor;

  for(unsigned int iTemp = 0; iTemp < singleTs.size(); ++iTemp) // loop on the single temperatures and prepare merged plots
    {
      single = singleTs.at(iTemp);
      gr = single->GetCalGraph();
      sprintf(title, "%.00fC", single->GetTemperature());
      gr->SetTitle(title);
      gr->SetFillColor(kWhite);
      iColor = iTemp % 9 + 1;
      if(iColor == 5) ++iColor; // skip yellow
      gr->SetLineColor(iColor); // set line color and style
      gr->SetMarkerColor(iColor);
      fit = gr->GetFunction("calFit");
      fit->SetLineColor(iColor);
      allTcal->Add(gr);
    }

  outFile->cd();

  TCanvas* allCalCan = new TCanvas("allCalCan");
  allCalCan->SetGridx();
  allCalCan->SetGridy();
  allTcal->Draw("APL");
  allTcal->GetXaxis()->SetTitle("Measured ionization [Vcal]");
  allTcal->GetYaxis()->SetTitle("Teoretical charge deposit [e]");
  allTcal->Write();
  TLegend* leg = allCalCan->BuildLegend();
  leg->SetFillColor(kWhite);
  allCalCan->Modified();
  allCalCan->Update();
  allCalCan->Write();

  outFile->Close();

  return 0;
}
