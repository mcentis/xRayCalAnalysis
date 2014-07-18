/*
 * program that performs one calibration
 */

#include "calSingleT.hh"

#include "iostream"
#include "vector"
#include "stdlib.h"

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TSystem.h"

int main(int argc, char* argv[])
{
  if(argc != 4)
    {
      std::cout << "\tusage: calOneT inFile temp outFile" << std::endl;
      return 1;
    }

  gSystem->Load("libTree"); // to avoid warnings

  calSingleT* cal = new calSingleT(argv[1], atof(argv[2]));
  cal->DoCal();

  TFile* outFile = new TFile(argv[3], "RECREATE");
  std::vector<TH1*> hiVec = cal->GetHistos();

  for(unsigned int i = 0; i < hiVec.size(); ++i)
    hiVec.at(i)->Write();

  cal->GetCalGraph()->Write();
      
  outFile->Close();

  return 0;
}
