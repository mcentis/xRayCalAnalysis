#ifndef CALSINGLET_HH
#define CALSINGLET_HH

#include "vector"

class TF1;
class TH1;
class TGraphErrors;
class TFile;
class TDirectory;

class calSingleT
{
public:
  calSingleT(const char* inFilePath, float temperature);
  ~calSingleT();

  void DoCal();

  double GetTemperature() {return temp;};
  TF1* GetCalFit() {return calFit;};
  TGraphErrors* GetCalGraph() {return calGr;};
  std::vector<TH1*> GetHistos() {return histos;};
private:
  TGraphErrors* calGr;
  TF1* calFit;
  TDirectory* dir;
  TFile* inFile;
  float temp;
  std::vector<TH1*> histos;
};

#endif // #ifndef CALSINGLET_HH
