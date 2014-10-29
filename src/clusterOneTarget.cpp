/*
 * program to cluster the data from the tree of the run of one target
 */

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
//#include "TSystem.h"

#include "iostream"
#include "math.h"
#include "string"

struct pixel
{
  int roc;
  int col;
  int row;
  int ph;
  double q;
};

struct cluster
{
  std::vector<pixel> pixels;
  double qtot;
};

int main(int argc, char* argv[])
{
  if(argc != 3)
    {
      std::cout << "\tusage: clusterOneTarget element inputFile" << std::endl;
      return 1;
    }

  //gSystem->Load("libTree"); // to avoid warnings

  TFile* inFile = TFile::Open(argv[2]);
  if(!inFile)
    {
      std::cout << "Can not open " << argv[2] << std::endl;
      return 1;
    }
  TDirectory * dir = (TDirectory*) inFile->Get("Xray");
  if(!dir)
    {
      std::cout << "Can not find TDirectory Xray" << std::endl;
      return 1;
    }
  TTree* tree = (TTree*) dir->Get("events");
  if(!tree)
    {
      std::cout << "Can not find TTree events" << std::endl;
      return 1;
    }

  //Declaration of leaves types
  UShort_t        header;
  UShort_t        trailer;
  UShort_t        npix;
  UChar_t         proc[2];
  UChar_t         pcol[2];
  UChar_t         prow[2];
  Double_t        pval[2];
  Double_t        pq[2];

  // Set branch addresses.
  tree->SetBranchAddress("header",&header);
  tree->SetBranchAddress("trailer",&trailer);
  tree->SetBranchAddress("npix",&npix);
  tree->SetBranchAddress("proc",proc);
  tree->SetBranchAddress("pcol",pcol);
  tree->SetBranchAddress("prow",prow);
  tree->SetBranchAddress("pval",pval);
  tree->SetBranchAddress("pq",pq);

  std::string outFileName(argv[2]);
  outFileName.insert(outFileName.size() - 5, "_clustered");
  TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

  char histName[50]; // name of the histo to look for
  sprintf(histName, "q_%s_C0_V0", argv[1]);
  TH1* nonClusQ = (TH1*) dir->Get(histName); // get the non clustered hist
  TAxis* xAx = nonClusQ->GetXaxis();
  TH1D* clusQ = new TH1D(histName, "Clustered charge distribution", xAx->GetNbins(), xAx->GetXmin(), xAx->GetXmax());

  sprintf(histName, "q_%s_C0_V0_notClustered", argv[1]);
  nonClusQ->SetName(histName);

  outFile->cd();
  nonClusQ->Write(); // write the old distr for comparison

  // std::vector<pixel> hitVec;
  // std::vector<cluster> cluVec;
  // pixel hit;
  // cluster clu;
  double dist;
  double distCut = 1; // distance cut for the clustering

  for(long int i = 0; i < tree->GetEntries(); ++i)
    {
      tree->GetEntry(i);

      // for(int n = 0; n < npix; ++n) //store all the hits
      // 	{
      // 	  hit.roc = proc[n];
      // 	  hit.col = pcol[n];
      // 	  hit.row = prow[n];
      // 	  hit.ph = pval[n];
      // 	  hit.q = pq[n];
      // 	  hitVec.push_back(hit);
      // 	}

      // for(std::iterator::vector<pixel> it; it == hitVec.end(); it++)

      // megadirty solution!!! (the array is at most of size 2)
      if(npix >= 1)
	{  
	  if(npix == 2)
	    {
	      dist = sqrt(pow(pcol[1] - pcol[0], 2) + pow(prow[1] - prow[0], 2));
	      if(dist <= distCut)
		clusQ->Fill(pq[0] + pq[1]);
	      else
		{
		  clusQ->Fill(pq[0]);
		  clusQ->Fill(pq[1]);
		}
	    }
	  else
	    clusQ->Fill(pq[0]);
	}
    }

  clusQ->Write();

  outFile->Close();

  return 0;
}
