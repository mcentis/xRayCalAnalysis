#ifndef CLUSTERDEFANDFUNCTION_HH
#define CLUSTERDEFANDFUNCTION_HH

struct cluster // define cluster structure
{
  std::vector<UChar_t> col;
  std::vector<UChar_t> row;
  std::vector<Double_t> charge;
  int nPix;
  double totCharge;

  cluster()
  {
    nPix = 0;
    totCharge = 0;
  };
};

std::vector<cluster> clusterEvent(UShort_t nPix, UChar_t* col, UChar_t* row, Double_t* charge, double distCut) // function to cluster an event, the distCut is the clustering distance
{
  std::vector<cluster> clustVec;

  if(nPix == 0) return clustVec; // no pixels in this event

  std::vector<bool> isUsed(nPix, false); // vector with info if pixel of event is already used (clustered)
  int totUsed = 0; // number of used pixels
  int firstNotUsed = 0; // number of used pixels

  cluster* clu;
  bool growing; // cluster is growing
  double distance;

  while(totUsed < nPix)
    {
      firstNotUsed = -1;
      for(int i = 0; i < nPix; i++) // find first not used pixel
	{
	  if(isUsed[i] == false)
	    {
	      firstNotUsed = i;
	      break;
	    }
	}

      if(firstNotUsed == -1)
	{
	  std::cout << "\tWARNING: Unexpected end of clustering" << std::cout;
	  break;
	}

      clu = new cluster(); // start cluster

      clu->col.push_back(col[firstNotUsed]); // add pixel to cluster
      clu->row.push_back(row[firstNotUsed]);
      clu->charge.push_back(charge[firstNotUsed]);
      clu->nPix++;
      clu->totCharge += charge[firstNotUsed];
      isUsed[firstNotUsed] = true;
      totUsed++;

      do
	{
	  growing = false;
	  
	  for(int i = 0; i < nPix; ++i)
	    {
	      if(isUsed[i] == false) // not used pixels
		{
		  for(int j = 0; j < clu->nPix; ++j)
		    {
		      distance = sqrt(pow(clu->col[j] - col[i], 2) + pow(clu->row[j] - row[i], 2));
		      if(distance <= distCut)
			{
			  clu->col.push_back(col[i]); // add pixel to cluster
			  clu->row.push_back(row[i]);
			  clu->charge.push_back(charge[i]);
			  clu->nPix++;
			  clu->totCharge += charge[i];
			  isUsed[i] = true;
			  totUsed++;
			  growing = true;
			  break; // 1 pix at time, and preserve growing information
			}
		    }
		}
	    }
	} while(growing);

      clustVec.push_back(*clu);
      delete clu;
    }

  return clustVec;
}

#endif
