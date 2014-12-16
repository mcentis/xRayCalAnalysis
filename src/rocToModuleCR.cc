#include "rocToModuleCR.hh"

void moduleColRow(int roc, int col, int row, int* modCR) // calculates module col and row
{
  modCR[0] = -1;
  modCR[1] = -1;

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
