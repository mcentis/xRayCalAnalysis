#ifndef ROCTOMODULECR_HH
#define ROCTOMODULECR_HH

// the constants should be moved elsewere !!!!
// roc dimensions
const int nCol = 52; // number of columns in a roc
const int nRow = 80; // number of rows in a roc

void moduleColRow(int roc, int col, int row, int* modCR); // calculates module col and row

#endif //#ifndef ROCTOMODULECR_HH
