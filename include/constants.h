/*
 * file with definition of constants for the programs
 */
#ifndef CONSTANTS_H
#define CONSTANTS_H

const double CuIon = (8027.84 + 8047.82) / (2 * 3.6); // eh pairs produced (K_alpha1 + K_alpha2) / (2 * 3.6)
const double MoIon = (17374.29 + 17479.372) / (2 * 3.6);
const double AgIon = (21990.30 + 22162.917) / (2 * 3.6);
const double SnIon = (25044.04 + 25271.36) / (2 * 3.6);
const double TeIon = (27201.99 + 27472.57) / (2 * 3.6);
const double BaIon = (31816.615 + 32193.262) / (2 * 3.6);

// roc dimensions
const int nCol = 52; // number of columns in a roc
const int nRow = 80; // number of rows in a roc

const int nRoc = 16;

#endif //#ifndef CONSTANTS_H
