/*
 * file with definition of constants for the program
 */
#ifndef CONSTANTS_H
#define CONSTANTS_H

// FIXME correct copper energy!!!
const double CuIon = (8000) / (2 * 3.6); // eh pairs produced (K_alpha1 + K_alpha2) / (2 * 3.6)
const double MoIon = (17374.29 + 17479.372) / (2 * 3.6);
const double AgIon = (21990.30 + 22162.917) / (2 * 3.6);
const double TeIon = (27201.99 + 27472.57) / (2 * 3.6);
const double BaIon = (31816.615 + 32193.262) / (2 * 3.6);

const int nHist = 3; // number of histos
const char* histNames[nHist] = {"q_Cu_C0_V0", "q_Mo_C0_V0", "q_Ag_C0_V0"}; // histo names
const double peakIon[nHist] = {CuIon, MoIon, AgIon}; // ionization assigned to the peacks

#endif //#ifndef CONSTANTS_H
