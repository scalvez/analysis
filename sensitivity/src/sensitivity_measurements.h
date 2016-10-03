#ifndef SENSITIVITY_MEASUREMENTS_H
#define SENSITIVITY_MEASUREMENTS_H 1

#include <map>
#include <TH1.h>

std::map < TString , double > quantity_efficiency;
std::map < TString , TH1F* > quantity_pdf;

#endif
