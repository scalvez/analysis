#ifndef CHANNEL_SELECTION_H
#define CHANNEL_SELECTION_H 1

#include <TString.h>
#include <vector>
#include <map>

void channel_selection(TString isotope, std::vector<TString> quantities_pdf, bool normalize = true);
// void channel_selection(TString isotope, std::vector<TString> quantities_pdf, std::map < TString , double > & quantity_efficiency, bool normalize = true);

#endif
