#ifndef PSEUDO_GENERATOR_H
#define PSEUDO_GENERATOR_H 1

#include <TString.h>
#include <vector>
#include <map>

void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity, double seed);
// void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity);

#endif
