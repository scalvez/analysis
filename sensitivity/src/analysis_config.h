#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H 1

#include "TCut.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>

const bool generate_pdf = true;
const bool generate_pseudo = true;
const bool poisson_pseudo = false;
const unsigned int number_of_pseudo_experiments = 1;
const bool fit = true;
const bool print_fits = true;

// const double bi214_channel_1e1a_efficiency = 10733./2500000;

// const double tl208_channel_1e1g_efficiency = 0.139637;
// const double bi214_channel_1e1g_efficiency = 0.123019;

// const TString internal_probability_min = "0.04";
// const double internal_probability_max = 1;
// const double vertices_probability_min = 0;
// const double vertices_probability_max = 0.01;

// // 2e_int channel beta-beta-like events
// const std::string good_internal_probability_min = "2e_electrons_internal_probability > 0.04 ";

// // const char prob_cut[200] = good_internal_probability_min;
// const std::string test_string = "2e_electrons_internal_probability > 0.04 ";
// const char *prob_cut = test_string.c_str();

// // const TString prob_cut = "2e_electrons_internal_probability > " + internal_probability_min;

// const TCut good_internal_probability_cut = prob_cut;

const std::map < std::string, double > isotope_activity = {
  {"2nu",9},
  {"bi214",10e-6},
  {"radon",150e-6},
  {"tl208",2e-6}
};

// const std::map < std::string, double > isotope_activity = {
//   {"2nu",9},
//   {"bi214",0},
//   {"radon",0},
//   {"tl208",0}
// };

// //radon
// const std::vector <TString> quantities = {
//   //not understood "1e1a_alpha_angle",
//   "1e1a_alpha_track_length",
//   "1e1a_electron_alpha_angle",
//   "1e1a_electron_alpha_vertices_probability",
// };

// // 2nu
// const std::vector <TString> quantities = {
//   "1e_electron_energy",
//   "2e_electron_maximal_energy",
//   "2e_electrons_energy_sum",
//   "2e_electrons_internal_probability",
//   "2e_electrons_cos_angle",
//   //maybe "2e_electrons_vertices_probability",
// };

// tl and bi
const std::vector <TString> quantities = {
  // "1e_electron_energy",
  "1e1g_gamma_energy",
  "1e1g_electron_gamma_energy_sum",
  "1e2g_gamma_max_energy",
  "1e2g_electron_gammas_energy_sum",
  "1e3g_electron_gammas_energy_sum",
  // "2e1g_electrons_gammas_energy_sum"
};

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax);

int get_isotope_mc_size(TString isotope);

#endif
