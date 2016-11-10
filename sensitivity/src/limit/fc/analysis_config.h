#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H 1

const int n_mc_0nu = 1e6;
const int n_mc_2nu = 1e6;

const double mass = 7.;
const double exposure_sec = 2.5 * 3.14e7;
const double exposure_y = 2.5;
// const double exposure_sec = 1./12 * 3.14e7;
// const double exposure_y = 1./12;
const double tracker_volume = 15.3; //m3
const double halflife_2nu = 9e19; // years;
const double Na = 6.022e23;
const double M_Se = 0.082; //kg/mol more like 81.6 if we enrich at 90%
const double const_se = log(2) * Na / M_Se; // years;

const double halflife_0nu = 1e25; // years;

// N = m * t * epsilon * Na * log(2) / M / T_1/2      //ignore enrichment

#endif
