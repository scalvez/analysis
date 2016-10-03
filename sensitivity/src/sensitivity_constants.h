#ifndef SENSITIVITY_CONSTANTS_H
#define SENSITIVITY_CONSTANTS_H 1

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

#endif
