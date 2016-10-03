#include "analysis_config.h"

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax) {

  // hardcoded for now, see if it is parametrized in the configuration file
  if(quantity.Contains("energy")) {
    nbins = 100;
    xmin = 0;
    xmax = 5;
  }
  else if(quantity.Contains("probability")){
    nbins = 100;
    xmin = 0;
    xmax = 1;
  }
  else if(quantity.Contains("cos_angle")) {
    nbins = 100;
    xmin = -1;
    xmax = 1;
  }
  else if(quantity.Contains("angle")) {
    nbins = 100;
    xmin = 0;
    xmax = 3.14159;
  }
  else if (quantity.Contains("track_length")) {
    nbins = 100;
    xmin = 0;
    xmax = 2000;
  }
  else {
    //Also maybe the alpha delayed time
    //for now, dummy values
    nbins = 100;
    xmin = 0;
    xmax = 1000;
  }

  return;
}

// TODO: to move in the analysis config and get MC dataset size from conf
int get_isotope_mc_size(TString isotope) {
  if(isotope.Contains("2nu"))
    return 1e8;
  else if(isotope.Contains("tl208"))
    return 1e8;
  else if(isotope.Contains("bi214"))
    return 1e8;
  else if(isotope.Contains("radon"))
    return 1e8;
  else
    return 1;
    //tmp dirty
}
