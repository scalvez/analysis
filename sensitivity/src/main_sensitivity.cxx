#include <iostream>
#include <fstream>
#include <string>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include "channel_selection.h"
#include "pseudo_generator.h"
#include "multi_fit.h"
#include "sensitivity_measurements.h"
#include "analysis_config.h"
#include "sensitivity_constants.h"

int main(int argc, char* argv[]) {

  if(argc != 2) {
    std::cout << " ERROR : A seed must be provided" << std::endl;
    return 1;
  }

  double seed = atof(argv[1]);

  std::vector <double> se_2nu_measurements;
  std::vector <double> tl_measurements;
  std::vector <double> bi_measurements;
  std::vector <double> radon_measurements;

  // std::vector <double> tl_measurements = {2.05429,1.94539,1.73072,2.20376,2.06746,1.86915,2.25492,1.97681,1.84659,2.01918,2.09664,2.17444,1.78004,1.92657};
  // std::vector <double> bi_measurements = {9.89487,10.0564,10.3757,9.6727,9.87552,10.1698,9.59648,10.0101,10.2031,9.94669,9.83203,9.71597,10.3024,10.0845};

  std::cout << " Starting the sensitivity measurement..." << std::endl;

  std::cout << " Channel selection..." << std::endl;
  if (generate_pdf) {
    for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
      const std::string & a_isotope = i->first;
      std::cout << "      Isotope : " << a_isotope << std::endl;
      channel_selection(a_isotope, quantities);
    }
  }

  std::cout << " Channel selection operated" << std::endl;

  if (generate_pseudo) {

    gSystem->Exec("rm ../pseudo/*_pseudo.root");

    std::cout << " Pseudo experiments generation " << std::endl;
    // unsigned int number_of_pseudo_experiments = 100;
    for(unsigned int n_pseudo = 0; n_pseudo < number_of_pseudo_experiments; ++n_pseudo) {
      std::cout << " [] Pseudo experiment n°" << n_pseudo+1 << " / "
                << number_of_pseudo_experiments << std::endl;
      std::map < std::string, std::vector<double> > activity_measurement;
      for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
        const std::string & a_isotope = i->first;
        const double & a_activity = i->second;
        pseudo_generator(a_isotope,quantities,a_activity,seed);
      }
      gSystem->Exec("hadd -f ../pseudo/pseudo.root ../pseudo/*_pseudo.root");
      if(fit) {
        multi_fit(activity_measurement);
        std::cout << " - Reconstructed activities - " << std::endl;
        for (auto i = activity_measurement.begin(); i != activity_measurement.end(); ++i) {
          std::string isotope = i->first;
          double activity = i->second.at(0);
          std::cout << " Isotope : " << isotope << std::endl;
          std::cout << " Activity : " << activity * 1e6  << " uBq/kg " << std::endl;

          if(isotope == "2nu")
            se_2nu_measurements.push_back(activity);
          if(isotope == "tl208")
            tl_measurements.push_back(activity * 1e6);
          // h_tl_meas->Fill(activity * 1e6);
          if(isotope == "bi214")
            bi_measurements.push_back(activity * 1e6);
          // h_bi_meas->Fill(activity * 1e6);
          if(isotope == "radon")
            radon_measurements.push_back(activity * 1e6);
        }

        // TH1F *h_2nu_meas = new TH1F("h_2nu_meas","h_2nu_meas",100,7,11);
        // TH1F *h_tl_meas = new TH1F("h_tl_meas","h_tl_meas",100,0,4);
        // TH1F *h_bi_meas = new TH1F("h_bi_meas","h_bi_meas",100,8,12);
        // TH1F *h_radon_meas = new TH1F("h_radon_meas","h_radon_meas",100,140,160);
        // TFile *f_output= new TFile("../measurements.root", "RECREATE");

        // for(auto i = 0; i < tl_measurements.size(); ++i) {
        //   h_2nu_meas->Fill(se_2nu_measurements.at(i));
        //   h_tl_meas->Fill(tl_measurements.at(i));
        //   h_bi_meas->Fill(bi_measurements.at(i));
        //   h_radon_meas->Fill(radon_measurements.at(i));
        // }
        // h_2nu_meas->Write();
        // h_tl_meas->Write();
        // h_bi_meas->Write();
        // h_radon_meas->Write();
      }
    }

    ofstream myfile;
    myfile.open ("measurements.txt", std::ios::app);
    for(unsigned int i = 0; i < se_2nu_measurements.size(); ++i) {
      myfile << "2nu " << se_2nu_measurements.at(i) << "  "
             << "tl " << tl_measurements.at(i) << "  "
             << "bi " << bi_measurements.at(i) << "  "
             << "radon " << radon_measurements.at(i) << "\n";
    }
  }

  quantity_efficiency.clear();
  quantity_pdf.clear();

  return 0;
}

// void main_sensitivity() {
//   main();
//   return;
// }
