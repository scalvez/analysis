#include <iostream>
#include <fstream>
#include <string>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include "TRandom.h"

#include "stats_tools.h"

int main(int argc, char* argv[]) {

  // if(argc != 2) {
  //   std::cout << " ERROR : A seed must be provided" << std::endl;
  //   return 1;
  // }

  // double seed = atof(argv[1]);

  // double test = atof(argv[1]);

  double p0_H = 0.5;
  double p0_Hbar = 0.5;

  // double p0_S = 0.5;
  // double p0_B = 0.5;

  // std::vector<double> pdf_B = {0.2,0.2,0.2,0.2,0.2};
  // std::vector<double> pdf_S = {0,0,1,0,0};
  // std::vector<double> data = {1,1,1,1,1};


  // double B = 9;
  // double S = 0;

  // std::cout << " test " << proba_bayes(S,B,pdf_S,pdf_B,pdf_B) << std::endl;

  // double tmp =0;
  //   for(double i_s = 0; i_s<=100; ++i_s)
  //     for(double i_b = 0; i_b<=100; ++i_b) {
  //       tmp += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
  //     }

    // std::cout << " integral is " << tmp << std::endl;
  std::vector<double> pdf_B;
  std::vector<double> pdf_S;

  TFile *f= new TFile("output.root", "RECREATE");
  TH1F* h_s90 = new TH1F("h","h",10,0,10);
  for(unsigned int i=0; i<100; ++i) {
    pdf_B.push_back(0.01);
    pdf_S.push_back( TMath::Gaus(i,48,2,true));
  }

  TRandom *rdm = new TRandom();
  rdm->SetSeed(2);

  double p = 0;
  // double s = 0;
  double tot_prob = 0;

  std::vector<double> s_90;
  for(unsigned int i_pseudo =0; i_pseudo <1; ++i_pseudo) {
    std::cout << "pseudo_experiment n°" << i_pseudo << std::endl;

  std::vector<double> data;

  for(unsigned int i=0; i<100; ++i)
    data.push_back(0);

  int n_random =rdm->Poisson(10);
  std::cout << "n_random " << n_random << std::endl;

  for(unsigned int j=0; j<n_random; ++j) {
    std::cout << "uniform " << int(rdm->Uniform()*100) << std::endl;
    data[int(rdm->Uniform()*100)]++;
  }

  for(double s = 0; s<10; ++s) {
      for(double b=0; b<=100; ++b) {
        double num = proba_bayes(s,b,pdf_S,pdf_B,data)*p0_S(s)*p0_B(b);
        double denum = 0;
        for(double i_s = 0; i_s<=100; ++i_s)
          for(double i_b = 0; i_b<=100; ++i_b) {
            denum += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
          }
        p+=num/denum;
      }
      tot_prob += p;
      if(p>=0.9) {
        std::cout << "---" << "s_90 =" << s << std::endl;
        s_90.push_back(s);
        h_s90->Fill(s);
        break;
      }
    }
  }
  h_s90->Write();
  f->Close();
  // std::cout << "proba is " << p << std::endl;
  return 0;
}
