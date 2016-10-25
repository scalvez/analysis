#include <iostream>
#include <fstream>
#include <string>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TSystem.h>
#include "TRandom.h"

#include "stats_tools.h"
#include "analysis_config.h"
#include "pseudo_generator.h"

int main(int argc, char* argv[]) {

  // if(argc != 2) {
  //   std::cout << " ERROR : A seed must be provided" << std::endl;
  //   return 1;
  // }

  // double seed = atof(argv[1]);

  // double test = atof(argv[1]);

  // Priors on H and Hbar
  double p0_H = 0.5;
  double p0_Hbar = 0.5;

  TH1F* h_s90 = new TH1F("h","h",10,0,10);
  TGraph* g_p = new TGraph();

  // Background and signal PDFs
  std::vector<double> pdf_B;
  std::vector<double> pdf_S;

  //Generate PDFs
  for(unsigned int i=0; i<100; ++i) {
    pdf_B.push_back(0.01);
    pdf_S.push_back( TMath::Gaus(i,48,1,true));
  }

  TFile *f= new TFile("output.root", "RECREATE");

  TRandom *rdm = new TRandom();
  rdm->SetSeed(1);

  double p = 0;

  std::vector<double> s_90;

  TH1 *h_data = new TH1I("data","data",100,0,100);

  // for(unsigned int i_pseudo =0; i_pseudo<200; ++i_pseudo) {
  for(unsigned int i_pseudo =0; i_pseudo<1; ++i_pseudo) {
    std::cout << "pseudo_experiment n°" << i_pseudo << std::endl;

    //rescale signal if pseudo signal

    //Generate data from PDFs
    std::vector<double> data;

    for(unsigned int i=0; i<100; ++i)
      data.push_back(0);

    // int n_bkg_random =8;
    int n_bkg_random =rdm->Poisson(10);
    std::cout << "n_bkg_random " << n_bkg_random << std::endl;

    for(unsigned int j=0; j<n_bkg_random; ++j) {
      data[int(rdm->Uniform()*100)]++;
    }

    // int n_sig_random =rdm->Poisson(20);
    int n_sig_random = 0;
    std::cout << "n_sig_random " << n_sig_random << std::endl;

    for(unsigned int j=0; j<n_sig_random; ++j) {
      int i_sig = rdm->Gaus(48,1);
      data[i_sig]++;
    }

    for(auto i_data = 1; i_data != data.size(); ++i_data) {
      h_data->Fill(i_data,data[i_data]);
    }

    // int test_s = 10;
    // int test_b = 108;

    // std::cout << " test proba " << proba_bayes(test_s,test_b,pdf_S,pdf_B,data)*p0_S(test_s)*p0_B(test_b) << std::endl;

    // return 0;

    double tot_prob = 0;
    //Loop over S90
    for(double s = 0; s<50; ++s) {
    // for(double s = 0; s<=0; ++s) {
      std::cout << " Computing proba for signal hypothesis S = " << s << std::endl;
      p = 0;
      // for(double b=0; b<=1000; b += 10) {
        for(double b=0; b<=100; ++b) {
        // std::cout << "    ... Integrating proba for background B=" << b << std::endl;

        double num = proba_bayes(s,b,pdf_S,pdf_B,data)*p0_S(s)*p0_B(b);
        // std::cout << " num proba bayes  s="<< s << "  b=" << b << "    " << num << std::endl;
        double denum = 0;
        for(double i_s = 0; i_s<=50; ++i_s)
          // for(double i_b = 0; i_b<=1000; i_b += 10) {
          for(double i_b = 0; i_b<=100; ++i_b) {

              denum += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
          }
        // std::cout << " denum proba bayes" << denum << std::endl;
        p+=num/denum;
      }
      g_p->SetPoint(s,s,p);
      std::cout << "signal events : " << s << "   ;  proba = " << p << std::endl;
      tot_prob += p;

      //For limit setting when signal is zero
      if(tot_prob>=0.9) {
        std::cout << "---" << "s_90 =" << s << std::endl;
        s_90.push_back(s);
        h_s90->Fill(s);
        break;
      }

    }

  } //end pseudo experiments loop

  h_data->Write();
  h_s90->Write();
  g_p->Write();
  f->Close();
  // std::cout << "proba is " << p << std::endl;
  return 0;
}

/** Piece of code 1 **/
    // Generate data from PDFs

    // for(unsigned int i=0; i<100; ++i)
    //   data.push_back(0);

    // int n_bkg_random =rdm->Poisson(8);
    // // int n_bkg_random =8;
    // std::cout << "n_bkg_random " << n_bkg_random << std::endl;

    // for(unsigned int j=0; j<n_bkg_random; ++j) {
    //   data[int(rdm->Uniform()*100)]++;
    // }

    // // int n_sig_random =rdm->Poisson(20);
    // int n_sig_random =20;
    // std::cout << "n_sig_random " << n_sig_random << std::endl;

    // for(unsigned int j=0; j<n_sig_random; ++j) {
    //   int i_sig = rdm->Gaus(48,1);
    //   data[i_sig]++;
    //   std::cout << " signal random gauss " << i_sig << std::endl;
    // }
/** End of piece of code 1 **/
