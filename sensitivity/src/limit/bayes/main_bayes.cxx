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

  // //Generate PDFs
  // for(unsigned int i=0; i<100; ++i) {
  //   pdf_B.push_back(0.01);
  //   pdf_S.push_back( TMath::Gaus(i,48,1,true));
  // }

  TFile *f_sig = TFile::Open("../0nu_pdf.root");
  TH1 *h_pdf_sig = (TH1*)f_sig->Get("2e_electrons_energy_sum");

  TFile *f_bkg = TFile::Open("../2nu_pdf.root");
  TH1 *h_pdf_bkg = (TH1*)f_bkg->Get("2e_electrons_energy_sum");

  //efficiency = subrange of pdf * eff MC * MC>2MeV correction
  double eff_2nu =  h_pdf_bkg->Integral(41,80) * h_pdf_bkg->GetEntries()/1e6 / 25;
  int n_2nu_events = eff_2nu * mass * exposure_y * Na * log(2) / M_Se / halflife_2nu;

  std::cout << "Expected number of background events " << n_2nu_events << std::endl;

  TFile *f= new TFile("output.root", "RECREATE");

  TH1F *h_bkg_pdf_trunc = new TH1F("bkg_pdf_trunc","bkg_pdf_trunc",40,2,4);
  TH1F *h_sig_pdf_trunc = new TH1F("sig_pdf_trunc","sig_pdf_trunc",40,2,4);

  for(auto i_bkg_pdf = 41; i_bkg_pdf <= 80; ++i_bkg_pdf) {
    h_bkg_pdf_trunc->Fill(h_pdf_bkg->GetXaxis()->GetBinCenter(i_bkg_pdf),
                          h_pdf_bkg->GetBinContent(i_bkg_pdf));
    h_sig_pdf_trunc->Fill(h_pdf_sig->GetXaxis()->GetBinCenter(i_bkg_pdf),
                          h_pdf_sig->GetBinContent(i_bkg_pdf));
  }

  // Background and signal PDFs
  std::vector<double> pdf_B;
  std::vector<double> pdf_S;

  double norm_factor_2nu = h_pdf_bkg->GetEntries()/1e6 / 25;
  double norm_factor_0nu = h_pdf_sig->GetEntries()/1e6;

  for(auto i_pdf = 1; i_pdf <=h_bkg_pdf_trunc->GetNbinsX();++i_pdf) {
    pdf_B.push_back(norm_factor_2nu*h_bkg_pdf_trunc->GetBinContent(i_pdf));
    pdf_S.push_back(norm_factor_0nu*h_sig_pdf_trunc->GetBinContent(i_pdf));
    std::cout << " i_pdf " << i_pdf << "  " << pdf_B.back() << "  " << pdf_S.back() << std::endl;
  }

  h_bkg_pdf_trunc->Scale(1./h_bkg_pdf_trunc->Integral(1,40));

  h_pdf_sig->SetAxisRange(2,4);
  h_pdf_bkg->SetAxisRange(2,4);

  TRandom *rdm = new TRandom();
  rdm->SetSeed(1);

  double p = 0;

  std::vector<double> s_90;

  TH1 *h_data = new TH1I("data","data",40,2,4);

  // for(unsigned int i_pseudo =0; i_pseudo<200; ++i_pseudo) {
  for(unsigned int i_pseudo =0; i_pseudo<1; ++i_pseudo) {
    std::cout << "pseudo_experiment n°" << i_pseudo << std::endl;

    h_data->Reset();

    //rescale signal if pseudo signal

    pseudo_generator(h_bkg_pdf_trunc,h_data,n_2nu_events,1);

    std::vector<double> data;
    for(auto i_h = 1; i_h <=h_data->GetNbinsX();++i_h)
      data.push_back(h_data->GetBinContent(i_h));

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

  double tot_prob = 0;
  for(double s = 0; s<50; ++s) {
    p = 0;
    for(double b=0; b<=10000; b += 100) {
      double num = proba_bayes(s,b,pdf_S,pdf_B,data)*p0_S(s)*p0_B(b);
      double denum = 0;
      for(double i_s = 0; i_s<=100; ++i_s)
        for(double i_b = 0; i_b<=10000; i_b += 100) {
          denum += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
        }
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
  }

  h_pdf_bkg->Write();
  h_bkg_pdf_trunc->Write();
  // h_data->Write();

  h_s90->Write();
  g_p->Write();
  f->Close();
  // std::cout << "proba is " << p << std::endl;
  return 0;
}
