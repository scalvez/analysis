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

  // Get the signal pdf
  TFile *f_sig = TFile::Open("../0nu_pdf.root");
  TH1 *h_pdf_sig = (TH1*)f_sig->Get("2e_electrons_energy_sum");

  // Get the background pdf
  TFile *f_bkg = TFile::Open("../2nu_pdf.root");
  TH1 *h_pdf_bkg = (TH1*)f_bkg->Get("2e_electrons_energy_sum");

  //efficiency = subrange of pdf [2;4]MeV * eff MC * MC>2MeV correction
  double eff_2nu =  h_pdf_bkg->Integral(49,80) * h_pdf_bkg->GetEntries()/1e6 / 25;
  // int n_2nu_events = eff_2nu * mass * exposure_y * Na * log(2) / M_Se / halflife_2nu;
  int n_2nu_events = 50;

  int n_0nu_events = 0;

  std::cout << "Expected number of background events " << n_2nu_events << std::endl;

  TFile *f= new TFile("output.root", "RECREATE");

  TH1F *h_bkg_pdf_trunc = new TH1F("bkg_pdf_trunc","bkg_pdf_trunc",32,2.4,4);
  TH1F *h_sig_pdf_trunc = new TH1F("sig_pdf_trunc","sig_pdf_trunc",32,2.4,4);

  for(auto i_bkg_pdf = 49; i_bkg_pdf <= 80; ++i_bkg_pdf) {
    h_bkg_pdf_trunc->Fill(h_pdf_bkg->GetXaxis()->GetBinCenter(i_bkg_pdf),
                          h_pdf_bkg->GetBinContent(i_bkg_pdf));
    h_sig_pdf_trunc->Fill(h_pdf_sig->GetXaxis()->GetBinCenter(i_bkg_pdf),
                          h_pdf_sig->GetBinContent(i_bkg_pdf));
  }

  // Background and signal PDFs
  std::vector<double> pdf_B;
  std::vector<double> pdf_S;

  // double norm_factor_2nu = 1. / 25;
  // double norm_factor_0nu = 1;
  double norm_factor_2nu = h_pdf_bkg->GetEntries()/1e6 / 25;
  double norm_factor_0nu = h_pdf_sig->GetEntries()/1e6;

  for(auto i_pdf = 1; i_pdf <= h_bkg_pdf_trunc->GetNbinsX(); ++i_pdf) {
    pdf_B.push_back(norm_factor_2nu*h_bkg_pdf_trunc->GetBinContent(i_pdf));
    pdf_S.push_back(norm_factor_0nu*h_sig_pdf_trunc->GetBinContent(i_pdf));
  }

  h_bkg_pdf_trunc->Scale(1./h_bkg_pdf_trunc->Integral(1,h_bkg_pdf_trunc->GetNbinsX()));
  h_sig_pdf_trunc->Scale(1./h_sig_pdf_trunc->Integral(1,h_sig_pdf_trunc->GetNbinsX()));

  h_bkg_pdf_trunc->Write();
  h_sig_pdf_trunc->Write();

  TRandom *rdm = new TRandom();
  rdm->SetSeed(0);

  double p = 0;

  std::vector<double> s_90;

  TH1 *h_data = new TH1I("data","data",32,2.4,4);

  // for(unsigned int i_pseudo =0; i_pseudo<200; ++i_pseudo) {
  for(unsigned int i_pseudo =0; i_pseudo<1; ++i_pseudo) {
    std::cout << "pseudo_experiment n°" << i_pseudo << std::endl;

    h_data->Reset();

    //rescale signal if pseudo signal

    pseudo_generator(h_bkg_pdf_trunc,h_data,n_2nu_events,0);
    // pseudo_generator(h_sig_pdf_trunc,h_data,n_0nu_events,0);

    h_data->Write();

    std::vector<double> data;

    for(auto i_h = 1; i_h <= h_data->GetNbinsX(); ++i_h)
      data.push_back(h_data->GetBinContent(i_h));

    // std::cout << " checking proba " << proba_bayes(0,1978,pdf_S,pdf_B,data)*p0_S(0)*p0_B(1978) << std::endl;

    // proba_bayes(0,2,pdf_S,pdf_B,data)*p0_S(0)*p0_B(2);

    // return 0;

    double tot_prob = 0;
    //Loop over S90
    for(double s = 0; s<50; ++s) {
    // for(double s = 0; s<=0; ++s) {
      std::cout << " Computing proba for signal hypothesis S = " << s << std::endl;
      p = 0;
      // for(double b=0; b<=10000; b += 100) {
      // for(double b=1978; b<=1978; b += 100) {
      for(double b=0; b<=100; ++b) {
      // for(double b=0; b<=0; ++b) {
        // std::cout << "    ... Integrating proba for background B=" << b << std::endl;

        double num = proba_bayes(s,b,pdf_S,pdf_B,data)*p0_S(s)*p0_B(b);

        // continue;

        // if(std::isinf(num))
        //   continue;

        // std::cout << " num proba bayes  s="<< s << "  b=" << b << "    " << num << std::endl;
        double denum = 0;
        for(double i_s = 0; i_s<=50; ++i_s)
          // for(double i_b = 0; i_b<=10000; i_b += 100) {
          for(double i_b = 0; i_b<=100; ++i_b) {
              denum += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
          }
        // std::cout << " denum proba bayes " << denum << std::endl;
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

  h_pdf_sig->SetAxisRange(2.4,4);
  h_pdf_bkg->SetAxisRange(2.4,4);

  h_pdf_bkg->Write();
  // h_data->Write();

  h_s90->Write();
  g_p->Write();
  f->Close();
  // std::cout << "proba is " << p << std::endl;
  return 0;
}
