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

  // Getting the expected number of background and signal events
  //efficiency = subrange of pdf [2.45;4]MeV * eff MC * MC>2MeV correction
  double eff_2nu =  h_pdf_bkg->Integral(50,80) * h_pdf_bkg->GetEntries()/1e6 / 25;
  int n_2nu_events = eff_2nu * mass * exposure_y * Na * log(2) / M_Se / halflife_2nu;
  // int n_2nu_events = 10;

  int n_0nu_events = 0;

  std::cout << "Expected number of background events " << n_2nu_events << std::endl;

  TFile *f= new TFile("output.root", "RECREATE");

  // Use extra histograms to hold the pdfs normalized on the specified subrange
  TH1F *h_bkg_pdf_trunc = new TH1F("bkg_pdf_trunc","bkg_pdf_trunc",31,2.45,4);
  TH1F *h_sig_pdf_trunc = new TH1F("sig_pdf_trunc","sig_pdf_trunc",31,2.45,4);

  for(auto i_pdf = 50; i_pdf <= 80; ++i_pdf) {
    h_bkg_pdf_trunc->SetBinContent(i_pdf-49,
                                   h_pdf_bkg->GetBinContent(i_pdf));
    h_bkg_pdf_trunc->SetBinError(i_pdf-49,
                                 h_pdf_bkg->GetBinError(i_pdf));
    h_sig_pdf_trunc->SetBinContent(i_pdf-49,
                                   h_pdf_sig->GetBinContent(i_pdf));
    h_sig_pdf_trunc->SetBinError(i_pdf-49,
                                 h_pdf_sig->GetBinError(i_pdf));
  }

  h_bkg_pdf_trunc->Scale(1./h_bkg_pdf_trunc->Integral(1,h_bkg_pdf_trunc->GetNbinsX()));
  h_sig_pdf_trunc->Scale(1./h_sig_pdf_trunc->Integral(1,h_sig_pdf_trunc->GetNbinsX()));

  h_bkg_pdf_trunc->Write();
  h_sig_pdf_trunc->Write();

  // Background and signal PDFs
  std::vector<double> pdf_B;
  std::vector<double> pdf_S;

  for(auto i_pdf = 1; i_pdf <= h_bkg_pdf_trunc->GetNbinsX(); ++i_pdf) {
    pdf_B.push_back(h_bkg_pdf_trunc->GetBinContent(i_pdf));
    pdf_S.push_back(h_sig_pdf_trunc->GetBinContent(i_pdf));
  }

  int seed = 0;
  TRandom *rdm = new TRandom();
  rdm->SetSeed(seed);

  // Vector holding the 90% rejected number of signal events for all pseudo-experiments
  std::vector<double> s_90;

  // Histogram holding the different pseudo-experiment (reset for each one)
  TH1 *h_data = new TH1I("data","data",31,2.45,4);

  // for(unsigned int i_pseudo =0; i_pseudo<200; ++i_pseudo) {
  for(unsigned int i_pseudo =0; i_pseudo<1; ++i_pseudo) {
    std::cout << "pseudo_experiment n°" << i_pseudo << std::endl;

    // Reset the pseudo-data from one pseudo-experiment to the other
    h_data->Reset();

    // Generate pseudo experiments according to the signal and background pdfs
    pseudo_generator(h_bkg_pdf_trunc,h_data,n_2nu_events,seed);
    pseudo_generator(h_sig_pdf_trunc,h_data,n_0nu_events,seed);

    //Save the first pseudo experiment
    if(i_pseudo==0)
      h_data->Write();

    std::vector<double> data;

    //Save the pseudo-experiment in a vector for convenience
    for(auto i_h = 1; i_h <= h_data->GetNbinsX(); ++i_h)
      data.push_back(h_data->GetBinContent(i_h));

    // For limit setting : integrate the total probability
    // for an ascending assumed number of signal events (then check when it exceeds 90%)
    double tot_prob = 0;

    // Loop over the assumed number of signal events. Max is set to 49
    for(double s = 0; s<=50; ++s) {

      std::cout << " -- Computing probability for signal hypothesis S = " << s << std::endl;

      //The probability for this signal hypothesis
      double p = 0;

      //Integrate over all relevant background values (<=100 given the range studied)
      for(double b=0; b<=100; ++b) {

        //Numerator of the probability
        double num = proba_bayes(s,b,pdf_S,pdf_B,data)*p0_S(s)*p0_B(b);
        // std::cout << " num proba bayes  s="<< s << "  b=" << b << "    " << num << std::endl;

        //Denominator of the probability
        double denom = 0;
        for(double i_s = 0; i_s<=50; ++i_s)
          for(double i_b = 0; i_b<=100; ++i_b) {
              denom += proba_bayes(i_s,i_b,pdf_S,pdf_B,data)*p0_S(i_s)*p0_B(i_b);
          }
        // std::cout << " denom proba bayes " << denom << std::endl;

        if(denom!=0)
          p+=num/denom;
      }

      //Save the probability as a function of the number of signal events
      g_p->SetPoint(s,s,p);
      std::cout << "For the hypothesis  S = " << s << " :  Probability = " << p << std::endl;
      tot_prob += p;

      //For limit setting when signal is zero
      if(tot_prob>=0.9) {
        std::cout << "---" << "Lower limit on signal s_90 =" << s << std::endl;
        s_90.push_back(s);
        h_s90->Fill(s);
        break;
      }

    }

  } //end pseudo experiments loop

  h_pdf_sig->SetAxisRange(2.45,4);
  h_pdf_bkg->SetAxisRange(2.45,4);

  h_pdf_bkg->Write();

  h_s90->Write();
  g_p->Write();
  f->Close();

  return 0;
}
