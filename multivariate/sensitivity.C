#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLine.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "config_sensitivity.h"

void sensitivity()
{
  TFile * f_bdt = TFile::Open("bdt_scores/source/B/bdt_scores.root");

  TH1F *h_0nu_bdt = (TH1F*)f_bdt->Get("0nu");
  TH1F *h_2nu_bdt = (TH1F*)f_bdt->Get("2nu");
  TH1F *h_tl208_bdt = (TH1F*)f_bdt->Get("tl208");
  TH1F *h_bi214_bdt = (TH1F*)f_bdt->Get("bi214");
  TH1F *h_radon_bdt = (TH1F*)f_bdt->Get("radon");
  THStack *h_stack_bdt_norm = (THStack*)f_bdt->Get("hs");

  TFile * f_roi = TFile::Open("spectra/source/A/spectra_Ecut.root");

  TH1F *h_0nu_roi = (TH1F*)f_roi->Get("0nu");
  TH1F *h_2nu_roi = (TH1F*)f_roi->Get("2nu");
  TH1F *h_tl208_roi = (TH1F*)f_roi->Get("tl208");
  TH1F *h_bi214_roi = (TH1F*)f_roi->Get("bi214");
  TH1F *h_radon_roi = (TH1F*)f_roi->Get("radon");
  THStack *h_stack_roi_norm = (THStack*)f_roi->Get("hs");

  TGraphErrors * g_eff_0nu_roi = new TGraphErrors();
  TGraphErrors * g_eff_2nu_roi = new TGraphErrors();
  TGraphErrors * g_eff_tl208_roi = new TGraphErrors();
  TGraphErrors * g_eff_bi214_roi = new TGraphErrors();
  TGraphErrors * g_eff_radon_roi = new TGraphErrors();

  TGraphErrors * g_eff_0nu_bdt = new TGraphErrors();
  TGraphErrors * g_eff_2nu_bdt = new TGraphErrors();
  TGraphErrors * g_eff_tl208_bdt = new TGraphErrors();
  TGraphErrors * g_eff_bi214_bdt = new TGraphErrors();
  TGraphErrors * g_eff_radon_bdt = new TGraphErrors();

  TFile *f_output= new TFile("sensitivity.root","RECREATE");

  TGraphErrors * bdt_significance = new TGraphErrors();
  bdt_significance->SetMarkerSize(2);
  bdt_significance->SetMarkerColor(kBlue);
  bdt_significance->SetDrawOption("AP");
  bdt_significance->SetTitle("BDT;BDT score;Significance");

  TGraphErrors * bdt_halflife = new TGraphErrors();
  bdt_halflife->SetMarkerSize(2);
  bdt_halflife->SetMarkerColor(kBlue);
  bdt_halflife->SetDrawOption("AP");
  bdt_halflife->SetTitle("BDT;BDT score;Halflife");

  unsigned int nbins = h_0nu_bdt->GetNbinsX();

  double best_halflife_limit_bdt = 0;
  double Ncount_bdt = 0;
  double N_2nu_bdt = 0;
  double N_tl_bdt = 0;
  double N_bi_bdt = 0;
  double N_radon_bdt = 0;

  int best_bdt_cut_bin = 0;
  double best_bdt_cut_score = 0;

  double tmp_eff = 0;

  for (unsigned int i = 1; i <= nbins; ++i) {

    double eff_0nu = h_0nu_bdt->Integral(i,nbins); // normalized to 1

    double eff_2nu = h_2nu_bdt->Integral(i,nbins); // normalized to 1
    double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
    double eff_tl208 = h_tl208_bdt->Integral(i,nbins);
    double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
    double eff_bi214 = h_bi214_bdt->Integral(i,nbins);
    double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
    double eff_radon = h_radon_bdt->Integral(i,nbins);
    double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

    double S = 0;
    double B = 0;

    S = eff_0nu;
    B += N_2nu;
    B += N_tl208;
    B += N_bi214;
    B += N_radon;

    g_eff_0nu_bdt->SetPoint(i-1,h_0nu_bdt->GetBinLowEdge(i),h_0nu_bdt->Integral(i,nbins)*conf_sens::eff_0nu_2e);
    g_eff_2nu_bdt->SetPoint(i-1,h_2nu_bdt->GetBinLowEdge(i),h_2nu_bdt->Integral(i,nbins)*conf_sens::eff_2nu_2e);
    g_eff_tl208_bdt->SetPoint(i-1,h_tl208_bdt->GetBinLowEdge(i),h_tl208_bdt->Integral(i,nbins)*conf_sens::eff_tl208_2e);
    g_eff_bi214_bdt->SetPoint(i-1,h_bi214_bdt->GetBinLowEdge(i),h_bi214_bdt->Integral(i,nbins)*conf_sens::eff_bi214_2e);
    g_eff_radon_bdt->SetPoint(i-1,h_radon_bdt->GetBinLowEdge(i),h_radon_bdt->Integral(i,nbins)*conf_sens::eff_radon_2e);

    if(S+B != 0)
      bdt_significance->SetPoint(i-1,h_0nu_bdt->GetBinLowEdge(i),S/sqrt(S+B));
    else
      bdt_significance->SetPoint(i-1,h_0nu_bdt->GetBinLowEdge(i),0);

    double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
    // double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214);
    // double N_excluded = get_number_of_excluded_events(N_2nu);
    double halflife = eff_0nu * conf_sens::eff_0nu_2e * conf_sens::k_sens / N_excluded;

    if(best_halflife_limit_bdt<halflife) {
      Ncount_bdt = N_excluded;
      N_2nu_bdt = N_2nu;
      N_tl_bdt = N_tl208;
      N_bi_bdt = N_bi214;
      N_radon_bdt = N_radon;
      best_bdt_cut_bin = i;
      best_bdt_cut_score = h_0nu_bdt->GetBinLowEdge(i);

    }

    best_halflife_limit_bdt = std::max(best_halflife_limit_bdt, halflife);
    bdt_halflife->SetPoint(i-1,h_0nu_bdt->GetBinLowEdge(i),halflife);
  }

  TGraphErrors * roi_significance = new TGraphErrors();
  roi_significance->SetMarkerSize(2);
  roi_significance->SetMarkerColor(kRed);
  roi_significance->SetDrawOption("AP");
  roi_significance->SetTitle("ROI;Energy;Significance");

  TGraphErrors * roi_halflife = new TGraphErrors();
  roi_halflife->SetMarkerSize(2);
  roi_halflife->SetMarkerColor(kBlue);
  roi_halflife->SetDrawOption("AP");
  roi_halflife->SetTitle("ROI;ROI energy;Halflife");

  double best_halflife_limit_roi = 0;
  double Ncount_roi = 0;
  double N_2nu_roi = 0;
  double N_tl_roi = 0;
  double N_bi_roi = 0;
  double N_radon_roi = 0;
  int lower_window_bin = 0;
  int upper_window_bin = 0;
  double lower_window_energy = 0;
  double upper_window_energy = 0;

  unsigned int nbins_roi = h_0nu_roi->GetNbinsX();

  for (unsigned int i = 1; i <= h_0nu_roi->GetNbinsX(); ++i) {

    double eff_0nu = h_0nu_roi->Integral(i,nbins_roi); // normalized to 1

    double eff_2nu = h_2nu_roi->Integral(i,nbins_roi); // normalized to 1
    double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
    double eff_tl208 = h_tl208_roi->Integral(i,nbins_roi);
    double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
    double eff_bi214 = h_bi214_roi->Integral(i,nbins_roi);
    double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
    double eff_radon = h_radon_roi->Integral(i,nbins_roi);
    double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

    double S = 0;
    double B = 0;
    S  = eff_0nu;
    B += N_2nu;
    B += N_tl208;
    B += N_bi214;
    B += N_radon;

    g_eff_0nu_roi->SetPoint(i-1,h_0nu_roi->GetBinLowEdge(i),h_0nu_roi->Integral(i,nbins_roi)*conf_sens::eff_0nu_2e);
    g_eff_2nu_roi->SetPoint(i-1,h_2nu_roi->GetBinLowEdge(i),h_2nu_roi->Integral(i,nbins_roi)*conf_sens::eff_2nu_2e);
    g_eff_tl208_roi->SetPoint(i-1,h_tl208_roi->GetBinLowEdge(i),h_tl208_roi->Integral(i,nbins_roi)*conf_sens::eff_tl208_2e);
    g_eff_bi214_roi->SetPoint(i-1,h_bi214_roi->GetBinLowEdge(i),h_bi214_roi->Integral(i,nbins_roi)*conf_sens::eff_bi214_2e);
    g_eff_radon_roi->SetPoint(i-1,h_radon_roi->GetBinLowEdge(i),h_radon_roi->Integral(i,nbins_roi)*conf_sens::eff_radon_2e);

    if(S+B != 0)
      roi_significance->SetPoint(i-1,h_0nu_roi->GetBinLowEdge(i),S/sqrt(S+B));
    else
      roi_significance->SetPoint(i-1,h_0nu_roi->GetBinLowEdge(i),0);

    double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
    // double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214);
    // double N_excluded = get_number_of_excluded_events(N_2nu);
    double halflife = eff_0nu * conf_sens::eff_0nu_2e * conf_sens::k_sens / N_excluded;

    if(best_halflife_limit_roi<halflife) {
      Ncount_roi = N_excluded;
      N_2nu_roi = N_2nu;
      N_tl_roi = N_tl208;
      N_bi_roi = N_bi214;
      N_radon_roi = N_radon;

      // std::cout << " debug " << h_0nu_roi->GetBinLowEdge(i) << "  " <<  N_2nu_roi << std::endl;
    }
    best_halflife_limit_roi = std::max(best_halflife_limit_roi, halflife);

    roi_halflife->SetPoint(i-1,h_0nu_roi->GetBinLowEdge(i),halflife);
  }

  //2D window optimization

  TH2F * halflife_2d_roi = new TH2F("halflife_2d_roi","halflife_2d_roi",h_0nu_roi->GetNbinsX(),0,4000,h_0nu_roi->GetNbinsX(),0,4000);

  double best_halflife_limit_roi_2d = 0;
  double Ncount_roi_2d = 0;
  double N_2nu_roi_2d = 0;
  double N_tl_roi_2d = 0;
  double N_bi_roi_2d = 0;
  double N_radon_roi_2d = 0;

  for (unsigned int i = 1; i < h_0nu_roi->GetNbinsX(); ++i) {
    for(unsigned int j = i+1; j <= h_0nu_roi->GetNbinsX(); ++j) {
      double eff_0nu = h_0nu_roi->Integral(i,j); // normalized to 1
      double eff_2nu = h_2nu_roi->Integral(i,j); // normalized to 1
      double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
      double eff_tl208 = h_tl208_roi->Integral(i,j);
      double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
      double eff_bi214 = h_bi214_roi->Integral(i,j);
      double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
      double eff_radon = h_radon_roi->Integral(i,j);
      double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

      double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
      // double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214);
      // double N_excluded = get_number_of_excluded_events(N_2nu);
      double halflife = eff_0nu * conf_sens::eff_0nu_2e * conf_sens::k_sens / N_excluded;

      if(best_halflife_limit_roi_2d<halflife) {
        Ncount_roi_2d = N_excluded;
        N_2nu_roi_2d = N_2nu;
        N_tl_roi_2d = N_tl208;
        N_bi_roi_2d = N_bi214;
        N_radon_roi_2d = N_radon;
        lower_window_bin = i;
        upper_window_bin = j;
        lower_window_energy = h_0nu_roi->GetBinLowEdge(i);
        upper_window_energy = h_0nu_roi->GetBinLowEdge(j);
      }
      best_halflife_limit_roi_2d = std::max(best_halflife_limit_roi_2d, halflife);
      halflife_2d_roi->SetBinContent(i,j,halflife);
    }
  }
  tmp_eff = h_0nu_roi->Integral(lower_window_bin,upper_window_bin)*conf_sens::eff_0nu_2e;

  TH2F * halflife_2d_bdt = new TH2F("halflife_2d_bdt","halflife_2d_bdt",h_0nu_bdt->GetNbinsX(), h_0nu_bdt->GetBinLowEdge(1), h_0nu_bdt->GetBinLowEdge(h_0nu_bdt->GetNbinsX()+1),
                                    h_0nu_bdt->GetNbinsX(), h_0nu_bdt->GetBinLowEdge(1), h_0nu_bdt->GetBinLowEdge(h_0nu_bdt->GetNbinsX()+1));

  double best_halflife_limit_bdt_2d = 0;
  double Ncount_bdt_2d = 0;
  double N_2nu_bdt_2d = 0;
  double N_tl_bdt_2d = 0;
  double N_bi_bdt_2d = 0;
  double N_radon_bdt_2d = 0;
  int lower_window_bin_bdt = 0;
  int upper_window_bin_bdt = 0;
  double lower_window_energy_bdt = 0;
  double upper_window_energy_bdt = 0;

  for (unsigned int i = 1; i < h_0nu_bdt->GetNbinsX(); ++i) {
    for(unsigned int j = i+1; j <= h_0nu_bdt->GetNbinsX(); ++j) {
      double eff_0nu = h_0nu_bdt->Integral(i,j); // normalized to 1
      double eff_2nu = h_2nu_bdt->Integral(i,j); // normalized to 1
      double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
      double eff_tl208 = h_tl208_bdt->Integral(i,j); // normalized to 1
      double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
      double eff_bi214 = h_bi214_bdt->Integral(i,j);
      double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
      double eff_radon = h_radon_bdt->Integral(i,j);
      double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

      // double N_excluded = get_number_of_excluded_events(N_2nu);
      double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
      // double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214);
      double halflife = eff_0nu * conf_sens::eff_0nu_2e * conf_sens::k_sens / N_excluded;

      if(best_halflife_limit_bdt_2d<halflife) {
        Ncount_bdt_2d = N_excluded;
        N_2nu_bdt_2d = N_2nu;
        N_tl_bdt_2d = N_tl208;
        N_bi_bdt_2d = N_bi214;
        N_radon_bdt_2d = N_radon;
        lower_window_bin_bdt = i;
        upper_window_bin_bdt = j;
        lower_window_energy_bdt = h_0nu_bdt->GetBinLowEdge(i);
        upper_window_energy_bdt = h_0nu_bdt->GetBinLowEdge(j);
      }
      best_halflife_limit_bdt_2d = std::max(best_halflife_limit_bdt_2d, halflife);
      halflife_2d_bdt->SetBinContent(i,j,halflife);
    }
  }

  // bdt_significance->Draw("AP");
  bdt_significance->SetName("bdt_significance_norm");
  roi_significance->SetName("roi_significance_norm");
  bdt_significance->Write();
  roi_significance->Write();

  h_stack_bdt_norm->SetName("hs_bdt_norm");
  h_stack_roi_norm->SetName("hs_roi_norm");
  h_stack_bdt_norm->Write();
  h_stack_roi_norm->Write();

  TDirectory *isotope = f_output->mkdir("isotope");
  isotope->cd();

  h_0nu_bdt->SetName("h_0nu_bdt_norm");
  h_2nu_bdt->SetName("h_2nu_bdt_norm");
  h_tl208_bdt->SetName("h_tl208_bdt_norm");
  h_bi214_bdt->SetName("h_bi214_bdt_norm");
  h_radon_bdt->SetName("h_radon_bdt_norm");
  h_0nu_bdt->Write();
  h_2nu_bdt->Write();
  h_tl208_bdt->Write();
  h_bi214_bdt->Write();
  h_radon_bdt->Write();

  h_0nu_roi->SetName("h_0nu_roi_norm");
  h_2nu_roi->SetName("h_2nu_roi_norm");
  h_tl208_roi->SetName("h_tl208_roi_norm");
  h_bi214_roi->SetName("h_bi214_roi_norm");
  h_radon_roi->SetName("h_radon_roi_norm");
  h_0nu_roi->Write();
  h_2nu_roi->Write();
  h_tl208_roi->Write();
  h_bi214_roi->Write();
  h_radon_roi->Write();

  THStack *hs_bdt_count = new THStack("hs_bdt_count","Stacked spectra");
  h_0nu_bdt->Scale(10);
  h_2nu_bdt->Scale(conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu );
  h_tl208_bdt->Scale(conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity );
  h_bi214_bdt->Scale(conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity );
  h_radon_bdt->Scale(conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity );

  hs_bdt_count->Add(h_tl208_bdt);
  hs_bdt_count->Add(h_bi214_bdt);
  hs_bdt_count->Add(h_radon_bdt);
  hs_bdt_count->Add(h_2nu_bdt);
  hs_bdt_count->Add(h_0nu_bdt);

  f_output->cd();
  hs_bdt_count->Write();

  THStack *hs_roi_count = new THStack("hs_roi_count","Stacked spectra");
  h_0nu_roi->Scale(10);
  h_2nu_roi->Scale(conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu);
  h_tl208_roi->Scale(conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity);
  h_bi214_roi->Scale(conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity);
  h_radon_roi->Scale(conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity);

  hs_roi_count->Add(h_tl208_roi);
  hs_roi_count->Add(h_bi214_roi);
  hs_roi_count->Add(h_radon_roi);
  hs_roi_count->Add(h_2nu_roi);
  hs_roi_count->Add(h_0nu_roi);

  hs_roi_count->Write();

  isotope->cd();

  h_0nu_bdt->SetName("h_0nu_bdt_count");
  h_2nu_bdt->SetName("h_2nu_bdt_count");
  h_tl208_bdt->SetName("h_tl208_bdt_count");
  h_bi214_bdt->SetName("h_bi214_bdt_count");
  h_radon_bdt->SetName("h_radon_bdt_count");
  h_0nu_bdt->Write();
  h_2nu_bdt->Write();
  h_tl208_bdt->Write();
  h_bi214_bdt->Write();
  h_radon_bdt->Write();

  h_0nu_roi->SetName("h_0nu_roi_count");
  h_2nu_roi->SetName("h_2nu_roi_count");
  h_tl208_roi->SetName("h_tl208_roi_count");
  h_bi214_roi->SetName("h_bi214_roi_count");
  h_radon_roi->SetName("h_radon_roi_count");
  h_0nu_roi->Write();
  h_2nu_roi->Write();
  h_tl208_roi->Write();
  h_bi214_roi->Write();
  h_radon_roi->Write();

  std::cout << " Best halflives " << std::endl;
  std::cout << "     BDT    :  "  << best_halflife_limit_bdt << "  ,   N_bg = " << Ncount_bdt
            << "  (N_2nu = " << N_2nu_bdt << ", N_tl = " << N_tl_bdt << ", N_bi = " << N_bi_bdt <<", N_radon = " << N_radon_bdt << ")" << std::endl;
  // std::cout << "     BDT    :  "  << best_halflife_limit_bdt << "  ,   N_bg = " << Ncount_bdt
  //           << "  (N_2nu = " << N_2nu_bdt << std::endl; // << ", N_radon = " << N_radon_bdt << ")" << std::endl;
  std::cout << "      Cut is  " << best_bdt_cut_bin << "  or  " << best_bdt_cut_score << std::endl;
  std::cout << "      Signal efficiency : " << h_0nu_bdt->Integral(best_bdt_cut_bin,h_0nu_bdt->GetNbinsX())/h_0nu_bdt->Integral(0, h_0nu_bdt->GetNbinsX()) << std::endl;// << "   Background efficiency : " << h_2nu_bdt->Integral(best_bdt_cut_bin,h_2nu_bdt->GetNbinsX())/h_2nu_bdt->Integral(0, h_2nu_bdt->GetNbinsX()) << std::endl;
  std::cout << "     BDT 2D  :  "  << best_halflife_limit_bdt_2d << "  ,   N_bg = " << Ncount_bdt_2d
            << "  (N_2nu = " << N_2nu_bdt_2d << ", N_tl = " << N_tl_bdt_2d << ", N_bi = " << N_bi_bdt_2d << ", N_radon = " << N_radon_bdt_2d << ")" << std::endl;
  // std::cout << "     BDT 2D :  "  << best_halflife_limit_bdt_2d << "  ,   N_bg = " << Ncount_bdt_2d
  //           << "  (N_2nu = " << N_2nu_bdt_2d << std::endl; //", N_radon = " << N_radon_bdt_2d << ")" << std::endl;
  std::cout << "      Lower window is  " << lower_window_bin_bdt << "  or  " << lower_window_energy_bdt << std::endl;
  std::cout << "      Upper window is  " << upper_window_bin_bdt << "  or  " << upper_window_energy_bdt << std::endl;
  std::cout << "      Signal efficiency : " << h_0nu_bdt->Integral(lower_window_bin_bdt,upper_window_bin_bdt)/h_0nu_bdt->Integral(0, h_0nu_bdt->GetNbinsX()) << std::endl; // "   Background efficiency : " << h_2nu_bdt->I
  std::cout << "     ROI    :  "  << best_halflife_limit_roi << "  ,   N_bg = " << Ncount_roi
            << "  (N_2nu = " << N_2nu_roi << ", N_tl = " << N_tl_roi << ", N_bi = " << N_bi_roi << ", N_radon = " << N_radon_roi << ")" << std::endl;
  // std::cout << "     ROI    :  "  << best_halflife_limit_roi << "  ,   N_bg = " << Ncount_roi
  //           << "  (N_2nu = " << N_2nu_roi << std::endl;// ", N_radon = " << N_radon_roi << ")" << std::endl;
  std::cout << "     ROI 2D :  "  << best_halflife_limit_roi_2d << "  ,   N_bg = " << Ncount_roi_2d
            << "  (N_2nu = " << N_2nu_roi_2d << ", N_tl = " << N_tl_roi_2d << ", N_bi = " << N_bi_roi_2d << ", N_radon = " << N_radon_roi_2d << ")" << std::endl;
  // std::cout << "     ROI 2D :  "  << best_halflife_limit_roi_2d << "  ,   N_bg = " << Ncount_roi_2d
  //           << "  (N_2nu = " << N_2nu_roi_2d << std::endl; //", N_radon = " << N_radon_roi_2d << ")" << std::endl;
  std::cout << "      Lower window is  " << lower_window_bin << "  or  " << lower_window_energy << std::endl;
  std::cout << "      Upper window is  " << upper_window_bin << "  or  " << upper_window_energy << std::endl;
  std::cout << "      Signal efficiency : " << h_0nu_roi->Integral(lower_window_bin,upper_window_bin)/h_0nu_roi->Integral(0, h_0nu_roi->GetNbinsX()) << std::endl; // "   Background efficiency : " << h_2nu_roi->Integral(lower_window_bin,upper_window_bin)/h_2nu_roi->Integral(0, h_2nu_roi->GetNbinsX()) << std::endl;
  // std::cout << " eff " << tmp_eff <<std::endl;

  f_output->cd();
  bdt_halflife->SetName("bdt_halflife");
  roi_halflife->SetName("roi_halflife");
  halflife_2d_roi->SetName("roi_halflife_2d");
  bdt_halflife->Write();
  roi_halflife->Write();
  halflife_2d_roi->Write();

  TDirectory *efficiencies = f_output->mkdir("efficiencies");
  efficiencies->cd();

  g_eff_0nu_bdt->SetName("g_eff_0nu_bdt");
  g_eff_0nu_bdt->SetTitle("0#nu efficiency;Energy;Efficiency");
  g_eff_0nu_bdt->SetDrawOption("AL");
  g_eff_0nu_bdt->SetLineColor(kRed);
  g_eff_0nu_bdt->SetLineWidth(2);
  g_eff_0nu_bdt->Write();

  g_eff_2nu_bdt->SetName("g_eff_2nu_bdt");
  g_eff_2nu_bdt->SetTitle("2#nu efficiency;Energy;Efficiency");
  g_eff_2nu_bdt->SetDrawOption("AL");
  g_eff_2nu_bdt->SetLineColor(kBlue);
  g_eff_2nu_bdt->SetLineWidth(2);
  g_eff_2nu_bdt->Write();

  g_eff_tl208_bdt->SetName("g_eff_tl208_bdt");
  g_eff_tl208_bdt->SetTitle("^{208}Tl efficiency;Energy;Efficiency");
  g_eff_tl208_bdt->SetDrawOption("AL");
  g_eff_tl208_bdt->SetLineColor(kGreen+1);
  g_eff_tl208_bdt->SetLineWidth(2);
  g_eff_tl208_bdt->Write();

  g_eff_bi214_bdt->SetName("g_eff_bi214_bdt");
  g_eff_bi214_bdt->SetTitle("^{214}Bi efficiency;Energy;Efficiency");
  g_eff_bi214_bdt->SetDrawOption("AL");
  g_eff_bi214_bdt->SetLineColor(kOrange-3);
  g_eff_bi214_bdt->SetLineWidth(2);
  g_eff_bi214_bdt->Write();

  g_eff_radon_bdt->SetName("g_eff_radon_bdt");
  g_eff_radon_bdt->SetTitle("Radon efficiency;Energy;Efficiency");
  g_eff_radon_bdt->SetDrawOption("AL");
  g_eff_radon_bdt->SetLineColor(kMagenta);
  g_eff_radon_bdt->SetLineWidth(2);
  g_eff_radon_bdt->Write();

  g_eff_0nu_roi->SetName("g_eff_0nu_roi");
  g_eff_0nu_roi->SetTitle("0#nu efficiency;Energy;Efficiency");
  g_eff_0nu_roi->SetDrawOption("AL");
  g_eff_0nu_roi->SetLineColor(kRed);
  g_eff_0nu_roi->SetLineWidth(2);
  g_eff_0nu_roi->Write();

  g_eff_2nu_roi->SetName("g_eff_2nu_roi");
  g_eff_2nu_roi->SetTitle("2#nu efficiency;Energy;Efficiency");
  g_eff_2nu_roi->SetDrawOption("AL");
  g_eff_2nu_roi->SetLineColor(kBlue);
  g_eff_2nu_roi->SetLineWidth(2);
  g_eff_2nu_roi->Write();

  g_eff_tl208_roi->SetName("g_eff_tl208_roi");
  g_eff_tl208_roi->SetTitle("^{208}Tl efficiency;Energy;Efficiency");
  g_eff_tl208_roi->SetDrawOption("AL");
  g_eff_tl208_roi->SetLineColor(kGreen+1);
  g_eff_tl208_roi->SetLineWidth(2);
  g_eff_tl208_roi->Write();

  g_eff_bi214_roi->SetName("g_eff_bi214_roi");
  g_eff_bi214_roi->SetTitle("^{214}Bi efficiency;Energy;Efficiency");
  g_eff_bi214_roi->SetDrawOption("AL");
  g_eff_bi214_roi->SetLineColor(kOrange-3);
  g_eff_bi214_roi->SetLineWidth(2);
  g_eff_bi214_roi->Write();

  g_eff_radon_roi->SetName("g_eff_radon_roi");
  g_eff_radon_roi->SetTitle("Radon efficiency;Energy;Efficiency");
  g_eff_radon_roi->SetDrawOption("AL");
  g_eff_radon_roi->SetLineColor(kMagenta);
  g_eff_radon_roi->SetLineWidth(2);
  g_eff_radon_roi->Write();

  // h_0nu_bdt->Draw("");
  // h_2nu_bdt->Draw("same");
  // h_tl208_bdt->Draw("same");
  // h_bi214_bdt->Draw("same");
  // h_radon_bdt->Draw("same");
  // leg->Draw("same");

  // gStyle->SetTitleFontSize(0.08);
  // gPad->SetLogy();
  // // gPad->SetPad(0.1,0.2,0.9,0.95);
  // gPad->SetBottomMargin(0.155);
  // gPad->SetTopMargin(1);
  // gPad->SetRightMargin(1000);

  return;
}
