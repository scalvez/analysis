#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>

// Continuous FeldmanCousin functions computed by M. Bongrand
// <bongrand@lal.in2p3.fr>
double get_number_of_excluded_events(const double number_of_events_)
{
  double number_of_excluded_events = 0.0;
  if (number_of_events_ < 29.0)
    {
      double x = number_of_events_;
      number_of_excluded_events =
        2.5617 + 0.747661 * x - 0.0666176 * std::pow(x,2)
        + 0.00432457 * std::pow(x,3) - 0.000139343 * std::pow(x,4)
        + 1.71509e-06 * std::pow(x,5);
    }
  else
    {
      number_of_excluded_events = 1.64 * std::sqrt(number_of_events_);
    }
  return number_of_excluded_events;
}


void quick_fc() {

  const double mass = 7.;
  const double exposure_sec = 2.5 * 3.14e7;
  const double exposure_y = 2.5;
  // const double mass = 100.;
  // const double exposure_sec = 5 * 3.14e7;
  // const double exposure_y = 5;
  const double tracker_volume = 15.3; //m3
  const double halflife_2nu = 9e19; // years;
  const double Na = 6.022e23;
  const double M_Se = 0.082; //kg/mol more like 81.6 if we enrich at 90%
  const double const_se = log(2) * Na / M_Se; // years;

  double eff_0nu = 0.238627;
  double eff_2nu = 0.0077340;
  double eff_tl208 = 0.00092096;
  double eff_bi214 = 0.00141587;
  double eff_radon = 8.616e-05;

  TString output_file = "sensitivity.root";

  TFile *f_0nu = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/fred/0nu_pdf.root");
  TFile *f_2nu = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/fred/2nu_pdf.root");
  TFile *f_tl208 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/fred/tl208_pdf.root");
  TFile *f_bi214 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/fred/bi214_pdf.root");
  TFile *f_radon = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/fred/radon_pdf.root");

  TH1 *h_0nu = (TH1*)f_0nu->Get("2e_electrons_energy_sum");
  TH1 *h_2nu = (TH1*)f_2nu->Get("2e_electrons_energy_sum");
  TH1 *h_tl208 = (TH1*)f_tl208->Get("2e_electrons_energy_sum");
  TH1 *h_bi214 = (TH1*)f_bi214->Get("2e_electrons_energy_sum");
  TH1 *h_radon = (TH1*)f_radon->Get("2e_electrons_energy_sum");

  TFile *f_output= new TFile(output_file, "RECREATE");

  TH1 *h_eff_0nu = new TH1F("eff_0nu","eff_0nu",100,0,5);
  TH1 *h_eff_2nu = new TH1F("eff_2nu","eff_2nu",100,0,5);
  TH1 *h_eff_tl208 = new TH1F("eff_tl208","eff_tl208",100,0,5);
  TH1 *h_eff_bi214 = new TH1F("eff_bi214","eff_bi214",100,0,5);
  TH1 *h_eff_radon = new TH1F("eff_radon","eff_radon",100,0,5);


  for(unsigned int i = 1; i <= 100; ++i) {
    h_eff_0nu->Fill(h_0nu->GetBinCenter(i),eff_0nu*h_0nu->Integral(i,100));
    h_eff_2nu->Fill(h_2nu->GetBinCenter(i),eff_2nu*h_2nu->Integral(i,100));
    h_eff_tl208->Fill(h_tl208->GetBinCenter(i),eff_tl208*h_tl208->Integral(i,100));
    h_eff_bi214->Fill(h_bi214->GetBinCenter(i),eff_bi214*h_bi214->Integral(i,100));
    h_eff_radon->Fill(h_radon->GetBinCenter(i),eff_radon*h_radon->Integral(i,100));
  }

  h_eff_0nu->SetLineColor(kRed);
  h_eff_2nu->SetLineColor(kBlue);
  h_eff_tl208->SetLineColor(kGreen+1);
  h_eff_bi214->SetLineColor(kOrange);
  h_eff_radon->SetLineColor(kMagenta);

  h_eff_0nu->Write();
  h_eff_2nu->Write();
  h_eff_tl208->Write();
  h_eff_bi214->Write();
  h_eff_radon->Write();

  h_0nu->Scale(eff_0nu);
  h_2nu->Scale(eff_2nu * mass * exposure_y * Na * log(2) / M_Se / halflife_2nu);
  h_tl208->Scale(eff_tl208 * 2e-6 * mass * exposure_sec);
  h_bi214->Scale(eff_bi214 * 10e-6 * mass * exposure_sec);
  h_radon->Scale(eff_radon * 10e-6 * tracker_volume * exposure_sec);

  std::cout << " Counts : " << std::endl
            << " 2nu   : " << eff_2nu * mass * exposure_y * Na * log(2) / M_Se / halflife_2nu << std::endl
            << " tl208 : " << eff_tl208 * 2e-6 * mass * exposure_sec << std::endl
            << " bi214 : " << eff_bi214 * 10e-6 * mass * exposure_sec << std::endl
            << " radon : " << eff_radon * 150e-6 * tracker_volume * exposure_sec << std::endl;

  TGraph *g_hl = new TGraph();
  TGraph *g_hl_radon = new TGraph();
  TGraph *g_hl_radon_rel = new TGraph();
  TGraph *g_count_radon = new TGraph();
  TGraph *g_count_radon_window = new TGraph();

  double max_hl_no_radon = 1;
  for(unsigned int i_radon=0; i_radon<=50;++i_radon) {
    double max_hl = 0;
    double max_hl_radon_count = 0;
    for(unsigned int i = 40; i <= 64; ++i) {
      double n_bkg = h_2nu->Integral(i,64) + h_tl208->Integral(i,64) + h_bi214->Integral(i,64) + h_radon->Integral(i,64)*i_radon;
      double halflife = h_0nu->Integral(i,64) * mass * exposure_y * Na * log(2) / M_Se / get_number_of_excluded_events(n_bkg);
      if(halflife>max_hl) {
        max_hl=halflife;
        max_hl_radon_count = h_radon->Integral(i,64)*i_radon;
        if(i_radon==0)
          max_hl_no_radon = halflife;
      }
      // g_hl->SetPoint(i-40,h_2nu->GetBinLowEdge(i),halflife);
    }
    g_hl_radon->SetPoint(i_radon,i_radon*10,max_hl/1e24);
    g_hl_radon_rel->SetPoint(i_radon,i_radon*10,max_hl/max_hl_no_radon);
    g_count_radon->SetPoint(i_radon,i_radon*10,max_hl_radon_count);
    g_count_radon_window->SetPoint(i_radon,i_radon*10,i_radon*h_radon->Integral(57,64));
  }

  TGraph *g_hl_radon_exp_50 = new TGraph();
  TGraph *g_hl_radon_exp_150 = new TGraph();
  TGraph *g_hl_radon_exp_500 = new TGraph();

  vector <double> exposure_factor;
  exposure_factor.push_back(1/17.5);
  exposure_factor.push_back(5/17.5);
  exposure_factor.push_back(1);
  exposure_factor.push_back(50/17.5);
  exposure_factor.push_back(100/17.5);
  exposure_factor.push_back(150/17.5);
  exposure_factor.push_back(250/17.5);
  exposure_factor.push_back(350/17.5);
  exposure_factor.push_back(500/17.5);

  //50 muBq
  for(unsigned int i_exp=0;i_exp<exposure_factor.size();++i_exp) {
    double max_hl = 0;
    for(unsigned int i = 40; i <= 64; ++i) {
      double n_bkg = h_2nu->Integral(i,64) + h_tl208->Integral(i,64) + h_bi214->Integral(i,64) + h_radon->Integral(i,64)*5;
      n_bkg = n_bkg * exposure_factor[i_exp];
      double halflife = h_0nu->Integral(i,64) * mass * exposure_y * exposure_factor[i_exp] * Na * log(2) / M_Se / get_number_of_excluded_events(n_bkg);
      if(halflife>max_hl)
        max_hl=halflife;
    }
    g_hl_radon_exp_50->SetPoint(i_exp,17.5*exposure_factor[i_exp],max_hl/1e25);
  }

  //150 muBq
  for(unsigned int i_exp=0;i_exp<exposure_factor.size();++i_exp) {
    double max_hl = 0;
    for(unsigned int i = 40; i <= 64; ++i) {
      double n_bkg = h_2nu->Integral(i,64) + h_tl208->Integral(i,64) + h_bi214->Integral(i,64) + h_radon->Integral(i,64)*15;
      n_bkg = n_bkg * exposure_factor[i_exp];
      double halflife = h_0nu->Integral(i,64) * mass * exposure_y * exposure_factor[i_exp] * Na * log(2) / M_Se / get_number_of_excluded_events(n_bkg);
      if(halflife>max_hl)
        max_hl=halflife;
    }
    g_hl_radon_exp_150->SetPoint(i_exp,17.5*exposure_factor[i_exp],max_hl/1e25);
  }

  //500 muBq
  for(unsigned int i_exp=0;i_exp<exposure_factor.size();++i_exp) {
    double max_hl = 0;
    for(unsigned int i = 40; i <= 64; ++i) {
      double n_bkg = h_2nu->Integral(i,64) + h_tl208->Integral(i,64) + h_bi214->Integral(i,64) + h_radon->Integral(i,64)*50;
      n_bkg = n_bkg * exposure_factor[i_exp];
      double halflife = h_0nu->Integral(i,64) * mass * exposure_y * exposure_factor[i_exp] * Na * log(2) / M_Se / get_number_of_excluded_events(n_bkg);
      if(halflife>max_hl)
        max_hl=halflife;
    }
    g_hl_radon_exp_500->SetPoint(i_exp,17.5*exposure_factor[i_exp],max_hl/1e25);
  }

  g_hl_radon->SetTitle(";Radon contamination [#muBq/m^{3}];T_{1/2}^{0#nu} [10^{24} y]");
  g_hl_radon->SetLineColor(kBlue+1);
  g_hl_radon->SetLineWidth(2);

  g_hl_radon_rel->SetTitle(";Radon contamination [#muBq/m^{3}];Relative sensitivity loss [%]");
  g_hl_radon_rel->SetLineColor(kBlue+1);
  g_hl_radon_rel->SetLineWidth(2);

  g_count_radon->SetTitle(";Radon contamination [#muBq/m^{3}];Radon events in the ROI");
  g_count_radon->SetLineColor(kBlue+1);
  g_count_radon->SetLineWidth(2);

  g_count_radon_window->SetTitle(";Radon contamination [#muBq/m^{3}];Radon events in the [2.8,3.2]MeV window");
  g_count_radon_window->SetLineColor(kBlue+1);
  g_count_radon_window->SetLineWidth(2);

  h_0nu->SetName("0nu");
  h_2nu->SetName("2nu");
  h_tl208->SetName("tl208");
  h_bi214->SetName("bi214");
  h_radon->SetName("radon");

  g_hl_radon_exp_50->SetTitle(";Exposure [kg.y];T_{1/2}^{0#nu} [10^{25} y]");
  g_hl_radon_exp_50->SetLineColor(kGreen+1);
  g_hl_radon_exp_50->SetLineWidth(2);

  g_hl_radon_exp_150->SetTitle(";Exposure [kg.y];T_{1/2}^{0#nu} [10^{25} y]");
  g_hl_radon_exp_150->SetLineColor(kBlue);
  g_hl_radon_exp_150->SetLineWidth(2);

  g_hl_radon_exp_500->SetTitle(";Exposure [kg.y];T_{1/2}^{0#nu} [10^{25} y]");
  g_hl_radon_exp_500->SetLineColor(kRed);
  g_hl_radon_exp_500->SetLineWidth(2);

  g_hl->SetName("hl");
  g_hl->Write();

  g_hl_radon->SetName("hl_radon");
  g_hl_radon->Write();

  g_hl_radon_rel->SetName("hl_radon_rel");
  g_hl_radon_rel->Write();

  g_count_radon->SetName("radon_count");
  g_count_radon->Write();

  g_count_radon_window->SetName("radon_count_window");
  g_count_radon_window->Write();

  g_hl_radon_exp_50->SetName("hl_radon_exp_50");
  g_hl_radon_exp_50->Write();

  g_hl_radon_exp_150->SetName("hl_radon_exp_150");
  g_hl_radon_exp_150->Write();

  g_hl_radon_exp_500->SetName("hl_radon_exp_500");
  g_hl_radon_exp_500->Write();


  h_0nu->Write();
  h_2nu->Write();
  h_tl208->Write();
  h_bi214->Write();
  h_radon->Write();

  return;
}
