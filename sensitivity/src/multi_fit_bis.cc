#include "TH1.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TSystem.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "stats_tools.h"
#include "sensitivity_constants.h"
// #include "analysis_config.h"
// #include "multi_fit.h"

void fcn_to_minimize(int& npar, double* deriv, double& f, double par[], int flag)
{
  std::map < std::string, double > isotope_activity;
  isotope_activity.insert(std::pair<std::string,double>("bi214",100e-6));
  isotope_activity.insert( std::pair<std::string,double>("tl208",100e-6) );

  std::vector <TString> quantities;
  quantities.push_back("1e1g_electron_gamma_energy_sum");
  // quantities.push_back("1e2g_electron_gammas_energy_sum");

  TFile * f_pseudo = TFile::Open("pseudo.root");

  f = 0.;

  for(std::vector<TString>::iterator j = quantities.begin(); j != quantities.end(); ++j) {
  // for(auto j = quantities.begin(); j != quantities.end(); ++j) {
    // std::cout << "quantity " << *j << std::endl;

    for(int i = 0; i < 100; i++)
      {
        double b = 0;
        unsigned int count = 0;
        double efficiency[2] = {0.123019,0.139367};
        for(std::map<std::string, double>::iterator k = isotope_activity.begin(); k != isotope_activity.end(); ++k) {
        // for(auto k = isotope_activity.begin(); k != isotope_activity.end(); ++k) {
          TString qty = *j;
          TString isotope = k->first;
          // std::cout << "isotope " << isotope << std::endl;

          TString channel = isotope + "_" + qty;
          // double efficiency = quantity_efficiency.at(channel);
          TString pdf_file = isotope + "_pdf.root";
          TFile *file = TFile::Open(pdf_file);

          TH1F *h = (TH1F*)file->Get(qty);

          b += par[count] * efficiency[count] * h->GetBinContent(i) * mass * exposure;
          // std::cout << "i activity  efficiency  histo_eff   " << std::endl
          //           << i << "  " << par[count] << "  " << efficiency << "  " << h->GetBinContent(i) << std::endl;
          count++;
          file->Close();
        }
        TH1F *pseudo = (TH1F*)f_pseudo->Get(*j);
        double d = pseudo->GetBinContent(i);

        if(b==0)
          b=1e-15;
        // std::cout << "i b d " << i << "  " << b << "  " << d << std::endl;
        f += 2*(b-d*log(b)+log_factorial(d));
      }
  }

  f_pseudo->Close();
  return;
}

// double likelihood(double activity)
// {
//   TFile * f_bi214 = TFile::Open("bi214_pdf.root");

//   TH1F *bi214_alpha_track_length_pdf = (TH1F*)f_bi214->Get("h");

//   TFile * f_bi214_pseudo = TFile::Open("pseudo.root");

//   TH1F *bi214_alpha_track_length_pseudo = (TH1F*)f_bi214_pseudo->Get("h");

//   double f = 0.;

//   for(int i = 0; i < 100; i++)
//     {
//       double b = activity*bi214_channel_1e1a_efficiency*bi214_alpha_track_length_pdf->GetBinContent(i)*7*3.14e7*2.5;
//       int d = bi214_alpha_track_length_pseudo->GetBinContent(i);
//       // std::cout << "bin " << i << "  :  b = " << b << "  d = " << d << std::endl;
//       // std::cout << "bin likelihood " << 2*(b-d*log(b+0.1)+log(factorial(d))) << std::endl;
//       f += 2*(b-d*log(b+0.1)+log_factorial(d));
//     }

//   // std::cout << " DEBUG " << activity*1e6 << "  " << f << std::endl;
//   return f;
// }

// void multi_fit_bis(std::map < std::string, std::vector<double> > & activity_measurement)
void multi_fit_bis()
{
  // gSystem->Load("/home/calvez/nemo/work_dir/Cadfael/build/ports/root/root-prefix/src/root/lib/libMinuit.so");

  const int npar = 2;

  TMinuit minuit(npar);

  minuit.SetFCN(fcn_to_minimize);  //set the fcn to minimize

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes
  double minVal[npar];            // minimum bound on parameter
  double maxVal[npar];            // maximum bound on parameter
  std::string parName[npar];           // parameter name

   par[0] = 100e-6;
   stepSize[0] = 1e-7;
   minVal[0] = 1e-9;
   maxVal[0] = 1e-3;
   parName[0] = "activity in 214Bi uBq/kg";

   minuit.DefineParameter(0, parName[0].c_str(),
                            par[0], stepSize[0], minVal[0], maxVal[0]);

   par[1] = 20e-6;
   stepSize[1] = 1e-7;
   minVal[1] = 1e-9;
   maxVal[1] = 1e-3;
   parName[1] = "activity in 208Tl uBq/kg";

   minuit.DefineParameter(1, parName[1].c_str(),
                            par[1], stepSize[1], minVal[1], maxVal[1]);

  // unsigned int count = 0;
  // for(auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
  //   par[count] = 20./1e6;
  //   stepSize[count] = 1e-7;
  //   minVal[count] = 1e-9;
  //   maxVal[count] = 1e-3;
  //   parName[count] = "activity in " + i->first + "in uBq/kg";
  //   count++;

  //   minuit.DefineParameter(count, parName[count].c_str(),
  //                          par[count], stepSize[count], minVal[count], maxVal[count]);
  // }

  minuit.Migrad();

  // unsigned int count_bis = 0;
  // for(auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
  //   double activity = 0;
  //   double activity_err = 0;
  //   minuit.GetParameter(count,activity,activity_err);
  //   activity_measurement.insert(std::pair<std::string,std::vector<double>>(i->first,{activity,activity_err}));
  //   count_bis++;
  // }


  double activity_bi214 = 0;
  double activity_err_bi214 = 0;
  minuit.GetParameter(0,activity_bi214,activity_err_bi214);
  std::cout << "activity_bi214 " << activity_bi214 << " +/- " << activity_err_bi214 << std::endl;

  double activity_tl208 = 0;
  double activity_err_tl208 = 0;
  minuit.GetParameter(1,activity_tl208,activity_err_tl208);
  std::cout << "activity_tl208 " << activity_tl208 << " +/- " << activity_err_tl208 << std::endl;


  // std::cout << " Measured bi214 activity is : " << activity_bi214*1e6 << " +/- " << activity_bi214_err*1e6 << " uBq/kg" << std::endl;
  // std::cout << " Measured tl208 activity is : " << activity_tl208*1e6 << " +/- " << activity_tl208_err*1e6 << " uBq/kg" << std::endl;

  // TGraph *g_likelihood = new TGraph(100);
  // double activity_start = 1e-5;
  // double activity_end = 1e-3;

  // for (unsigned int i = 0; i<100;++i) {
  //   double activity = activity_start + i*(activity_end-activity_start)/100.;
  //   double Q = likelihood(activity);
  //   g_likelihood->SetPoint(i,activity*1e6,Q);
  // }
  // g_likelihood->Draw();


  //tmp
  // double efficiency[2] = {0.139367,0.123019};
  double tl208_channel_1e1g_efficiency = 0.139367;
  double bi214_channel_1e1g_efficiency = 0.123019;

  TFile * f_bi214 = TFile::Open("bi214_pdf.root");
  TH1F *bi214_pdf = (TH1F*)f_bi214->Get("1e1g_electron_gamma_energy_sum");

  bi214_pdf->SetLineColor(kOrange);
  bi214_pdf->SetFillColor(kOrange);

  TFile * f_tl208 = TFile::Open("tl208_pdf.root");
  TH1F *tl208_pdf = (TH1F*)f_tl208->Get("1e1g_electron_gamma_energy_sum");
  tl208_pdf->SetLineColor(kGreen+1);
  tl208_pdf->SetFillColor(kGreen+1);

  TFile * f_pseudo = TFile::Open("pseudo.root");
  TH1F *pseudo = (TH1F*)f_pseudo->Get("1e1g_electron_gamma_energy_sum");
  pseudo->SetLineColor(kBlack);
  // pseudo->SetMarkerStyle(1);

  THStack *hs = new THStack("hs","Energy [keV]");
  bi214_pdf->Scale(activity_bi214*bi214_channel_1e1g_efficiency*mass*exposure);
  tl208_pdf->Scale(activity_tl208*tl208_channel_1e1g_efficiency*mass*exposure);

  hs->Add(bi214_pdf);
  hs->Add(tl208_pdf);

  TCanvas *c1 = new TCanvas("c1","example",600,700);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  // pad1->SetBottomMargin(0.05);
  pad1->SetTopMargin(0.03);
  pad1->Draw();
  pad1->cd();

  hs->DrawClone();

  pseudo->DrawClone("samePE");
  c1->cd();

  // TH1 *h_sum = new TH1F("h_sum","h_sum");

  bi214_pdf->Sumw2();
  bi214_pdf->Add(tl208_pdf);
  //Error computation to check
  // bi214_pdf->Sumw2();
  pseudo->Sumw2();

  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0.0);
  // pad2->SetBottomMargin(0.05);
  pad2->Draw();
  pad2->cd();

  bi214_pdf->SetStats(0);
  pseudo->SetStats(0);

  pseudo->Divide(bi214_pdf);

  pseudo->GetYaxis()->SetRangeUser(0.6,1.4);
  pseudo->Draw("ep");

  TLine *line = new TLine(0,1,5,1);
  // line->SetLineColor(kBlack);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  // c1->cd();

  return;
}
