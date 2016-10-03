#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPad.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <iostream>

#include "stats_tools.h"
#include "sensitivity_constants.h"
// #include "sensitivity_measurements.h"
#include "analysis_config.h"
#include "multi_fit.h"

extern std::map < TString , double > quantity_efficiency;
extern std::map < TString , TH1F* > quantity_pdf;

void fcn_to_minimize(int& /*npar*/, double* /*deriv*/, double& f, double par[], int /*flag*/)
{
  TFile * f_pseudo = TFile::Open("../pseudo/pseudo.root");
  f = 0.;
  for(auto j = quantities.begin(); j != quantities.end(); ++j) {
    TString qty = *j;
    for(int i = 0; i < 100; i++)
      {
        double b = 0;
        unsigned int count = 0;
        for(auto k = isotope_activity.begin(); k != isotope_activity.end(); ++k) {
          TString isotope = k->first;
          TString channel = isotope + "_" + qty;
          double efficiency = quantity_efficiency.at(channel);
          TString pdf_file = "../pdf/" + isotope + "_pdf.root";
          TFile * f = TFile::Open(pdf_file);
          TH1F *h = (TH1F*)f->Get(qty);

          if(isotope.Contains("2nu"))
            b += const_se / (par[count] * 1e19) * efficiency * h->GetBinContent(i) * mass * exposure_y;
          else if(isotope.Contains("radon"))
            b += par[count] * efficiency * h->GetBinContent(i) * tracker_volume * exposure_sec;
          else
            b += par[count] * efficiency * h->GetBinContent(i) * mass * exposure_sec;

          count++;
          f->Close();
        }

        TH1F *pseudo = (TH1F*)f_pseudo->Get(*j);
        // std::cout << " getting pseudo experiment " << std::endl;
        // std::cout << " pseudo is " << quantity_pdf.size() << std::endl;
        // double d = quantity_pdf.at(pseudo_qty)->GetBinContent(i);

        // TString pseudo_qty = "pseudo_" + qty;
        // // TH1F *pseudo = quantity_pdf.at(pseudo_qty);
        // std::cout << "   test h   "  << pseudo_qty << std::endl;
        // TH1F *pseudo = quantity_pdf.at(pseudo_qty);

        // std::cout << "            " << pseudo->GetNbinsX() << std::endl;

        // std::cout << "            "  << quantity_pdf.at(pseudo_qty)->GetBinContent(20) << std::endl;

        // for(auto it = quantity_pdf.begin(); it != quantity_pdf.end(); ++it) {
        //   std::cout << " test "  << it->first << std::endl;
        //   // std::cout << "      "  << it->second->GetBinContent(1) << std::endl;
        // }

        // auto it = quantity_pdf.begin();
        // it++;
        // it++;
        // std::cout << "test first " << it->first() << std::end;
        // std::cout << "test second " << it->second->Integral(0,100) << std::end;
        // std::cout << "      " << pseudo->Integral(0,100) << std::endl;

        // double d = quantity_pdf.at(pseudo_qty)->GetBinContent(i);
        double d = pseudo->GetBinContent(i);
        // double d = 1;

        if(b==0)
          b=1e-15;
        // std::cout << "i b d " << i << "  " << b << "  " << d << std::endl;
        f += 2*(b-d*log(b)+log_factorial(d));
      }
  }
  // f_pseudo->Close();

  std::cout << " New minimization  " << f << "    ";
  // << "  par0 : " << par[0] << "     par1 : " << par[1] << "     par2 : " << par[2] << "     par3 " << par[3] << std::endl;

  unsigned int c = 0;
  for(auto  i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    std::cout << i->first << " : " << par[c] << "    ";
    ++c;
  }
  std::cout << std::endl;
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

void multi_fit(std::map < std::string, std::vector<double> > & activity_measurement)
{
  std::cout << " Starting the fit " << std::endl;

  const int npar = isotope_activity.size();

  TMinuit minuit(npar);

  minuit.SetFCN(fcn_to_minimize);  //set the fcn to minimize

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes
  double minVal[npar];            // minimum bound on parameter
  double maxVal[npar];            // maximum bound on parameter
  std::string parName[npar];           // parameter name

  unsigned int count = 0;
  for(auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    const double param_guess = i->second;
    par[count] = param_guess;
    stepSize[count] = param_guess/10.;
    minVal[count] = param_guess/10.;
    maxVal[count] = param_guess * 10.;
    parName[count] = i->first;
    minuit.DefineParameter(count, parName[count].c_str(),
                           par[count], stepSize[count], minVal[count], maxVal[count]);
    count++;
  }

  minuit.FixParameter(0); //2nu
  minuit.FixParameter(2); //radon
  // minuit.FixParameter(1);
  // minuit.FixParameter(3);

  minuit.Migrad();
  // minuit.mnsimp(); //shit, not converging
  // minuit.mnseek(); //shit

  unsigned int count_bis = 0;
  for(auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    double activity = 0;
    double activity_err = 0;
    minuit.GetParameter(count_bis,activity,activity_err);
    activity_measurement.insert(std::pair<std::string,std::vector<double>>(i->first,{activity,activity_err}));
    count_bis++;
  }

  // TGraph *g_likelihood = new TGraph(100);
  // double activity_start = 1e-5;
  // double activity_end = 1e-3;

  // for (unsigned int i = 0; i<100;++i) {
  //   double activity = activity_start + i*(activity_end-activity_start)/100.;
  //   double Q = likelihood(activity);
  //   g_likelihood->SetPoint(i,activity*1e6,Q);
  // }
  // g_likelihood->Draw();

  //-----
  /*
  if(print_fits && number_of_pseudo_experiments==1) {
    // TFile *f_fits = new TFile("../fits.root","RECREATE");
    TFile *f_fits = TFile::Open("../fits.root","RECREATE");

    for(auto i = quantities.begin(); i != quantities.end(); ++i)
      {
        TString quantity = *i;
        TString tl208_quantity = "tl208_" + quantity;

        std::cout << "test " << quantity_efficiency.size() << "  " <<  tl208_quantity << std::endl;

        double tl208_quantity_efficiency = quantity_efficiency.at(tl208_quantity);

        TString bi214_quantity = "bi214_" + quantity;
        double bi214_quantity_efficiency = quantity_efficiency.at(bi214_quantity);

        TFile * f_bi214 = TFile::Open("../bi214_pdf.root");
        TH1F *bi214_pdf = (TH1F*)f_bi214->Get(quantity);

        bi214_pdf->SetLineColor(kOrange);
        bi214_pdf->SetFillColor(kOrange);

        TFile * f_tl208 = TFile::Open("../tl208_pdf.root");
        TH1F *tl208_pdf = (TH1F*)f_tl208->Get(quantity);

        tl208_pdf->SetLineColor(kGreen+1);
        tl208_pdf->SetFillColor(kGreen+1);

        TFile * f_pseudo = TFile::Open("../pseudo.root");
        TH1F *pseudo = (TH1F*)f_pseudo->Get(quantity);
        pseudo->SetLineColor(kBlack);
        // pseudo->SetMarkerStyle(1);

        THStack *hs = new THStack(quantity,"Variable");
        bi214_pdf->Scale(isotope_activity.at("bi214")*bi214_quantity_efficiency*mass*exposure_sec);
        tl208_pdf->Scale(isotope_activity.at("tl208")*tl208_quantity_efficiency*mass*exposure_sec);

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

        bi214_pdf->Sumw2();
        bi214_pdf->Add(tl208_pdf);
        //Error computation to check
        // bi214_pdf->Sumw2();
        pseudo->Sumw2();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0.0);
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

        c1->cd();
        TString file = quantity + ".pdf";
        c1->Print(file);
      }
  }*/

  return;
}
