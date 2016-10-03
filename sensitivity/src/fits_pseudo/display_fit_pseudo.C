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

void display_fit_pseudo() {

  // TFile *f_fits = new TFile("../fits.root","RECREATE");
  TFile *f_fits = TFile::Open("./fits.root","RECREATE");

  // TString quantity = "1e_electron_energy";
  // TString quantity = "2e_electron_maximal_energy";
  // TString quantity = "2e_electrons_energy_sum";
  // // TString quantity = "2e_electrons_internal_probability";
  // TString quantity = "2e_electrons_cos_angle";

  // // TString quantity = "1e_electron_energy";
  TString quantity = "1e1g_gamma_energy";
  // TString quantity = "1e1g_electron_gamma_energy_sum";
  // TString quantity = "1e2g_gamma_max_energy";
  // TString quantity = "1e2g_electron_gammas_energy_sum";
  // TString quantity = "1e3g_electron_gammas_energy_sum";
  // // TString quantity = "2e1g_electrons_gammas_energy_sum";

  double se_hl_meas = 9e19;
  double tl_activity_meas = 2e-6;
  double bi_activity_meas = 10e-6;
  double radon_activity_meas = 150e-6;

  // double nu_quantity_efficiency = 0.210014; //1e
  // double tl208_quantity_efficiency = 0.0229171; //1e
  // double bi214_quantity_efficiency = 0.121589; //1e
  // double radon_quantity_efficiency = 0.0114585; //1e

  // double nu_quantity_efficiency = 0.10154; //2e
  // double tl208_quantity_efficiency = 0.00121521; //2e
  // double bi214_quantity_efficiency = 0.00165497; //2e
  // double radon_quantity_efficiency = 0.00024013; //2e

  double nu_quantity_efficiency = 0.00384706; //1e1g
  double tl208_quantity_efficiency = 0.137997; //1e1g
  double bi214_quantity_efficiency = 0.121921; //1e1g
  double radon_quantity_efficiency = 0.0112898; //1e1g

  // double nu_quantity_efficiency = 8.453e-05; //1e2g
  // double tl208_quantity_efficiency = 0.137908; //1e2g
  // double bi214_quantity_efficiency = 0.0635762; //1e2g
  // double radon_quantity_efficiency = 0.00576512; //1e2g

  // double nu_quantity_efficiency = 1.41e-06; //1e3g
  // double tl208_quantity_efficiency = 0.0350548; //1e3g
  // double bi214_quantity_efficiency = 0.00944775; //1e3g
  // double radon_quantity_efficiency = 0.00085126; //1e3g

  // double nu_quantity_efficiency = 0.00086855; //2e1g
  // double tl208_quantity_efficiency = 0.00389406; //2e1g
  // double bi214_quantity_efficiency = 0.00227784; //2e1g
  // double radon_quantity_efficiency = 0.0002913; //2e1g

  TFile * f_2nu = TFile::Open("../pdf/2nu_pdf.root");
  TH1F *nu_pdf = (TH1F*)f_2nu->Get(quantity);
  TFile * f_tl208 = TFile::Open("../pdf/tl208_pdf.root");
  TH1F *tl208_pdf = (TH1F*)f_tl208->Get(quantity);
  TFile * f_bi214 = TFile::Open("../pdf/bi214_pdf.root");
  TH1F *bi214_pdf = (TH1F*)f_bi214->Get(quantity);
  TFile * f_radon = TFile::Open("../pdf/radon_pdf.root");
  TH1F *radon_pdf = (TH1F*)f_radon->Get(quantity);

  nu_pdf->SetLineColor(kBlue);
  nu_pdf->SetFillColor(kBlue);
  nu_pdf->SetMarkerColor(kBlue);

  tl208_pdf->SetLineColor(kGreen+1);
  tl208_pdf->SetFillColor(kGreen+1);
  tl208_pdf->SetMarkerColor(kGreen+1);

  bi214_pdf->SetLineColor(kOrange-3);
  bi214_pdf->SetFillColor(kOrange-3);
  bi214_pdf->SetMarkerColor(kOrange-3);

  radon_pdf->SetLineColor(kMagenta);
  radon_pdf->SetFillColor(kMagenta);
  radon_pdf->SetMarkerColor(kMagenta);

  TFile * f_pseudo = TFile::Open("../pseudo/pseudo.root");
  TH1F *pseudo = (TH1F*)f_pseudo->Get(quantity);
  pseudo->SetLineColor(kBlack);
  // pseudo->SetMarkerStyle(1);

  THStack *hs = new THStack(quantity,quantity);

  nu_pdf->Scale(nu_quantity_efficiency*7*2.5*log(2)*6.022e23/se_hl_meas/0.082);
  tl208_pdf->Scale(tl_activity_meas*tl208_quantity_efficiency*7*2.5*3.14e7);
  bi214_pdf->Scale(bi_activity_meas*bi214_quantity_efficiency*7*2.5*3.14e7);
  radon_pdf->Scale(radon_activity_meas*radon_quantity_efficiency*15.3*2.5*3.14e7);

 std:cout << " n2nu = " << nu_quantity_efficiency*7*2.5*log(2)*6.022e23/se_hl_meas/0.082 << std::endl
          << " ntl = " << tl_activity_meas*tl208_quantity_efficiency*7*2.5*3.14e7 << std::endl
          << " nbi = " << bi_activity_meas*bi214_quantity_efficiency*7*2.5*3.14e7 << std::endl
          << " nradon = " << radon_activity_meas*radon_quantity_efficiency*15.3*2.5*3.14e7 << std::endl;

  hs->Add(tl208_pdf);
  hs->Add(bi214_pdf);
  hs->Add(radon_pdf);
  hs->Add(nu_pdf);

  TCanvas *c1 = new TCanvas("c1","example",600,700);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  // pad1->SetBottomMargin(0.05);
  pad1->SetTopMargin(0.03);
  pad1->Draw();
  pad1->cd();

  Double_t xl1=.75, yl1=0.65, xl2=0.9, yl2=0.97;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(nu_pdf,"2#nu");
  leg->AddEntry(tl208_pdf,"^{208}Tl");
  leg->AddEntry(bi214_pdf,"^{214}Bi");
  leg->AddEntry(radon_pdf,"Radon");
  // leg->AddEntry((TObject*)0,"#chi^{2}/ndf : 48.53/68", "");
  leg->AddEntry((TObject*)0,"#chi^{2}/ndf : 79.45/58", "");
  leg->SetFillColor(kWhite);

  // hs->SetTitle("1e1#gamma : Gamma energy; Gamma energy; Events");
  hs->SetTitle("1e1#gamma : Gamma energy;; Events");
  // hs->SetTitle("1e1#gamma : Energy sum;; Events");

  hs->DrawClone();
  leg->Draw("same");

  // hs->GetXaxis()->SetTitle("Gamma energy");
  // hs->GetXaxis()->SetTitleFont(62);
  // hs->GetXaxis()->SetTitleSize(0.05);
  // hs->GetXaxis()->SetTitleOffset(0.90);
  // hs->GetYaxis()->SetTitle("Events");
  // hs->GetYaxis()->SetTitleFont(62);
  // hs->GetYaxis()->SetTitleSize(0.05);
  // hs->GetYaxis()->SetTitleOffset(0.90);
  // hs->SetLineColor(kBlue);
  // hs->SetMarkerColor(kBlue);
  // hs->SetLineWidth(2);
  // hs->SetFillColor(kWhite);

  pseudo->DrawClone("samePE");

  c1->cd();

  ////// CHI2

  TH1F *h_tmp = new TH1F("h_tmp","h_tmp",100,0,5);

  double chi2 = 0;
  double nbin_data = 0;
  for (unsigned int i = 1; i <= 100; ++i) {
    h_tmp->Fill(i,nu_pdf->GetBinContent(i)+tl208_pdf->GetBinContent(i)+bi214_pdf->GetBinContent(i)+radon_pdf->GetBinContent(i));

    if(pseudo->GetBinError(i) != 0) {
      chi2 += pow((nu_pdf->GetBinContent(i)+tl208_pdf->GetBinContent(i)+bi214_pdf->GetBinContent(i)+radon_pdf->GetBinContent(i))-pseudo->GetBinContent(i),2)/pow(pseudo->GetBinError(i),2);
      ++nbin_data;
    }
    std::cout << "debug " <<  nu_pdf->GetBinContent(i)+tl208_pdf->GetBinContent(i)+bi214_pdf->GetBinContent(i)+radon_pdf->GetBinContent(i) << std::endl;
  }

  double chi;
  std::cout << "chi2 by hand " <<   chi2 << std::endl;
  std::cout << "   ndof " << nbin_data-1 << std::endl;
  // std::cout << "chi " <<   pseudo->Chi2TestX(h_tmp,&chi) << std::endl;
  std::cout << "     " <<   chi << std::endl;
  std::cout << "KS " <<   pseudo->KolmogorovTest(h_tmp) << std::endl;

  //////
  bi214_pdf->Sumw2();
  bi214_pdf->Add(tl208_pdf);
  bi214_pdf->Add(nu_pdf);
  bi214_pdf->Add(radon_pdf);
  //Error computation to check
  // bi214_pdf->Sumw2();

  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->Draw();
  pad2->cd();

  bi214_pdf->SetStats(0);
  pseudo->SetStats(0);

  // pseudo->Divide(bi214_pdf);

  for (unsigned int i =1;i<=pseudo->GetNbinsX();i++) {
    if(pseudo->GetBinContent(i)!=0) {
      double res = (pseudo->GetBinContent(i)- bi214_pdf->GetBinContent(i))/pseudo->GetBinError(i);
      pseudo->SetBinContent(i,res);
      pseudo->SetBinError(i,0);
      // std::cout << res << std::endl;
    }
  }
  pseudo->Sumw2();

  // pseudo->SetTitle("1e1#gamma : Gamma energy");
  pseudo->SetTitle("1e1#gamma : Energy sum");
  // pseudo->GetXaxis()->SetTitle("Gamma energy");
  pseudo->GetXaxis()->SetTitle("Energy sum");
  pseudo->GetXaxis()->SetTitleFont(62);
  pseudo->GetXaxis()->SetTitleSize(0.05);
  pseudo->GetXaxis()->SetTitleOffset(0.90);
  pseudo->GetYaxis()->SetTitle("Residuals");
  pseudo->GetYaxis()->SetTitleFont(62);
  pseudo->GetYaxis()->SetTitleSize(0.05);
  pseudo->GetYaxis()->SetTitleOffset(0.90);
  // pseudo->SetLineColor(kBlue);
  // pseudo->SetMarkerColor(kBlue);
  // pseudo->SetLineWidth(2);
  pseudo->SetFillColor(kWhite);

  pseudo->GetYaxis()->SetRangeUser(-5,5);
  pseudo->SetDrawOption("P");
  pseudo->Draw("P");

  TLine *line = new TLine(0,0,5,0);
  // line->SetLineColor(kBlack);
  line->SetLineStyle(kDashed);
  line->Draw("same");



  // c1->cd();
  // TString file = quantity + ".pdf";
  // c1->Print(file);

  return;
}
