#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void edit_pdf_distrib()
{

  TFile f_2nu("./2nu_pdf.root");
  TFile f_tl("./tl208_pdf.root");
  TFile f_bi("./bi214_pdf.root");
  TFile f_radon("./radon_pdf.root");

  TH1F *h_2nu_1e = (TH1F*)f_2nu.Get("1e_electron_energy");
  TH1F *h_2nu_1e1g_egamma = (TH1F*)f_2nu.Get("1e1g_gamma_energy");
  TH1F *h_2nu_1e1g_esum = (TH1F*)f_2nu.Get("1e1g_electron_gamma_energy_sum");
  TH1F *h_2nu_1e2g_egamma_max = (TH1F*)f_2nu.Get("1e2g_gamma_max_energy");
  TH1F *h_2nu_1e2g_esum = (TH1F*)f_2nu.Get("1e2g_electron_gammas_energy_sum");
  TH1F *h_2nu_1e3g_esum = (TH1F*)f_2nu.Get("1e3g_electron_gammas_energy_sum");
  TH1F *h_2nu_2e1g_esum = (TH1F*)f_2nu.Get("2e1g_electrons_gammas_energy_sum");

  TH1F *h_tl_1e = (TH1F*)f_tl.Get("1e_electron_energy");
  TH1F *h_tl_1e1g_egamma = (TH1F*)f_tl.Get("1e1g_gamma_energy");
  TH1F *h_tl_1e1g_esum = (TH1F*)f_tl.Get("1e1g_electron_gamma_energy_sum");
  TH1F *h_tl_1e2g_egamma_max = (TH1F*)f_tl.Get("1e2g_gamma_max_energy");
  TH1F *h_tl_1e2g_esum = (TH1F*)f_tl.Get("1e2g_electron_gammas_energy_sum");
  TH1F *h_tl_1e3g_esum = (TH1F*)f_tl.Get("1e3g_electron_gammas_energy_sum");
  TH1F *h_tl_2e1g_esum = (TH1F*)f_tl.Get("2e1g_electrons_gammas_energy_sum");

  TH1F *h_bi_1e = (TH1F*)f_bi.Get("1e_electron_energy");
  TH1F *h_bi_1e1g_egamma = (TH1F*)f_bi.Get("1e1g_gamma_energy");
  TH1F *h_bi_1e1g_esum = (TH1F*)f_bi.Get("1e1g_electron_gamma_energy_sum");
  TH1F *h_bi_1e2g_egamma_max = (TH1F*)f_bi.Get("1e2g_gamma_max_energy");
  TH1F *h_bi_1e2g_esum = (TH1F*)f_bi.Get("1e2g_electron_gammas_energy_sum");
  TH1F *h_bi_1e3g_esum = (TH1F*)f_bi.Get("1e3g_electron_gammas_energy_sum");
  TH1F *h_bi_2e1g_esum = (TH1F*)f_bi.Get("2e1g_electrons_gammas_energy_sum");

  TH1F *h_radon_1e = (TH1F*)f_radon.Get("1e_electron_energy");
  TH1F *h_radon_1e1g_egamma = (TH1F*)f_radon.Get("1e1g_gamma_energy");
  TH1F *h_radon_1e1g_esum = (TH1F*)f_radon.Get("1e1g_electron_gamma_energy_sum");
  TH1F *h_radon_1e2g_egamma_max = (TH1F*)f_radon.Get("1e2g_gamma_max_energy");
  TH1F *h_radon_1e2g_esum = (TH1F*)f_radon.Get("1e2g_electron_gammas_energy_sum");
  TH1F *h_radon_1e3g_esum = (TH1F*)f_radon.Get("1e3g_electron_gammas_energy_sum");
  TH1F *h_radon_2e1g_esum = (TH1F*)f_radon.Get("2e1g_electrons_gammas_energy_sum");

  h_2nu_1e->SetStats(0);
  h_2nu_1e1g_egamma->SetStats(0);
  h_2nu_1e1g_esum->SetStats(0);
  h_2nu_1e2g_egamma_max->SetStats(0);
  h_2nu_1e2g_esum->SetStats(0);
  h_2nu_1e3g_esum->SetStats(0);
  h_2nu_2e1g_esum->SetStats(0);

  h_2nu_1e->SetTitle("1e : Electron energy");
  h_2nu_1e->GetXaxis()->SetTitle("Electron energy");
  h_2nu_1e->GetXaxis()->SetTitleFont(62);
  h_2nu_1e->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e->GetYaxis()->SetTitle("Probability");
  h_2nu_1e->GetYaxis()->SetTitleFont(62);
  h_2nu_1e->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e->SetLineColor(kBlue);
  h_2nu_1e->SetMarkerColor(kBlue);
  h_2nu_1e->SetLineWidth(2);
  h_2nu_1e->SetFillColor(kWhite);

  h_2nu_1e1g_egamma->SetTitle("1e1g : Gamma energy");
  h_2nu_1e1g_egamma->GetXaxis()->SetTitle("Gamma energy");
  h_2nu_1e1g_egamma->GetXaxis()->SetTitleFont(62);
  h_2nu_1e1g_egamma->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e1g_egamma->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e1g_egamma->GetYaxis()->SetTitle("Probability");
  h_2nu_1e1g_egamma->GetYaxis()->SetTitleFont(62);
  h_2nu_1e1g_egamma->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e1g_egamma->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e1g_egamma->SetLineColor(kBlue);
  h_2nu_1e1g_egamma->SetMarkerColor(kBlue);
  h_2nu_1e1g_egamma->SetLineWidth(2);
  h_2nu_1e1g_egamma->SetFillColor(kWhite);

  h_2nu_1e1g_esum->SetTitle("1e1g : Energy sum");
  h_2nu_1e1g_esum->GetXaxis()->SetTitle("Energy sum");
  h_2nu_1e1g_esum->GetXaxis()->SetTitleFont(62);
  h_2nu_1e1g_esum->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e1g_esum->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e1g_esum->GetYaxis()->SetTitle("Probability");
  h_2nu_1e1g_esum->GetYaxis()->SetTitleFont(62);
  h_2nu_1e1g_esum->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e1g_esum->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e1g_esum->SetLineColor(kBlue);
  h_2nu_1e1g_esum->SetMarkerColor(kBlue);
  h_2nu_1e1g_esum->SetLineWidth(2);
  h_2nu_1e1g_esum->SetFillColor(kWhite);

  h_2nu_1e2g_egamma_max->SetTitle("1e2g : Max gamma energy");
  h_2nu_1e2g_egamma_max->GetXaxis()->SetTitle("Max gamma energy");
  h_2nu_1e2g_egamma_max->GetXaxis()->SetTitleFont(62);
  h_2nu_1e2g_egamma_max->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e2g_egamma_max->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e2g_egamma_max->GetYaxis()->SetTitle("Probability");
  h_2nu_1e2g_egamma_max->GetYaxis()->SetTitleFont(62);
  h_2nu_1e2g_egamma_max->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e2g_egamma_max->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e2g_egamma_max->SetLineColor(kBlue);
  h_2nu_1e2g_egamma_max->SetMarkerColor(kBlue);
  h_2nu_1e2g_egamma_max->SetLineWidth(2);
  h_2nu_1e2g_egamma_max->SetFillColor(kWhite);

  h_2nu_1e2g_esum->SetTitle("1e2g : Energy sum");
  h_2nu_1e2g_esum->GetXaxis()->SetTitle("Energy sum");
  h_2nu_1e2g_esum->GetXaxis()->SetTitleFont(62);
  h_2nu_1e2g_esum->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e2g_esum->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e2g_esum->GetYaxis()->SetTitle("Probability");
  h_2nu_1e2g_esum->GetYaxis()->SetTitleFont(62);
  h_2nu_1e2g_esum->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e2g_esum->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e2g_esum->SetLineColor(kBlue);
  h_2nu_1e2g_esum->SetMarkerColor(kBlue);
  h_2nu_1e2g_esum->SetLineWidth(2);
  h_2nu_1e2g_esum->SetFillColor(kWhite);

  h_2nu_1e3g_esum->SetTitle("1e3g : Energy sum");
  h_2nu_1e3g_esum->GetXaxis()->SetTitle("Energy sum");
  h_2nu_1e3g_esum->GetXaxis()->SetTitleFont(62);
  h_2nu_1e3g_esum->GetXaxis()->SetTitleSize(0.05);
  h_2nu_1e3g_esum->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_1e3g_esum->GetYaxis()->SetTitle("Probability");
  h_2nu_1e3g_esum->GetYaxis()->SetTitleFont(62);
  h_2nu_1e3g_esum->GetYaxis()->SetTitleSize(0.05);
  h_2nu_1e3g_esum->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_1e3g_esum->SetLineColor(kBlue);
  h_2nu_1e3g_esum->SetMarkerColor(kBlue);
  h_2nu_1e3g_esum->SetLineWidth(2);
  h_2nu_1e3g_esum->SetFillColor(kWhite);

  h_2nu_2e1g_esum->SetTitle("2e1g : Energy sum");
  h_2nu_2e1g_esum->GetXaxis()->SetTitle("Energy sum");
  h_2nu_2e1g_esum->GetXaxis()->SetTitleFont(62);
  h_2nu_2e1g_esum->GetXaxis()->SetTitleSize(0.05);
  h_2nu_2e1g_esum->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_2e1g_esum->GetYaxis()->SetTitle("Probability");
  h_2nu_2e1g_esum->GetYaxis()->SetTitleFont(62);
  h_2nu_2e1g_esum->GetYaxis()->SetTitleSize(0.05);
  h_2nu_2e1g_esum->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_2e1g_esum->SetLineColor(kBlue);
  h_2nu_2e1g_esum->SetMarkerColor(kBlue);
  h_2nu_2e1g_esum->SetLineWidth(2);
  h_2nu_2e1g_esum->SetFillColor(kWhite);

  h_tl_1e->SetLineColor(kGreen+1);
  h_tl_1e->SetMarkerColor(kGreen+1);
  h_tl_1e->SetLineWidth(2);
  h_tl_1e->SetFillColor(kWhite);
  h_tl_1e1g_egamma->SetLineColor(kGreen+1);
  h_tl_1e1g_egamma->SetMarkerColor(kGreen+1);
  h_tl_1e1g_egamma->SetLineWidth(2);
  h_tl_1e1g_egamma->SetFillColor(kWhite);
  h_tl_1e1g_esum->SetLineColor(kGreen+1);
  h_tl_1e1g_esum->SetMarkerColor(kGreen+1);
  h_tl_1e1g_esum->SetLineWidth(2);
  h_tl_1e1g_esum->SetFillColor(kWhite);
  h_tl_1e2g_egamma_max->SetLineColor(kGreen+1);
  h_tl_1e2g_egamma_max->SetMarkerColor(kGreen+1);
  h_tl_1e2g_egamma_max->SetLineWidth(2);
  h_tl_1e2g_egamma_max->SetFillColor(kWhite);
  h_tl_1e2g_esum->SetLineColor(kGreen+1);
  h_tl_1e2g_esum->SetMarkerColor(kGreen+1);
  h_tl_1e2g_esum->SetLineWidth(2);
  h_tl_1e2g_esum->SetFillColor(kWhite);
  h_tl_1e3g_esum->SetLineColor(kGreen+1);
  h_tl_1e3g_esum->SetMarkerColor(kGreen+1);
  h_tl_1e3g_esum->SetLineWidth(2);
  h_tl_1e3g_esum->SetFillColor(kWhite);
  h_tl_2e1g_esum->SetLineColor(kGreen+1);
  h_tl_2e1g_esum->SetMarkerColor(kGreen+1);
  h_tl_2e1g_esum->SetLineWidth(2);
  h_tl_2e1g_esum->SetFillColor(kWhite);

  h_bi_1e->SetLineColor(kOrange-3);
  h_bi_1e->SetMarkerColor(kOrange-3);
  h_bi_1e->SetLineWidth(2);
  h_bi_1e->SetFillColor(kWhite);
  h_bi_1e1g_egamma->SetLineColor(kOrange-3);
  h_bi_1e1g_egamma->SetMarkerColor(kOrange-3);
  h_bi_1e1g_egamma->SetLineWidth(2);
  h_bi_1e1g_egamma->SetFillColor(kWhite);
  h_bi_1e1g_esum->SetLineColor(kOrange-3);
  h_bi_1e1g_esum->SetMarkerColor(kOrange-3);
  h_bi_1e1g_esum->SetLineWidth(2);
  h_bi_1e1g_esum->SetFillColor(kWhite);
  h_bi_1e2g_egamma_max->SetLineColor(kOrange-3);
  h_bi_1e2g_egamma_max->SetMarkerColor(kOrange-3);
  h_bi_1e2g_egamma_max->SetLineWidth(2);
  h_bi_1e2g_egamma_max->SetFillColor(kWhite);
  h_bi_1e2g_esum->SetLineColor(kOrange-3);
  h_bi_1e2g_esum->SetMarkerColor(kOrange-3);
  h_bi_1e2g_esum->SetLineWidth(2);
  h_bi_1e2g_esum->SetFillColor(kWhite);
  h_bi_1e3g_esum->SetLineColor(kOrange-3);
  h_bi_1e3g_esum->SetMarkerColor(kOrange-3);
  h_bi_1e3g_esum->SetLineWidth(2);
  h_bi_1e3g_esum->SetFillColor(kWhite);
  h_bi_2e1g_esum->SetLineColor(kOrange-3);
  h_bi_2e1g_esum->SetMarkerColor(kOrange-3);
  h_bi_2e1g_esum->SetLineWidth(2);
  h_bi_2e1g_esum->SetFillColor(kWhite);

  h_radon_1e->SetLineColor(kMagenta);
  h_radon_1e->SetMarkerColor(kMagenta);
  h_radon_1e->SetLineWidth(2);
  h_radon_1e->SetFillColor(kWhite);
  h_radon_1e1g_egamma->SetLineColor(kMagenta);
  h_radon_1e1g_egamma->SetMarkerColor(kMagenta);
  h_radon_1e1g_egamma->SetLineWidth(2);
  h_radon_1e1g_egamma->SetFillColor(kWhite);
  h_radon_1e1g_esum->SetLineColor(kMagenta);
  h_radon_1e1g_esum->SetMarkerColor(kMagenta);
  h_radon_1e1g_esum->SetLineWidth(2);
  h_radon_1e1g_esum->SetFillColor(kWhite);
  h_radon_1e2g_egamma_max->SetLineColor(kMagenta);
  h_radon_1e2g_egamma_max->SetMarkerColor(kMagenta);
  h_radon_1e2g_egamma_max->SetLineWidth(2);
  h_radon_1e2g_egamma_max->SetFillColor(kWhite);
  h_radon_1e2g_esum->SetLineColor(kMagenta);
  h_radon_1e2g_esum->SetMarkerColor(kMagenta);
  h_radon_1e2g_esum->SetLineWidth(2);
  h_radon_1e2g_esum->SetFillColor(kWhite);
  h_radon_1e3g_esum->SetLineColor(kMagenta);
  h_radon_1e3g_esum->SetMarkerColor(kMagenta);
  h_radon_1e3g_esum->SetLineWidth(2);
  h_radon_1e3g_esum->SetFillColor(kWhite);
  h_radon_2e1g_esum->SetLineColor(kMagenta);
  h_radon_2e1g_esum->SetMarkerColor(kMagenta);
  h_radon_2e1g_esum->SetLineWidth(2);
  h_radon_2e1g_esum->SetFillColor(kWhite);

  TCanvas *c1 = new TCanvas();

  Double_t xl1=.75, yl1=0.65, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_2nu_1e,"2#nu");
  leg->AddEntry(h_tl_1e,"^{208}Tl");
  leg->AddEntry(h_bi_1e,"^{214}Bi");
  leg->AddEntry(h_radon_1e,"Radon");
  leg->SetFillColor(kWhite);

  h_2nu_2e1g_esum->Draw("");
  h_tl_2e1g_esum->Draw("same");
  h_bi_2e1g_esum->Draw("same");
  h_radon_2e1g_esum->Draw("same");

  leg->Draw("same");

  return;

}
