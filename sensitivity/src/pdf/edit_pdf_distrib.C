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
  TH1F *h_2nu_2e_emax = (TH1F*)f_2nu.Get("2e_electron_maximal_energy");
  TH1F *h_2nu_2e_esum = (TH1F*)f_2nu.Get("2e_electrons_energy_sum");
  TH1F *h_2nu_2e_pint = (TH1F*)f_2nu.Get("2e_electrons_internal_probability");
  TH1F *h_2nu_2e_cos = (TH1F*)f_2nu.Get("2e_electrons_cos_angle");

  TH1F *h_tl_1e = (TH1F*)f_tl.Get("1e_electron_energy");
  TH1F *h_tl_2e_emax = (TH1F*)f_tl.Get("2e_electron_maximal_energy");
  TH1F *h_tl_2e_esum = (TH1F*)f_tl.Get("2e_electrons_energy_sum");
  TH1F *h_tl_2e_pint = (TH1F*)f_tl.Get("2e_electrons_internal_probability");
  TH1F *h_tl_2e_cos = (TH1F*)f_tl.Get("2e_electrons_cos_angle");

  TH1F *h_bi_1e = (TH1F*)f_bi.Get("1e_electron_energy");
  TH1F *h_bi_2e_emax = (TH1F*)f_bi.Get("2e_electron_maximal_energy");
  TH1F *h_bi_2e_esum = (TH1F*)f_bi.Get("2e_electrons_energy_sum");
  TH1F *h_bi_2e_pint = (TH1F*)f_bi.Get("2e_electrons_internal_probability");
  TH1F *h_bi_2e_cos = (TH1F*)f_bi.Get("2e_electrons_cos_angle");

  TH1F *h_radon_1e = (TH1F*)f_radon.Get("1e_electron_energy");
  TH1F *h_radon_2e_emax = (TH1F*)f_radon.Get("2e_electron_maximal_energy");
  TH1F *h_radon_2e_esum = (TH1F*)f_radon.Get("2e_electrons_energy_sum");
  TH1F *h_radon_2e_pint = (TH1F*)f_radon.Get("2e_electrons_internal_probability");
  TH1F *h_radon_2e_cos = (TH1F*)f_radon.Get("2e_electrons_cos_angle");

  h_2nu_1e->SetStats(0);
  h_2nu_2e_emax->SetStats(0);
  h_2nu_2e_esum->SetStats(0);
  h_2nu_2e_pint->SetStats(0);
  h_2nu_2e_cos ->SetStats(0);

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

  h_2nu_2e_emax->SetTitle("2e : Electron maximal energy");
  h_2nu_2e_emax->GetXaxis()->SetTitle("Electron maximal energy");
  h_2nu_2e_emax->GetXaxis()->SetTitleFont(62);
  h_2nu_2e_emax->GetXaxis()->SetTitleSize(0.05);
  h_2nu_2e_emax->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_2e_emax->GetYaxis()->SetTitle("Probability");
  h_2nu_2e_emax->GetYaxis()->SetTitleFont(62);
  h_2nu_2e_emax->GetYaxis()->SetTitleSize(0.05);
  h_2nu_2e_emax->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_2e_emax->SetLineColor(kBlue);
  h_2nu_2e_emax->SetMarkerColor(kBlue);
  h_2nu_2e_emax->SetLineWidth(2);
  h_2nu_2e_emax->SetFillColor(kWhite);

  h_2nu_2e_esum->SetTitle("2e : Electrons energy sum");
  h_2nu_2e_esum->GetXaxis()->SetTitle("Electrons energy sum");
  h_2nu_2e_esum->GetXaxis()->SetTitleFont(62);
  h_2nu_2e_esum->GetXaxis()->SetTitleSize(0.05);
  h_2nu_2e_esum->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_2e_esum->GetYaxis()->SetTitle("Probability");
  h_2nu_2e_esum->GetYaxis()->SetTitleFont(62);
  h_2nu_2e_esum->GetYaxis()->SetTitleSize(0.05);
  h_2nu_2e_esum->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_2e_esum->SetLineColor(kBlue);
  h_2nu_2e_esum->SetMarkerColor(kBlue);
  h_2nu_2e_esum->SetLineWidth(2);
  h_2nu_2e_esum->SetFillColor(kWhite);

  h_2nu_2e_pint->SetTitle("2e : Electrons internal probability");
  h_2nu_2e_pint->GetXaxis()->SetTitleFont(62);
  h_2nu_2e_pint->GetXaxis()->SetTitleSize(0.05);
  h_2nu_2e_pint->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_2e_pint->GetYaxis()->SetTitle("Probability");
  h_2nu_2e_pint->GetYaxis()->SetTitleFont(62);
  h_2nu_2e_pint->GetYaxis()->SetTitleSize(0.05);
  h_2nu_2e_pint->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_2e_pint->SetLineColor(kBlue);
  h_2nu_2e_pint->SetMarkerColor(kBlue);
  h_2nu_2e_pint->SetLineWidth(2);
  h_2nu_2e_pint->SetFillColor(kWhite);

  h_2nu_2e_cos->GetXaxis()->SetTitle("2e : Electrons Cos#theta");
  h_2nu_2e_cos->GetXaxis()->SetTitle("Electrons Cos#theta");
  h_2nu_2e_cos->GetXaxis()->SetTitleFont(62);
  h_2nu_2e_cos->GetXaxis()->SetTitleSize(0.05);
  h_2nu_2e_cos->GetXaxis()->SetTitleOffset(0.90);
  h_2nu_2e_cos->GetYaxis()->SetTitle("Probability");
  h_2nu_2e_cos->GetYaxis()->SetTitleFont(62);
  h_2nu_2e_cos->GetYaxis()->SetTitleSize(0.05);
  h_2nu_2e_cos->GetYaxis()->SetTitleOffset(0.90);
  h_2nu_2e_cos->SetLineColor(kBlue);
  h_2nu_2e_cos->SetMarkerColor(kBlue);
  h_2nu_2e_cos->SetLineWidth(2);
  h_2nu_2e_cos->SetFillColor(kWhite);

  h_tl_1e->SetLineColor(kGreen+1);
  h_tl_1e->SetMarkerColor(kGreen+1);
  h_tl_1e->SetLineWidth(2);
  h_tl_1e->SetFillColor(kWhite);
  h_tl_2e_esum->SetLineColor(kGreen+1);
  h_tl_2e_esum->SetMarkerColor(kGreen+1);
  h_tl_2e_esum->SetLineWidth(2);
  h_tl_2e_esum->SetFillColor(kWhite);
  h_tl_2e_emax->SetLineColor(kGreen+1);
  h_tl_2e_emax->SetMarkerColor(kGreen+1);
  h_tl_2e_emax->SetLineWidth(2);
  h_tl_2e_emax->SetFillColor(kWhite);
  h_tl_2e_pint->SetLineColor(kGreen+1);
  h_tl_2e_pint->SetMarkerColor(kGreen+1);
  h_tl_2e_pint->SetLineWidth(2);
  h_tl_2e_pint->SetFillColor(kWhite);
  h_tl_2e_cos->SetLineColor(kGreen+1);
  h_tl_2e_cos->SetMarkerColor(kGreen+1);
  h_tl_2e_cos->SetLineWidth(2);
  h_tl_2e_cos->SetFillColor(kWhite);

  h_bi_1e->SetLineColor(kOrange-3);
  h_bi_1e->SetMarkerColor(kOrange-3);
  h_bi_1e->SetLineWidth(2);
  h_bi_1e->SetFillColor(kWhite);
  h_bi_2e_esum->SetLineColor(kOrange-3);
  h_bi_2e_esum->SetMarkerColor(kOrange-3);
  h_bi_2e_esum->SetLineWidth(2);
  h_bi_2e_esum->SetFillColor(kWhite);
  h_bi_2e_emax->SetLineColor(kOrange-3);
  h_bi_2e_emax->SetMarkerColor(kOrange-3);
  h_bi_2e_emax->SetLineWidth(2);
  h_bi_2e_emax->SetFillColor(kWhite);
  h_bi_2e_pint->SetLineColor(kOrange-3);
  h_bi_2e_pint->SetMarkerColor(kOrange-3);
  h_bi_2e_pint->SetLineWidth(2);
  h_bi_2e_pint->SetFillColor(kWhite);
  h_bi_2e_cos->SetLineColor(kOrange-3);
  h_bi_2e_cos->SetMarkerColor(kOrange-3);
  h_bi_2e_cos->SetLineWidth(2);
  h_bi_2e_cos->SetFillColor(kWhite);

  h_radon_1e->SetLineColor(kMagenta);
  h_radon_1e->SetMarkerColor(kMagenta);
  h_radon_1e->SetLineWidth(2);
  h_radon_1e->SetFillColor(kWhite);
  h_radon_2e_esum->SetLineColor(kMagenta);
  h_radon_2e_esum->SetMarkerColor(kMagenta);
  h_radon_2e_esum->SetLineWidth(2);
  h_radon_2e_esum->SetFillColor(kWhite);
  h_radon_2e_emax->SetLineColor(kMagenta);
  h_radon_2e_emax->SetMarkerColor(kMagenta);
  h_radon_2e_emax->SetLineWidth(2);
  h_radon_2e_emax->SetFillColor(kWhite);
  h_radon_2e_pint->SetLineColor(kMagenta);
  h_radon_2e_pint->SetMarkerColor(kMagenta);
  h_radon_2e_pint->SetLineWidth(2);
  h_radon_2e_pint->SetFillColor(kWhite);
  h_radon_2e_cos->SetLineColor(kMagenta);
  h_radon_2e_cos->SetMarkerColor(kMagenta);
  h_radon_2e_cos->SetLineWidth(2);
  h_radon_2e_cos->SetFillColor(kWhite);

  TCanvas *c1 = new TCanvas();

  Double_t xl1=.75, yl1=0.65, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_2nu_1e,"2#nu");
  leg->AddEntry(h_tl_1e,"^{208}Tl");
  leg->AddEntry(h_bi_1e,"^{214}Bi");
  leg->AddEntry(h_radon_1e,"Radon");
  leg->SetFillColor(kWhite);

  h_2nu_2e_cos->Draw("");
  h_tl_2e_cos->Draw("same");
  h_bi_2e_cos->Draw("same");
  h_radon_2e_cos->Draw("same");

  leg->Draw("same");


  return;

}
