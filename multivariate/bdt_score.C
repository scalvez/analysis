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

void bdt_score()
{

  bool counts = true;

  TFile * f_0nu = TFile::Open("./bdt_scores/test_consistency/test_0nu_Ecut.root");
  // TFile * f_2nu = TFile::Open("2nu_1M.root");
  TFile * f_2nu = TFile::Open("./bdt_scores/test_consistency/test_2nu_full_Ecut.root");
  TFile * f_tl208 = TFile::Open("./bdt_scores/test_consistency/test_tl208_Ecut.root");
  TFile * f_bi214 = TFile::Open("./bdt_scores/test_consistency/test_bi214_Ecut.root");
  TFile * f_radon = TFile::Open("./bdt_scores/test_consistency/test_radon_Ecut.root");

  TH1F *h_0nu_bdt = (TH1F*)f_0nu->Get("MVA_BDT");
  TH1F *h_2nu_bdt = (TH1F*)f_2nu->Get("MVA_BDT");
  TH1F *h_tl208_bdt = (TH1F*)f_tl208->Get("MVA_BDT");
  TH1F *h_bi214_bdt = (TH1F*)f_bi214->Get("MVA_BDT");
  TH1F *h_radon_bdt = (TH1F*)f_radon->Get("MVA_BDT");

  // TH1F *h_0nu_bdt = (TH1F*)f_0nu->Get("MVA_BDTG");
  // TH1F *h_2nu_bdt = (TH1F*)f_2nu->Get("MVA_BDTG");
  // TH1F *h_tl208_bdt = (TH1F*)f_tl208->Get("MVA_BDTG");
  // TH1F *h_bi214_bdt = (TH1F*)f_bi214->Get("MVA_BDTG");
  // TH1F *h_radon_bdt = (TH1F*)f_radon->Get("MVA_BDTG");

  h_0nu_bdt->Sumw2();
  h_2nu_bdt->Sumw2();
  h_tl208_bdt->Sumw2();
  h_bi214_bdt->Sumw2();
  h_radon_bdt->Sumw2();

  if(counts)
    TFile *f_output= new TFile("./bdt_scores/test_consistency/bdt_scores_counts_test_Ecut.root","RECREATE");
  else
    TFile *f_output= new TFile("./bdt_scores/test_consistency/bdt_scores_test_Ecut.root","RECREATE");

  if(counts)
    h_0nu_bdt->Scale(1./h_0nu_bdt->GetEntries());
  else
    h_0nu_bdt->Scale(1./h_0nu_bdt->GetEntries());
  h_0nu_bdt->SetLineColor(kRed);
  h_0nu_bdt->SetLineWidth(2);
  // h_0nu_bdt->SetFillColor(kRed);
  h_0nu_bdt->SetTitle("0#nu;BDT score; Probability");
  h_0nu_bdt->SetName("0nu");
  // h_0nu->Rebin();

  if(counts)
    h_2nu_bdt->Scale(conf_sens::N_source_2nu/h_2nu_bdt->GetEntries());
  else
    h_2nu_bdt->Scale(1./h_2nu_bdt->GetEntries());
  h_2nu_bdt->SetLineColor(kBlue);
  h_2nu_bdt->SetLineWidth(2);
  // h_2nu_bdt->SetFillColor(kBlue);
  h_2nu_bdt->SetTitle("2#nu;BDT score; Probability");
  h_2nu_bdt->SetName("2nu");
  // h_2nu->Rebin();

  if(counts)
    h_tl208_bdt->Scale(conf_sens::N_source_tl208/h_tl208_bdt->GetEntries());
  else
    h_tl208_bdt->Scale(1./h_tl208_bdt->GetEntries());
  h_tl208_bdt->SetLineColor(kGreen+1);
  h_tl208_bdt->SetLineWidth(2);
  // h_tl208_bdt->SetFillColor(kGreen+1);
  h_tl208_bdt->SetTitle("^{208}Tl;BDT score; Probability");
  h_tl208_bdt->SetName("tl208");
  // h_tl208->Rebin();

  if(counts)
    h_bi214_bdt->Scale(conf_sens::N_source_bi214/h_bi214_bdt->GetEntries());
  else
    h_bi214_bdt->Scale(1./h_bi214_bdt->GetEntries());
  h_bi214_bdt->SetLineColor(kOrange-3);
  h_bi214_bdt->SetLineWidth(2);
  // h_bi214_bdt->SetFillColor(kOrange-3);
  h_bi214_bdt->SetTitle("^{214}Bi;BDT score; Probability");
  h_bi214_bdt->SetName("bi214");
  // h_bi214->Rebin();

  if(counts)
    h_radon_bdt->Scale(conf_sens::N_source_radon/h_radon_bdt->GetEntries());
  else
    h_radon_bdt->Scale(1./h_radon_bdt->GetEntries());
  h_radon_bdt->SetLineColor(kMagenta);
  h_radon_bdt->SetLineWidth(2);
  // h_radon_bdt->SetFillColor(kMagenta);
  h_radon_bdt->SetTitle("Radon;BDT score; Probability");
  h_radon_bdt->SetName("radon");
  // h_radon->Rebin();

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu_bdt,"0#nu");
  leg->AddEntry(h_2nu_bdt,"2#nu");
  leg->AddEntry(h_tl208_bdt,"^{208}Tl");
  leg->AddEntry(h_bi214_bdt,"^{214}Bi");
  leg->AddEntry(h_radon_bdt,"Radon");
  leg->SetFillColor(kWhite);

  // h_0nu_bdt->Draw("");
  // h_2nu_bdt->Draw("same");
  // h_tl208_bdt->Draw("same");
  // h_bi214_bdt->Draw("same");
  // h_radon_bdt->Draw("same");
  // leg->Draw("same");

  h_0nu_bdt->Write();
  h_2nu_bdt->Write();
  h_tl208_bdt->Write();
  h_bi214_bdt->Write();
  h_radon_bdt->Write();

  THStack *hs = new THStack("hs","BDT scores");
  hs->Add(h_2nu_bdt);
  hs->Add(h_tl208_bdt);
  hs->Add(h_bi214_bdt);
  hs->Add(h_radon_bdt);
  hs->Add(h_0nu_bdt);

  hs->SetTitle("BDT;BDT score;Probability");
  hs->Write();

  TCanvas * c1 = new TCanvas();
  c1->cd();
  // hs->Draw();
  h_2nu_bdt->Draw("");
  h_tl208_bdt->Draw("same");
  h_bi214_bdt->Draw("same");
  h_radon_bdt->Draw("same");
  h_0nu_bdt->Draw("same");

  leg->Draw("same");

  c1->Write();

  TH1F *h_data = new TH1F("h_data","h_data",202,-1.01,1.01) ;

  for(unsigned int i = 1; i<=202; ++i) {

    h_data->SetBinContent(i,h_0nu_bdt->GetBinContent(i) +
                          h_2nu_bdt->GetBinContent(i)   +
                          h_tl208_bdt->GetBinContent(i) +
                          h_bi214_bdt->GetBinContent(i) +
                          h_radon_bdt->GetBinContent(i)
                          );
  }

  h_data->Write();

  // gStyle->SetTitleFontSize(0.08);
  // gPad->SetLogy();
  // // gPad->SetPad(0.1,0.2,0.9,0.95);
  // gPad->SetBottomMargin(0.155);
  // gPad->SetTopMargin(1);
  // gPad->SetRightMargin(1000);

  return;
}

/*
    // --- Example of simple scan

    TFile f("simple_scans.root");
    TGraphErrors *h_82 =(TGraphErrors*)f.Get("82");
    Double_t xl1=.05, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
    leg->AddEntry(h0nu,"0#nu2#beta no cuts");
    leg->AddEntry(h2nu,"2#nu2#beta no cuts");


  */
