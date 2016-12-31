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

void bdt_score_simple()
{
bool counts = true;

  TFile * f_0nu = TFile::Open("0nu.root");
  TFile * f_2nu = TFile::Open("2nu_1M.root");

  TH1F *h_0nu_bdt = (TH1F*)f_0nu->Get("MVA_BDT");
  TH1F *h_2nu_bdt = (TH1F*)f_2nu->Get("MVA_BDT");

  // TH1F *h_0nu_bdt = (TH1F*)f_0nu->Get("MVA_BDTG");
  // TH1F *h_2nu_bdt = (TH1F*)f_2nu->Get("MVA_BDTG");

  h_0nu_bdt->Sumw2();
  h_2nu_bdt->Sumw2();

  if(counts)
    TFile *f_output= new TFile("bdt_scores_counts_simple.root","RECREATE");
  else
    TFile *f_output= new TFile("bdt_scores_simple.root","RECREATE");

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

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu_bdt,"0#nu");
  leg->AddEntry(h_2nu_bdt,"2#nu");
  leg->SetFillColor(kWhite);

  // h_0nu_bdt->Draw("");
  // h_2nu_bdt->Draw("same");
  // h_tl208_bdt->Draw("same");
  // h_bi214_bdt->Draw("same");
  // h_radon_bdt->Draw("same");
  // leg->Draw("same");

  h_0nu_bdt->Write();
  h_2nu_bdt->Write();

  THStack *hs = new THStack("hs","BDT scores");
  hs->Add(h_2nu_bdt);
  hs->Add(h_0nu_bdt);

  hs->SetTitle("BDT;BDT score;Probability");
  hs->Write();

  TCanvas * c1 = new TCanvas();
  c1->cd();
  // hs->Draw();
  h_0nu_bdt->Draw("");
  h_2nu_bdt->Draw("same");

  leg->Draw("same");

  c1->Write();

  // TH1F *h_data = new TH1F("h_data","h_data",100,-0.8,0.8) ;
  TH1F *h_data = new TH1F("h_data","h_data",100,-1,1) ;

  for(unsigned int i = 1; i<=100; ++i) {

    h_data->SetBinContent(i,h_0nu_bdt->GetBinContent(i) +
                          h_2nu_bdt->GetBinContent(i));
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
