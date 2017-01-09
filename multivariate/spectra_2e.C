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
#include "TTree.h"
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

void spectra_2e()
{
  bool counts = false;

  // TFile * f_0nu = TFile::Open("$SW_WORK_DIR/multivariate_analysis/data_0nu_mm.root");
  // // TFile * f_2nu = TFile::Open("2nu_1M.root");
  // TFile * f_2nu = TFile::Open("$SW_WORK_DIR/multivariate_analysis/data_2nu.root");
  // TFile * f_tl208 = TFile::Open("$SW_WORK_DIR/multivariate_analysis/data_tl208.root");
  // TFile * f_bi214 = TFile::Open("$SW_WORK_DIR/multivariate_analysis/data_bi214.root");
  // // TFile * f_radon = TFile::Open("../data/trees_source_2e_application/radon.root");

  TFile * f_0nu = TFile::Open("../data/trees_source_2e_application/0nu_1M.root");
  // TFile * f_2nu = TFile::Open("../data/trees_source_2e_application/2nu_1M.root");
  TFile * f_2nu = TFile::Open("../data/trees_source_2e_application/2nu_full_2MeV.root");
  TFile * f_tl208 = TFile::Open("../data/trees_source_2e_application/tl208.root");
  TFile * f_bi214 = TFile::Open("../data/trees_source_2e_application/bi214.root");
  TFile * f_radon = TFile::Open("../data/trees_source_2e_application/radon.root");

  // TFile * f_0nu = TFile::Open("../data/trees_topo_cuts_2e_application/0nu.root");
  // TFile * f_2nu = TFile::Open("../data/trees_topo_cuts_2e_application/2nu_1M.root");
  // TFile * f_tl208 = TFile::Open("../data/trees_topo_cuts_2e_application/tl208.root");
  // TFile * f_bi214 = TFile::Open("../data/trees_topo_cuts_2e_application/bi214.root");
  // TFile * f_radon = TFile::Open("../data/trees_topo_cuts_2e_application/radon.root");

  // TFile * f_0nu = TFile::Open("../data/trees_source_2e/0nu.root");
  // TFile * f_2nu = TFile::Open("../data/trees_source_2e/2nu.root");
  // TFile * f_tl208 = TFile::Open("../data/trees_source_2e/tl208.root");
  // TFile * f_bi214 = TFile::Open("../data/trees_source_2e/bi214.root");
  // TFile * f_radon = TFile::Open("../data/trees_source_2e/radon.root");

  TTree *tree_0nu = (TTree*)f_0nu->Get("snemodata");
  TTree *tree_2nu = (TTree*)f_2nu->Get("snemodata");
  TTree *tree_tl208 = (TTree*)f_tl208->Get("snemodata");
  TTree *tree_bi214 = (TTree*)f_bi214->Get("snemodata");
  TTree *tree_radon = (TTree*)f_radon->Get("snemodata");

  TH1F *h_0nu = new TH1F("0nu","0nu",100,0,5);
  TH1F *h_2nu = new TH1F("2nu","2nu",100,0,5);
  TH1F *h_tl208 = new TH1F("tl208","tl208",100,0,5);
  TH1F *h_bi214 = new TH1F("bi214","bi214",100,0,5);
  TH1F *h_radon = new TH1F("radon","radon",100,0,5);

  double electrons_energy_sum_0nu = 0;
  tree_0nu->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum_0nu);
  int nentries_0nu = tree_0nu->GetEntriesFast();
  for(int i = 0; i< nentries_0nu; ++i) {
    tree_0nu->GetEntry(i);
    // std::cout << ergy_0nu << std::endl;
    h_0nu->Fill(electrons_energy_sum_0nu);
  }

  double electrons_energy_sum_2nu = 0;
  tree_2nu->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum_2nu);
  int nentries_2nu = tree_2nu->GetEntriesFast();
  for(int i = 0; i< nentries_2nu; ++i) {
  // for(int i = 0; i< 100000; ++i) {
    tree_2nu->GetEntry(i);
    // std::cout << electrons_energy_sum_2nu << std::endl;
    h_2nu->Fill(electrons_energy_sum_2nu);
  }

  double electrons_energy_sum_tl208 = 0;
  tree_tl208->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum_tl208);
  int nentries_tl208 = tree_tl208->GetEntriesFast();
  for(int i = 0; i< nentries_tl208; ++i) {
    tree_tl208->GetEntry(i);
    // std::cout << electrons_energy_sum_tl208 << std::endl;
    h_tl208->Fill(electrons_energy_sum_tl208);
  }

  double electrons_energy_sum_bi214 = 0;
  tree_bi214->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum_bi214);
  int nentries_bi214 = tree_bi214->GetEntriesFast();
  for(int i = 0; i< nentries_bi214; ++i) {
    tree_bi214->GetEntry(i);
    // std::cout << electrons_energy_sum_bi214 << std::endl;
    h_bi214->Fill(electrons_energy_sum_bi214);
  }

  double electrons_energy_sum_radon = 0;
  tree_radon->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum_radon);
  int nentries_radon = tree_radon->GetEntriesFast();
  for(int i = 0; i< nentries_radon; ++i) {
    tree_radon->GetEntry(i);
    // std::cout << electrons_energy_sum_radon << std::endl;
    h_radon->Fill(electrons_energy_sum_radon);
  }

  h_0nu->Sumw2();
  h_2nu->Sumw2();
  h_tl208->Sumw2();
  h_bi214->Sumw2();
  h_radon->Sumw2();

  if (counts)
    TFile *f_output= new TFile("./spectra/spectra_counts.root","RECREATE");
  else
    TFile *f_output= new TFile("./spectra/spectra.root","RECREATE");

  if(counts)
    h_0nu->Scale(1./h_0nu->GetEntries());
  else
    h_0nu->Scale(1./h_0nu->GetEntries());
  h_0nu->SetLineColor(kRed);
  h_0nu->SetLineWidth(2);
  // h_0nu->SetFillColor(kRed);
  h_0nu->SetTitle("0#nu;Energy [keV]; Probability");
  h_0nu->SetName("0nu");
  // h_0nu->Rebin();

  if(counts)
    h_2nu->Scale(conf_sens::N_source_2nu/h_2nu->GetEntries());
  else
    h_2nu->Scale(1./h_2nu->GetEntries());
  h_2nu->SetLineColor(kBlue);
  h_2nu->SetLineWidth(2);
  // h_2nu->SetFillColor(kBlue);
  h_2nu->SetTitle("2#nu;Energy [keV]; Probability");
  h_2nu->SetName("2nu");
  // h_2nu->Rebin();

  if(counts)
    h_tl208->Scale(conf_sens::N_source_tl208/h_tl208->GetEntries());
  else
    h_tl208->Scale(1./h_tl208->GetEntries());
  h_tl208->SetLineColor(kGreen+1);
  h_tl208->SetLineWidth(2);
  // h_tl208->SetFillColor(kGreen+1);
  h_tl208->SetTitle("^{208}Tl;Energy [keV]; Probability");
  h_tl208->SetName("tl208");
  // h_tl208->Rebin();

  if(counts)
    h_bi214->Scale(conf_sens::N_source_bi214/h_bi214->GetEntries());
  else
    h_bi214->Scale(1./h_bi214->GetEntries());
  h_bi214->SetLineColor(kOrange-3);
  h_bi214->SetLineWidth(2);
  // h_bi214->SetFillColor(kOrange-3);
  h_bi214->SetTitle("^{214}Bi;Energy [keV]; Probability");
  h_bi214->SetName("bi214");
  // h_bi214->Rebin();

  if(counts)
    h_radon->Scale(conf_sens::N_source_radon/h_radon->GetEntries());
  else
    h_radon->Scale(1./h_radon->GetEntries());
  h_radon->SetLineColor(kMagenta);
  h_radon->SetLineWidth(2);
  // h_radon->SetFillColor(kMagenta);
  h_radon->SetTitle("Radon;Energy [keV]; Probability");
  h_radon->SetName("radon");
  // h_radon->Rebin();

  TH1F *h_data = new TH1F("h_data","h_data",100,0,5) ;
  for(unsigned int i = 1; i<=100; ++i) {

    h_data->SetBinContent(i,h_0nu->GetBinContent(i) +
                          h_2nu->GetBinContent(i)   +
                          h_tl208->GetBinContent(i) +
                          h_bi214->GetBinContent(i)
                          + h_radon->GetBinContent(i)
                          );
  }

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu,"0#nu");
  leg->AddEntry(h_2nu,"2#nu");
  leg->AddEntry(h_tl208,"^{208}Tl");
  leg->AddEntry(h_bi214,"^{214}Bi");
  leg->AddEntry(h_radon,"Radon");
  leg->SetFillColor(kWhite);

  h_0nu->Draw("");
  h_2nu->Draw("same");
  h_tl208->Draw("same");
  h_bi214->Draw("same");
  h_radon->Draw("same");
  leg->Draw("same");

  h_0nu->Write();
  h_2nu->Write();
  h_tl208->Write();
  h_bi214->Write();
  h_radon->Write();

  h_data->Write();

  THStack *hs = new THStack("hs","Energy [keV]");
  hs->Add(h_2nu);
  hs->Add(h_tl208);
  hs->Add(h_bi214);
  hs->Add(h_radon);
  hs->Add(h_0nu);

  hs->SetTitle("Spectrum;Energy;Probability");
  hs->Write();

  // TCanvas * c1 = new TCanvas();
  // c1->cd();
  // // hs->Draw();
  // // leg->Draw("same");
  // c1->Write();


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
