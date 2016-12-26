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

void spectra_2e_simple()
{
  bool counts = true;

  // TFile * f_0nu = TFile::Open("../data/trees_source_2e_application/0nu.root");
  // TFile * f_2nu = TFile::Open("../data/trees_source_2e_application/2nu.root");

  TFile * f_0nu = TFile::Open("../data/trees_topo_cuts_2e_application/0nu.root");
  TFile * f_2nu = TFile::Open("../data/trees_topo_cuts_2e_application/2nu.root");

  // TFile * f_0nu = TFile::Open("../data/trees_source_2e/0nu.root");
  // TFile * f_2nu = TFile::Open("../data/trees_source_2e/2nu.root");

  TTree *tree_0nu = (TTree*)f_0nu->Get("snemodata");
  TTree *tree_2nu = (TTree*)f_2nu->Get("snemodata");

  TH1F *h_0nu = new TH1F("0nu","0nu",100,0,5);
  TH1F *h_2nu = new TH1F("2nu","2nu",100,0,5);

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

  h_0nu->Sumw2();
  h_2nu->Sumw2();

  if (counts)
    TFile *f_output= new TFile("spectra_counts_simple.root","RECREATE");
  else
    TFile *f_output= new TFile("spectra_counts_simple.root","RECREATE");

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

  TH1F *h_data = new TH1F("h_data","h_data",100,0,5) ;
  for(unsigned int i = 1; i<=100; ++i) {

    h_data->SetBinContent(i,h_0nu->GetBinContent(i) +
                          h_2nu->GetBinContent(i));
  }

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu,"0#nu");
  leg->AddEntry(h_2nu,"2#nu");
  leg->SetFillColor(kWhite);

  h_0nu->Draw("");
  h_2nu->Draw("same");
  leg->Draw("same");

  h_0nu->Write();
  h_2nu->Write();

  h_data->Write();

  THStack *hs = new THStack("hs","Energy [keV]");
  hs->Add(h_2nu);
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
