#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
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
#include "TRandom.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void activity_meas_accuracy()
{
  TH1F *h_2nu = new TH1F("h_2nu","",100,8.9,9.1);
  // TH1F *h_2nu = new TH1F("h_2nu","",100,8,10);
  TH1F *h_tl = new TH1F("h_tl","",100,0,4);
  TH1F *h_bi = new TH1F("h_bi","",100,0,20);

  char channel[200];

  TFile *f= new TFile("activity_measurement.root","RECREATE");
  TCanvas *c= new TCanvas();
  char file_name[200];

  sprintf(file_name,"./new_activity_measurements/measurements_tl_bi_1.75y.txt");

  h_2nu->SetName("activity_measurement_2nu");
  h_tl->SetName("activity_measurement_tl");
  h_bi->SetName("activity_measurement_bi");

  std::fstream data;
  data.open(file_name, std::fstream::in);

  int count = 1;
  std::string dump;
  double hl_2nu, act_tl, act_bi, act_radon;

  while(data >> dump >> hl_2nu >> dump >> act_tl >> dump >> act_bi >> dump >> act_radon)
    {
      h_2nu->Fill(hl_2nu);
      h_tl->Fill(act_tl);
      h_bi->Fill(act_bi);
      // std::cout << " 2nu " << hl_2nu << std::endl;
      count++;
    }

  // TRandom rdm;
  // unsigned int i_rdm = 0;
  // for(; i_rdm<800; ++i_rdm) {
  //   h_2nu->Fill(rdm.Gaus(9,3.94e-3));
  // }

  data.close();

  h_2nu->GetXaxis()->SetTitleFont(62);
  h_2nu->GetXaxis()->SetTitleSize(0.05);
  h_2nu->GetXaxis()->SetTitleOffset(0.88);
  h_2nu->GetYaxis()->SetTitle("Pseudo-experiments");
  h_2nu->GetYaxis()->SetTitleFont(62);
  h_2nu->GetYaxis()->SetTitleSize(0.05);
  h_2nu->GetYaxis()->SetTitleOffset(0.9);
  h_2nu->GetYaxis()->SetTitleOffset(0.9);
  h_2nu->SetMarkerStyle(20);
  h_2nu->SetMarkerSize(1);
  h_2nu->SetMarkerColor(kBlue);
  h_2nu->SetLineColor(kBlue);
  h_2nu->SetFillColor(kBlue);
  h_2nu->SetDrawOption("A");

  h_tl->GetXaxis()->SetTitleFont(62);
  h_tl->GetXaxis()->SetTitleSize(0.05);
  h_tl->GetXaxis()->SetTitleOffset(0.88);
  h_tl->GetYaxis()->SetTitle("Pseudo-experiments");
  h_tl->GetYaxis()->SetTitleFont(62);
  h_tl->GetYaxis()->SetTitleSize(0.05);
  h_tl->GetYaxis()->SetTitleOffset(0.9);
  h_tl->GetYaxis()->SetTitleOffset(0.9);
  h_tl->SetMarkerStyle(20);
  h_tl->SetMarkerSize(1);
  h_tl->SetMarkerColor(kBlue);
  h_tl->SetLineColor(kBlue);
  h_tl->SetFillColor(kBlue);
  h_tl->SetDrawOption("A");

  h_bi->GetXaxis()->SetTitleFont(62);
  h_bi->GetXaxis()->SetTitleSize(0.05);
  h_bi->GetXaxis()->SetTitleOffset(0.88);
  h_bi->GetYaxis()->SetTitle("Pseudo-experiments");
  h_bi->GetYaxis()->SetTitleFont(62);
  h_bi->GetYaxis()->SetTitleSize(0.05);
  h_bi->GetYaxis()->SetTitleOffset(0.9);
  h_bi->GetYaxis()->SetTitleOffset(0.9);
  h_bi->SetMarkerStyle(20);
  h_bi->SetMarkerSize(1);
  h_bi->SetMarkerColor(kBlue);
  h_bi->SetLineColor(kBlue);
  h_bi->SetFillColor(kBlue);
  h_bi->SetDrawOption("A");


  h_2nu->Write();
  h_tl->Write();
  h_bi->Write();

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
