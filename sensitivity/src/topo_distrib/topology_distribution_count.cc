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

void topology_distribuion_count()
{

  const double mass = 7.;
  const double exposure_sec = 2.5 * 3.14e7;
  const double exposure_y = 2.5;
  const double tracker_volume = 15.3; //m3
  const double halflife_2nu = 9e19; // years;
  const double Na = 6.022e23;
  const double M_Se = 0.082; //kg/mol more like 81.6 if we enrich at 90%
  const double const_se = log(2) * Na / M_Se; // years;

  const double activity_tl = 2e-6 ;
  const double activity_bi = 10e-6;
  const double activity_radon = 150e-6;

  TCanvas *c= new TCanvas();

  TFile f_2nu("2nu_topo_distrib.root");
  TFile f_tl("tl208_topo_distrib.root");
  TFile f_bi("bi214_topo_distrib.root");
  TFile f_radon("radon_topo_distrib.root");

  TH1F *h_2nu = (TH1F*)f_2nu.Get("h_topology_distribution");
  TH1F *h_tl = (TH1F*)f_tl.Get("h_topology_distribution");
  TH1F *h_bi = (TH1F*)f_bi.Get("h_topology_distribution");
  TH1F *h_radon = (TH1F*)f_radon.Get("h_topology_distribution");

  TH1F *new_h_2nu = new TH1F("h_2nu","h_2nu",78,0,78);
  TH1F *new_h_tl = new TH1F("h_tl","h_tl",78,0,78);
  TH1F *new_h_bi = new TH1F("h_bi","h_bi",78,0,78);
  TH1F *new_h_radon = new TH1F("h_radon","h_radon",78,0,78);

  h_2nu->Scale(mass*exposure_y/halflife_2nu*const_se);
  h_tl->Scale(activity_tl*exposure_sec*mass);
  h_bi->Scale(activity_bi*exposure_sec*mass);
  h_radon->Scale(activity_radon*exposure_sec*tracker_volume);

  h_2nu->GetXaxis()->SetTitle("Topology");
  h_2nu->GetXaxis()->SetTitleFont(62);
  h_2nu->GetXaxis()->SetTitleSize(0.05);
  h_2nu->GetXaxis()->SetTitleOffset(0.88);
  h_2nu->GetXaxis()->SetRangeUser(0.,13.);
  h_2nu->GetYaxis()->SetTitle("Counts");
  h_2nu->GetYaxis()->SetTitleFont(62);
  h_2nu->GetYaxis()->SetTitleSize(0.05);
  h_2nu->GetYaxis()->SetTitleOffset(0.9);
  h_2nu->GetYaxis()->SetTitleOffset(0.9);
  // h_2nu->GetYaxis()->SetRangeUser(1e-8,1.);
  h_2nu->SetMarkerStyle(20);
  h_2nu->SetMarkerSize(1);
  h_2nu->SetMarkerColor(kBlue);
  h_2nu->SetLineColor(kBlue);
  h_2nu->SetFillColor(kBlue);
  h_2nu->SetDrawOption("PEB");

  h_tl->GetXaxis()->SetTitle("Topology");
  h_tl->GetXaxis()->SetTitleFont(62);
  h_tl->GetXaxis()->SetTitleSize(0.05);
  h_tl->GetXaxis()->SetTitleOffset(0.88);
  h_tl->GetXaxis()->SetRangeUser(0.,13.);
  h_tl->GetYaxis()->SetTitle("Counts");
  h_tl->GetYaxis()->SetTitleFont(62);
  h_tl->GetYaxis()->SetTitleSize(0.05);
  h_tl->GetYaxis()->SetTitleOffset(0.9);
  h_tl->GetYaxis()->SetTitleOffset(0.9);
  // h_tl->GetYaxis()->SetRangeUser(1e-8,1.);
  h_tl->SetMarkerStyle(20);
  h_tl->SetMarkerSize(1);
  h_tl->SetMarkerColor(kGreen+1);
  h_tl->SetLineColor(kGreen+1);
  h_tl->SetFillColor(kGreen+1);
  h_tl->SetDrawOption("PEB");

  h_bi->GetXaxis()->SetTitle("Topology");
  h_bi->GetXaxis()->SetTitleFont(62);
  h_bi->GetXaxis()->SetTitleSize(0.05);
  h_bi->GetXaxis()->SetTitleOffset(0.88);
  h_bi->GetXaxis()->SetRangeUser(0.,13.);
  h_bi->GetYaxis()->SetTitle("Counts");
  h_bi->GetYaxis()->SetTitleFont(62);
  h_bi->GetYaxis()->SetTitleSize(0.05);
  h_bi->GetYaxis()->SetTitleOffset(0.9);
  h_bi->GetYaxis()->SetTitleOffset(0.9);
  // h_bi->GetYaxis()->SetRangeUser(1e-8,1.);
  h_bi->SetMarkerStyle(20);
  h_bi->SetMarkerSize(1);
  h_bi->SetMarkerColor(kOrange-3);
  h_bi->SetLineColor(kOrange-3);
  h_bi->SetFillColor(kOrange-3);
  h_bi->SetDrawOption("PEB");

  h_radon->GetXaxis()->SetTitle("Topology");
  h_radon->GetXaxis()->SetTitleFont(62);
  h_radon->GetXaxis()->SetTitleSize(0.05);
  h_radon->GetXaxis()->SetTitleOffset(0.88);
  h_radon->GetXaxis()->SetRangeUser(0.,13.);
  h_radon->GetYaxis()->SetTitle("Counts");
  h_radon->GetYaxis()->SetTitleFont(62);
  h_radon->GetYaxis()->SetTitleSize(0.05);
  h_radon->GetYaxis()->SetTitleOffset(0.9);
  h_radon->GetYaxis()->SetTitleOffset(0.9);
  // h_radon->GetYaxis()->SetRangeUser(1e-8,1.);
  h_radon->SetMarkerStyle(20);
  h_radon->SetMarkerSize(1);
  h_radon->SetMarkerColor(kMagenta);
  h_radon->SetLineColor(kMagenta);
  h_radon->SetFillColor(kMagenta);
  h_radon->SetDrawOption("PEB");


  new_h_2nu->GetXaxis()->SetTitle("Topology");
  new_h_2nu->GetXaxis()->SetTitleFont(62);
  new_h_2nu->GetXaxis()->SetTitleSize(0.05);
  new_h_2nu->GetXaxis()->SetTitleOffset(0.88);
  new_h_2nu->GetXaxis()->SetRangeUser(0.,78.);
  new_h_2nu->GetYaxis()->SetTitle("Counts");
  new_h_2nu->GetYaxis()->SetTitleFont(62);
  new_h_2nu->GetYaxis()->SetTitleSize(0.05);
  new_h_2nu->GetYaxis()->SetTitleOffset(0.9);
  new_h_2nu->GetYaxis()->SetTitleOffset(0.9);
  // new_h_2nu->GetYaxis()->SetRangeUser(1e-8,1.);
  new_h_2nu->SetMarkerStyle(20);
  new_h_2nu->SetMarkerSize(1);
  new_h_2nu->SetMarkerColor(kBlue);
  new_h_2nu->SetLineColor(kBlue);
  new_h_2nu->SetFillColor(kBlue);
  new_h_2nu->SetDrawOption("PEB");

  new_h_tl->GetXaxis()->SetTitle("Topology");
  new_h_tl->GetXaxis()->SetTitleFont(62);
  new_h_tl->GetXaxis()->SetTitleSize(0.05);
  new_h_tl->GetXaxis()->SetTitleOffset(0.88);
  new_h_tl->GetXaxis()->SetRangeUser(0.,78.);
  new_h_tl->GetYaxis()->SetTitle("Counts");
  new_h_tl->GetYaxis()->SetTitleFont(62);
  new_h_tl->GetYaxis()->SetTitleSize(0.05);
  new_h_tl->GetYaxis()->SetTitleOffset(0.9);
  new_h_tl->GetYaxis()->SetTitleOffset(0.9);
  // new_h_tl->GetYaxis()->SetRangeUser(1e-8,1.);
  new_h_tl->SetMarkerStyle(20);
  new_h_tl->SetMarkerSize(1);
  new_h_tl->SetMarkerColor(kGreen+1);
  new_h_tl->SetLineColor(kGreen+1);
  new_h_tl->SetFillColor(kGreen+1);
  new_h_tl->SetDrawOption("PEB");

  new_h_bi->GetXaxis()->SetTitle("Topology");
  new_h_bi->GetXaxis()->SetTitleFont(62);
  new_h_bi->GetXaxis()->SetTitleSize(0.05);
  new_h_bi->GetXaxis()->SetTitleOffset(0.88);
  new_h_bi->GetXaxis()->SetRangeUser(0.,78.);
  new_h_bi->GetYaxis()->SetTitle("Counts");
  new_h_bi->GetYaxis()->SetTitleFont(62);
  new_h_bi->GetYaxis()->SetTitleSize(0.05);
  new_h_bi->GetYaxis()->SetTitleOffset(0.9);
  new_h_bi->GetYaxis()->SetTitleOffset(0.9);
  // new_h_bi->GetYaxis()->SetRangeUser(1e-8,1.);
  new_h_bi->SetMarkerStyle(20);
  new_h_bi->SetMarkerSize(1);
  new_h_bi->SetMarkerColor(kOrange-3);
  new_h_bi->SetLineColor(kOrange-3);
  new_h_bi->SetFillColor(kOrange-3);
  new_h_bi->SetDrawOption("PEB");

  new_h_radon->GetXaxis()->SetTitle("Topology");
  new_h_radon->GetXaxis()->SetTitleFont(62);
  new_h_radon->GetXaxis()->SetTitleSize(0.05);
  new_h_radon->GetXaxis()->SetTitleOffset(0.88);
  new_h_radon->GetXaxis()->SetRangeUser(0.,78.);
  new_h_radon->GetYaxis()->SetTitle("Counts");
  new_h_radon->GetYaxis()->SetTitleFont(62);
  new_h_radon->GetYaxis()->SetTitleSize(0.05);
  new_h_radon->GetYaxis()->SetTitleOffset(0.9);
  new_h_radon->GetYaxis()->SetTitleOffset(0.9);
  // new_h_radon->GetYaxis()->SetRangeUser(1e-8,1.);
  new_h_radon->SetMarkerStyle(20);
  new_h_radon->SetMarkerSize(1);
  new_h_radon->SetMarkerColor(kMagenta);
  new_h_radon->SetLineColor(kMagenta);
  new_h_radon->SetFillColor(kMagenta);
  new_h_radon->SetDrawOption("PEB");

  gPad->SetLogy();

  h_2nu->SetStats(0);
  new_h_2nu->SetStats(0);

  for(unsigned int i = 0; i<13; ++i) {
    new_h_2nu->Fill(6*i+1,h_2nu->GetBinContent(i));
    new_h_tl->Fill(6*i+2,h_tl->GetBinContent(i));
    new_h_bi->Fill(6*i+3,h_bi->GetBinContent(i));
    new_h_radon->Fill(6*i+4,h_radon->GetBinContent(i));
    new_h_2nu->GetXaxis()->SetBinLabel(6*i+3,h_2nu->GetXaxis()->GetBinLabel(i));
  }

  // Double_t xl1=.65, yl1=0.7, xl2=0.9, yl2=0.9;
  // TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  // leg->AddEntry(h_2nu,"2#nu");
  // leg->AddEntry(h_tl,"^{208}Tl");
  // leg->AddEntry(h_bi,"^{214}Bi");
  // leg->AddEntry(h_radon,"Radon");
  // leg->SetFillColor(kWhite);

  // h_2nu->Draw("PEB");
  // h_tl->Draw("same");
  // h_bi->Draw("same");
  // h_radon->Draw("same");

  // leg->Draw("same");


  Double_t xl1=.65, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(new_h_2nu,"2#nu");
  leg->AddEntry(new_h_tl,"^{208}Tl");
  leg->AddEntry(new_h_bi,"^{214}Bi");
  leg->AddEntry(new_h_radon,"Radon");
  leg->SetFillColor(kWhite);

  new_h_2nu->Draw("B");
  new_h_tl->Draw("same");
  new_h_bi->Draw("same");
  new_h_radon->Draw("same");

  leg->Draw("same");

  TLine *line = new TLine(12,0,12,5e5);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  TLine *line2 = new TLine(2*6,0,2*6,5e5);
  line2->SetLineStyle(kDashed);
  line2->Draw("same");

  TLine *line3 = new TLine(3*6,0,3*6,5e5);
  line3->SetLineStyle(kDashed);
  line3->Draw("same");

  TLine *line4 = new TLine(4*6,0,4*6,5e5);
  line4->SetLineStyle(kDashed);
  line4->Draw("same");

  TLine *line5 = new TLine(5*6,0,5*6,5e5);
  line5->SetLineStyle(kDashed);
  line5->Draw("same");

  TLine *line6 = new TLine(6*6,0,6*6,3e3);
  line6->SetLineStyle(kDashed);
  line6->Draw("same");

  TLine *line7 = new TLine(7*6,0,7*6,3e3);
  line7->SetLineStyle(kDashed);
  line7->Draw("same");

  TLine *line8 = new TLine(8*6,0,8*6,3e3);
  line8->SetLineStyle(kDashed);
  line8->Draw("same");

  TLine *line9 = new TLine(9*6,0,9*6,3e3);
  line9->SetLineStyle(kDashed);
  line9->Draw("same");

  TLine *line10 = new TLine(10*6,0,10*6,5e5);
  line10->SetLineStyle(kDashed);
  line10->Draw("same");

  TLine *line10 = new TLine(10*6,0,10*6,5e5);
  line10->SetLineStyle(kDashed);
  line10->Draw("same");

  TLine *line11 = new TLine(11*6,0,11*6,5e5);
  line11->SetLineStyle(kDashed);
  line11->Draw("same");

  // TLine *line12 = new TLine(12*6,0,12*6,10e5);
  // line12->SetLineStyle(kDashed);
  // line12->Draw("same");


  return;
}
