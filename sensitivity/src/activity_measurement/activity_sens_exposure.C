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
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void activity_sens_exposure () {

  TGraph *g_bi = new TGraph(6);

  g_bi->SetPoint(0,0.05*12,60);
  g_bi->SetPoint(1,0.15*12,20);
  g_bi->SetPoint(2,0.3*12,10);
  g_bi->SetPoint(3,0.6*12,5);
  g_bi->SetPoint(4,1*12,3);
  g_bi->SetPoint(5,2.5*12,3/2.5);

  g_bi->GetXaxis()->SetTitleFont(62);
  g_bi->GetXaxis()->SetTitleSize(0.05);
  g_bi->GetXaxis()->SetTitleOffset(0.88);
  g_bi->GetXaxis()->SetTitle("Time [months]");
  // g_bi->GetYaxis()->SetTitle("Relative uncertainty on T_{1/2}^{2#nu} in %");
  // g_bi->GetYaxis()->SetTitle("Relative uncertainty on A(^{208}Tl) in %");
  g_bi->GetYaxis()->SetTitle("Activity measured @ 10% [#muBq/kg]");
  g_bi->GetYaxis()->SetTitleFont(62);
  g_bi->GetYaxis()->SetTitleSize(0.05);
  g_bi->GetYaxis()->SetTitleOffset(0.9);
  g_bi->GetYaxis()->SetTitleOffset(0.9);
  g_bi->SetMarkerStyle(20);
  g_bi->SetMarkerSize(1);
  // g_bi->SetMarkerColor(kBlue);
  // g_bi->SetLineColor(kBlue);
  // g_bi->SetFillColor(kBlue);
  g_bi->SetMarkerColor(kOrange-3);
  g_bi->SetLineColor(kOrange-3);
  g_bi->SetFillColor(kOrange-3);
  g_bi->SetDrawOption("P");

  TGraph *g_tl = new TGraph(6);

  g_tl->SetPoint(0,0.05*12,1.3/0.05);
  g_tl->SetPoint(1,0.15*12,1.3/0.15);
  g_tl->SetPoint(2,0.3*12,1.3/0.3);
  g_tl->SetPoint(3,0.6*12,1.3/0.6);
  g_tl->SetPoint(4,1*12,1.3);
  g_tl->SetPoint(5,2.5*12,1.3/2.5);

  g_tl->GetXaxis()->SetTitleFont(62);
  g_tl->GetXaxis()->SetTitleSize(0.05);
  g_tl->GetXaxis()->SetTitleOffset(0.88);
  g_tl->GetXaxis()->SetTitle("Time [months]");
  // g_tl->GetYaxis()->SetTitle("Relative uncertainty on T_{1/2}^{2#nu} in %");
  // g_tl->GetYaxis()->SetTitle("Relative uncertainty on A(^{208}Tl) in %");
  g_tl->GetYaxis()->SetTitle("Activity measured @ 10% [#muBq/kg]");
  g_tl->GetYaxis()->SetTitleFont(62);
  g_tl->GetYaxis()->SetTitleSize(0.05);
  g_tl->GetYaxis()->SetTitleOffset(0.9);
  g_tl->GetYaxis()->SetTitleOffset(0.9);
  g_tl->SetMarkerStyle(20);
  g_tl->SetMarkerSize(1);
  // g_tl->SetMarkerColor(kBlue);
  // g_tl->SetLineColor(kBlue);
  // g_tl->SetFillColor(kBlue);
  g_tl->SetMarkerColor(kGreen);
  g_tl->SetLineColor(kGreen);
  g_tl->SetFillColor(kGreen);
  g_tl->SetDrawOption("P");

  Double_t xl1=.75, yl1=0.65, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(g_bi,"^{214}Bi");
  leg->AddEntry(g_tl,"^{208}Tl");
  leg->SetFillColor(kWhite);

  g_tl->Draw("AP");
  g_bi->Draw("same");
  leg->Draw("same");

  return;
}
