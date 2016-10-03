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

void activity_accuracy_exposure () {

  TGraph *g_tl = new TGraph(7);
  TGraph *g_bi = new TGraph(7);

  //2nu
  // g->SetPoint(0,1./365.25,1e-1);
  // g->SetPoint(0,1./12,1e-1);
  // g->SetPoint(0,1.,5.9e-2);
  // // g->SetPoint(0,2.5,);

  //tl
  g_tl->SetPoint(0,1.,28.3);
  g_tl->SetPoint(1,2.,21.6);
  g_tl->SetPoint(2,3.,17.8);
  g_tl->SetPoint(3,6.,11.9);
  g_tl->SetPoint(4,12,8);
  g_tl->SetPoint(5,21,6.73);
  g_tl->SetPoint(6,30,5.5);

  // //bi
  g_bi->SetPoint(0,1.,18.2);
  g_bi->SetPoint(1,2.,10.9);
  g_bi->SetPoint(2,3.,8.29);
  g_bi->SetPoint(3,6.,5.8);
  g_bi->SetPoint(4,12,3.8);
  g_bi->SetPoint(5,21,3.1);
  g_bi->SetPoint(6,30,2.4);

  g_bi->GetXaxis()->SetTitleFont(62);
  g_bi->GetXaxis()->SetTitleSize(0.05);
  g_bi->GetXaxis()->SetTitleOffset(0.88);
  g_bi->GetXaxis()->SetTitle("Time [months]");
  // g_bi->GetYaxis()->SetTitle("Relative uncertainty on T_{1/2}^{2#nu} in %");
  // g_bi->GetYaxis()->SetTitle("Relative uncertainty on A(^{208}Tl) in %");
  g_bi->GetYaxis()->SetTitle("Relative uncertainty on the activity in %");
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
  g_bi->SetDrawOption("AP");

  g_tl->GetXaxis()->SetTitleFont(62);
  g_tl->GetXaxis()->SetTitleSize(0.05);
  g_tl->GetXaxis()->SetTitleOffset(0.88);
  g_tl->GetXaxis()->SetTitle("Time [months]");
  // g_tl->GetYaxis()->SetTitle("Relative uncertainty on T_{1/2}^{2#nu} in %");
  // g_tl->GetYaxis()->SetTitle("Relative uncertainty on A(^{208}Tl) in %");
  g_tl->GetYaxis()->SetTitle("Relative uncertainty on the activity in %");
  g_tl->GetYaxis()->SetTitleFont(62);
  g_tl->GetYaxis()->SetTitleSize(0.05);
  g_tl->GetYaxis()->SetTitleOffset(0.9);
  g_tl->GetYaxis()->SetTitleOffset(0.9);
  g_tl->SetMarkerStyle(20);
  g_tl->SetMarkerSize(1);
  // g_tl->SetMarkerColor(kBlue);
  // g_tl->SetLineColor(kBlue);
  // g_tl->SetFillColor(kBlue);
  g_tl->SetMarkerColor(kGreen+1);
  g_tl->SetLineColor(kGreen+1);
  g_tl->SetFillColor(kGreen+1);
  g_tl->SetDrawOption("AP");

  Double_t xl1=.75, yl1=0.65, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(g_bi,"^{214}Bi : 10 #muBq/kg");
  leg->AddEntry(g_tl,"^{208}Tl : 2 #muBq/kg");
  leg->SetFillColor(kWhite);

  g_tl->SetTitle("");
  g_tl->GetXaxis()->SetRangeUser(0,31);
  g_tl->GetYaxis()->SetRangeUser(0,31);
  g_tl->Draw("");
  g_bi->Draw("PLsame");
  leg->Draw("same");

  return;
}
