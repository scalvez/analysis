#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>


void edit_sensitivity() {

  TFile *f = TFile::Open("sensitivity.root");

  TGraph *g_50 = (TGraph*)f->Get("hl_radon_exp_50");
  TGraph *g_150 = (TGraph*)f->Get("hl_radon_exp_150");
  TGraph *g_500 = (TGraph*)f->Get("hl_radon_exp_500");

  Double_t xl1=.65, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetHeader("Radon contamination :");
  leg->AddEntry(g_50,"50 #muBq/m^{3}");
  leg->AddEntry(g_150,"150 #muBq/m^{3}");
  leg->AddEntry(g_500,"500 #muBq/m^{3}");
  leg->SetFillColor(kWhite);

  g_50->SetFillColor(kWhite);
  g_150->SetFillColor(kWhite);
  g_500->SetFillColor(kWhite);

  g_50->Draw("");
  g_150->Draw("same");
  g_500->Draw("same");
  leg->Draw("same");

  return;
}
