#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TFile.h"
#include "RooProfileLL.h"
#include "RooAbsPdf.h"
#include "RooStats/HypoTestResult.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooTreeDataStore.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStopwatch.h"

#include "RooHistPdf.h"
#include "RooAddPdf.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/UniformProposal.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/NumberCountingPdfFactory.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalFunction.h"
#include "RooStats/ProposalHelper.h"
#include "RooFitResult.h"
#include "TGraph2D.h"

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;

void limit_nemo() {

  TFile *infile_sig = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/0nu_pdf_trunc.root");
  TH1D* sig = (TH1D*)infile_sig->Get("2e_electrons_energy_sum");

  TFile *infile_bkg1 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/2nu_pdf_trunc.root");
  TH1D* bkg1 = (TH1D*)infile_bkg1->Get("2e_electrons_energy_sum");

  // std::cout << " integral signal " <<  bkg1->Integral(1,15) << std::endl;
  // return;

  RooRealVar E("E","E",2.45,3.2);
  E.setBins(15);

  RooDataHist sig_hist("sig_hist","0nu histogram ",RooArgList(E),sig);
  RooDataHist bkg1_hist("bkg1_hist","2nu histogram ",RooArgList(E),bkg1);

  RooHistPdf hist_pdf_sig("hist_pdf_sig","hist_pdf_sig",E,sig_hist,0);
  RooHistPdf hist_pdf_bkg1("hist_pdf_bkg1","hist_pdf_bkg1",E,bkg1_hist,0);

  RooRealVar nsig("nsig","fitted number of 0nu events",100,0,1000);
  RooRealVar nbkg1("nbkg","fitted number of 2nu events",100,0,10000);

  RooAddPdf model("model","model",RooArgList(hist_pdf_sig,hist_pdf_bkg1),RooArgList(nsig,nbkg1));

  RooDataSet* data = model.generate(E);

  RooPlot* frame = E.frame(Title("Energy sum"),Bins(15));

  // sig_hist.plotOn(frame);
  // bkg1_hist.plotOn(frame);
  data->plotOn(frame);
  model.fitTo(*data);

  model.plotOn(frame);
  frame->Draw();

  return;
}
