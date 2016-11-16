#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

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

void test_hist_pdf() {

  RooWorkspace w("w");
  w.factory("Uniform:bkg_pdf(x[0,100])");
  w.factory("Gaussian:sig_pdf(x, mass[50], sigma[5])");

  w.factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[0,10000]*bkg_pdf)");  // for extended model

  RooAbsPdf * pdf = w.pdf("model");
  RooRealVar * x = w.var("x");  // the observable

  w.var("nsig")->setVal(10000);
  w.var("nbkg")->setVal(0);

  RooDataSet* data_sig_pdf = pdf->generate(*x) ;

  w.var("nsig")->setVal(0);
  w.var("nbkg")->setVal(10000);

  RooDataSet* data_bkg_pdf = pdf->generate(*x) ;

  RooPlot* frame = x->frame(Title("Test"),Bins(100)) ;

  // data_sig_pdf->plotOn(frame);
  // data_bkg_pdf->plotOn(frame);

  RooDataHist* hist_sig = data_sig_pdf->binnedClone() ;
  RooDataHist* hist_bkg = data_bkg_pdf->binnedClone() ;

  // hist_sig->plotOn(frame);
  // hist_bkg->plotOn(frame);

  RooHistPdf hist_pdf_sig("hist_pdf_sig","hist_pdf_sig",*x,*hist_sig,0) ;
  RooHistPdf hist_pdf_bkg("hist_pdf_bkg","hist_pdf_bkg",*x,*hist_bkg,0) ;

  // hist_pdf_sig.plotOn(frame);
  // hist_pdf_bkg.plotOn(frame);

  w.var("nsig")->setVal(1000);
  w.var("nbkg")->setVal(10000);

  RooDataSet* data = pdf->generate(*x) ;

  RooRealVar nsig("nsig","fitted number of signal events",1000, 0, 1000) ;
  RooRealVar nbkg("nbkg","fitted number of bkg events",500,0,10000);

  RooAddPdf model_hist("model_hist","model_hist",RooArgList(hist_pdf_sig,hist_pdf_bkg),RooArgList(nsig,nbkg)) ;

  model_hist.fitTo(*data);
  data->plotOn(frame);
  model_hist.plotOn(frame);

  frame->Draw();

  return;
}
