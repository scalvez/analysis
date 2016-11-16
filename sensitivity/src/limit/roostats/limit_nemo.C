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

  RooRealVar nsig("nsig","fitted number of 0nu events",0,0,1000);
  RooRealVar nbkg1("nbkg","fitted number of 2nu events",34,0,10000);

  RooAddPdf model("model","model",RooArgList(hist_pdf_sig,hist_pdf_bkg1),RooArgList(nsig,nbkg1));

  // RooRandom::randomGenerator()->SetSeed(4);
  RooDataSet* data = model.generate(E);

  RooPlot* frame = E.frame(Title("Energy sum"),Bins(15));

  // sig_hist.plotOn(frame);
  // bkg1_hist.plotOn(frame);
  data->plotOn(frame);
  data->plotOn(frame);
  model.plotOn(frame, RooFit::LineColor(kGray));
  model.plotOn(frame, RooFit::Components("hist_pdf_bkg1"), RooFit::LineStyle(kDashed) );
  model.plotOn(frame, RooFit::Components("hist_pdf_sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  RooFitResult * fit_res = model.fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));

  // fit_res->Print();

  // std::cout << " test " << nsig << std::endl;
  // model.plotOn(frame);
  // frame->Draw();

  RooWorkspace *w = new RooWorkspace("w","workspace") ;
  w->import(model) ;
  w->import(*data) ;
  w->Print() ;

  // RooRealVar * test = w->var("nsig");
  // test->plotOn(frame);
  // std::cout << "test " << test << std::endl;

  // RooAbsPdf * pdf = w.pdf("model");

  RooStats::ModelConfig mc("ModelConfig",w);
  mc.SetPdf(model);
  mc.SetParametersOfInterest(*w->var("nsig"));
  mc.SetObservables(*w->var("E"));
  w->defineSet("nuisParams","nbkg");

  mc.SetNuisanceParameters(*w->set("nuisParams"));

  w->import(mc);

  double confidenceLevel = 0.9;

  /* // Bayesian calculator
    BayesianCalculator bayesianCalc(*data,mc);
  bayesianCalc.SetConfidenceLevel(confidenceLevel);

  // set the type of interval (not really needed for central which is the default)
  // bayesianCalc.SetLeftSideTailFraction(0.5); // for central interval
  bayesianCalc.SetLeftSideTailFraction(0.); // for upper limit
  //bayesianCalc.SetShortestInterval(); // for shortest interval


  // set the integration type (not really needed for the default ADAPTIVE)
  // possible alternative values are  "VEGAS" , "MISER", or "PLAIN"  (MC integration from libMathMore)
  // "TOYMC" (toy MC integration, work when nuisances exist and they have a constraints pdf)
  TString integrationType = "";

  // this is needed if using TOYMC
  if (integrationType.Contains("TOYMC") ) {
    RooAbsPdf * nuisPdf = RooStats::MakeNuisancePdf(mc, "nuisance_pdf");
    if (nuisPdf) bayesianCalc.ForceNuisancePdf(*nuisPdf);
  }

  bayesianCalc.SetIntegrationType(integrationType);

  // compute interval by scanning the posterior function
  // it is done by default when computing shortest intervals
  bayesianCalc.SetScanOfPosterior(100);

  RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

  SimpleInterval* interval = bayesianCalc.GetInterval();
  if (!interval) {
     cout << "Error computing Bayesian interval - exit " << endl;
     return;
  }


  double lowerLimit = interval->LowerLimit();
  double upperLimit = interval->UpperLimit();

  cout << "\n90% interval on " <<firstPOI->GetName()<<" is : ["<<
    lowerLimit << ", "<<
    upperLimit <<"] "<<endl;

  // draw plot of posterior function

  TCanvas* c_bc = new TCanvas("IntervalPlot");

  RooPlot * plot_bc = bayesianCalc.GetPosteriorPlot();
  if (plot_bc) plot_bc->Draw();

  // // lim += interval->UpperLimit(*firstPOI);
  // //lim += interval->UpperLimit();
  // lim += upperLimit;

  */  //  --- End Bayesian Calculator

  // ---- Profile Likelihood Calculator
  ProfileLikelihoodCalculator pl(*data,mc);
  pl.SetConfidenceLevel(confidenceLevel); // 68.3% interval
  LikelihoodInterval* interval = pl.GetInterval();

  // find the interval on the first Parameter of Interest (nsig)
  RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

  double lowerLimit = interval->LowerLimit(*firstPOI);
  double upperLimit = interval->UpperLimit(*firstPOI);

  cout << "\n90% interval on " <<firstPOI->GetName()<<" is : ["<<
    lowerLimit << ", "<<
    upperLimit <<"] "<<endl;
  //  --- End Profile Likelihood Calculator

  frame->Draw();

  return;
}
