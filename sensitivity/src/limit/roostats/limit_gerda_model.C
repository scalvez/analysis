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

void limit_gerda_model() {

  double lim = 0;
  for(unsigned int i=1; i<=100;++i) {
    // if(!(i%10))
    //   std::cout << " i " << i << std::endl;
  RooWorkspace w("w");
  w.factory("Uniform:bkg_pdf(x[0,100])");
  w.factory("Gaussian:sig_pdf(x, mass[50], sigma[1])");

  w.factory("SUM:model(nsig[0,100]*sig_pdf, nbkg[0,1000]*bkg_pdf)");  // for extended model

  RooAbsPdf * pdf = w.pdf("model");
  RooRealVar * x = w.var("x");  // the observable

  // generate the data (nsig = 30, nbkg=1000)
  w.var("nsig")->setVal(0);
  w.var("nbkg")->setVal(10);
  // use fixed random numbers for reproducibility
  // RooRandom::randomGenerator()->SetSeed(111);
  RooRandom::randomGenerator()->SetSeed(i);
  RooDataSet * data = pdf->generate(*x);  // will generate according to total S+B events
  data->SetName("data");
  w.import(*data);

  data->Print();
  x->setBins(100);
  RooPlot * plot = x->frame();
  data->plotOn(plot);
  // plot->Draw();

  pdf->fitTo(*data);

  // RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  // r->Print();

  pdf->plotOn(plot);
  pdf->plotOn(plot, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed) );
  pdf->plotOn(plot, RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  // add also the fit parameters
  pdf->paramOn(plot);
  // plot->Draw();

  RooStats::ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*pdf);
  mc.SetParametersOfInterest(*w.var("nsig"));
  mc.SetObservables(*w.var("x"));
  // define set of nuisance parameters
  // w.defineSet("nuisParams","a,nbkg");
  w.defineSet("nuisParams","nbkg");

  mc.SetNuisanceParameters(*w.set("nuisParams"));
  // import model in the workspace
  w.import(mc);

  double confidenceLevel = 0.9;

  /*
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

  // LikelihoodIntervalPlot * plot = new LikelihoodIntervalPlot(interval);
  // //plot->SetRange(0,100);  // possible eventually to change ranges
  // plot->SetNPoints(50);  // do not use too many points, it could become very slow for some models
  // plot->Draw("");  // use option TF1 if too slow (plot.Draw("tf1")
  */

  // Feldman-Cousing didn't appear to work

  /*
  //Bayesian MCMC
  // this proposal function seems fairly robust
  SequentialProposal sp(0.1);

  MCMCCalculator mcmc(*data,mc);
  mcmc.SetConfidenceLevel(confidenceLevel);
  //  mcmc.SetProposalFunction(*pf);
  mcmc.SetProposalFunction(sp);
  mcmc.SetNumIters(100000);         // Metropolis-Hastings algorithm iterations
  mcmc.SetNumBurnInSteps(50);       // first N steps to be ignored as burn-in

  // default is the shortest interval.  here use central
  // mcmc.SetLeftSideTailFraction(0.5); // for central Bayesian interval
  mcmc.SetLeftSideTailFraction(0); // for one-sided Bayesian interval

  RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();
  //firstPOI->setMax(100);

  MCMCInterval* interval = mcmc.GetInterval();

  // print out the iterval on the first Parameter of Interest
  cout << "\n90% interval on " <<firstPOI->GetName()<<" is : ["<<
    interval->LowerLimit(*firstPOI) << ", "<<
    interval->UpperLimit(*firstPOI) <<"] "<<endl;

  // // make a plot of posterior function
  // TCanvas* c1 = new TCanvas("IntervalPlot");
  // MCMCIntervalPlot plot2(*interval);
  // plot2.Draw();
  */


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

  // lim += interval->UpperLimit(*firstPOI);
  lim += interval->UpperLimit();
  }
  std::cout << std::endl << "--- Average limit is " << lim/100 << std::endl;
  return;
}
