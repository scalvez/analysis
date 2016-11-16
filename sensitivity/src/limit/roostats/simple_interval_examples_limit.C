// Example showing confidence intervals with four techniques.
/*
IntervalExamples

Author Kyle Cranmer
date   Sep. 2010

An example that shows confidence intervals with four techniques.
The model is a Normal Gaussian G(x|mu,sigma) with 100 samples of x.
The answer is known analytically, so this is a good example to validate
the RooStats tools.

expected interval is [-0.162917, 0.229075]
plc  interval is     [-0.162917, 0.229075]
fc   interval is     [-0.17    , 0.23]        // stepsize is 0.01
bc   interval is     [-0.162918, 0.229076]
mcmc interval is     [-0.166999, 0.230224]


*/

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"

#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooPoisson.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include"Math/DistFunc.h"
#include "TH1F.h"
#include "TMarker.h"
#include "TStopwatch.h"

#include <iostream>

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;

void simple_interval_examples_limit()
{

  // Time this macro
  TStopwatch t;
  t.Start();


  // set RooFit random seed for reproducible results
  RooRandom::randomGenerator()->SetSeed(30011);

  // make a simple model via the workspace factory
  RooWorkspace* wspace = new RooWorkspace();
  wspace->factory("Gaussian::signal(x[0,100],mu[50],sigma[1])");
  wspace->factory("Uniform::background(x[0,100])");
  wspace->defineSet("poi","mu");
  wspace->defineSet("obs","x");

  // specify components of model for statistical tools
  ModelConfig* modelConfig = new ModelConfig("Example G(x|mu,1)");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf( *wspace->pdf("normal") );
  modelConfig->SetParametersOfInterest( *wspace->set("poi") );
  modelConfig->SetObservables( *wspace->set("obs") );

  // create a toy dataset
  RooDataSet* data = wspace->pdf("normal")->generate(*wspace->set("obs"),100);
  std::cout << "--- Printing Data " << std::endl;
  data->Print();
  std::cout << "--- Data printed " << std::endl;


  // for convenience later on
  RooRealVar* x = wspace->var("x");
  RooRealVar* mu = wspace->var("mu");

  double tmp_mean_mu = data->mean(*x);
  std::cout << "--- Mu : " <<  data->mean(*x) << std::endl;


  // some plots
  TCanvas* canvas = new TCanvas("canvas");
  canvas->Divide(2,2);

  // plot the data
  canvas->cd(1);
  RooPlot* frame = x->frame();
  data->plotOn(frame);
  data->statOn(frame);
  frame->Draw();

  // set confidence level
  // double confidenceLevel = 0.95;
  double confidenceLevel = 0.9;

  // example use profile likelihood calculator
  ProfileLikelihoodCalculator plc(*data, *modelConfig);
  plc.SetConfidenceLevel( confidenceLevel);
  LikelihoodInterval* plInt = plc.GetInterval();
  std::cout << "--- after plc Mu : " <<  data->mean(*x) << std::endl;

  // example use of Feldman-Cousins
  FeldmanCousins fc(*data, *modelConfig);
  fc.SetConfidenceLevel( confidenceLevel);
  fc.SetNBins(100); // number of points to test per parameter
  fc.UseAdaptiveSampling(true); // make it go faster
  std::cout << "--- after fc Mu : " <<  data->mean(*x) << std::endl;

  // Here, we consider only ensembles with 100 events
  // The PDF could be extended and this could be removed
  fc.FluctuateNumDataEntries(false);
  std::cout << "--- after fc fluctuate : " <<  data->mean(*x) << std::endl;

  // Proof
  //  ProofConfig pc(*wspace, 4, "workers=4", kFALSE);    // proof-lite
  //ProofConfig pc(w, 8, "localhost");    // proof cluster at "localhost"
  //  ToyMCSampler* toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler();
  //  toymcsampler->SetProofConfig(&pc);     // enable proof

  PointSetInterval* interval = (PointSetInterval*) fc.GetInterval();
  std::cout << "--- after fc interval : " <<  data->mean(*x) << std::endl;

  // example use of BayesianCalculator
  // now we also need to specify a prior in the ModelConfig
  wspace->factory("Uniform::prior(mu)");
  modelConfig->SetPriorPdf(*wspace->pdf("prior"));

  // example usage of BayesianCalculator
  BayesianCalculator bc(*data, *modelConfig);
  bc.SetConfidenceLevel( confidenceLevel);
  SimpleInterval* bcInt = bc.GetInterval();
  std::cout << "--- after bc Mu : " <<  data->mean(*x) << std::endl;

  // example use of MCMCInterval
  MCMCCalculator mc(*data, *modelConfig);
  mc.SetConfidenceLevel( confidenceLevel);
  // special options
  mc.SetNumBins(200);        // bins used internally for representing posterior
  mc.SetNumBurnInSteps(500); // first N steps to be ignored as burn-in
  mc.SetNumIters(100000);    // how long to run chain
  mc.SetLeftSideTailFraction(0.5); // for central interval
  MCMCInterval* mcInt = mc.GetInterval();
  std::cout << "--- after mc Mu : " <<  data->mean(*x) << std::endl;

  // // for this example we know the expected intervals
  // double expectedLL = data->mean(*x)
  //   + ROOT::Math::normal_quantile(  (1-confidenceLevel)/2,1)
  //   / sqrt(data->numEntries());
  // double expectedUL = data->mean(*x)
  //   + ROOT::Math::normal_quantile_c((1-confidenceLevel)/2,1)
  //   / sqrt(data->numEntries()) ;
  // for this example we know the expected intervals
  double expectedLL = tmp_mean_mu
    + ROOT::Math::normal_quantile(  (1-confidenceLevel)/2,1)
    / sqrt(data->numEntries());
  double expectedUL = tmp_mean_mu
    + ROOT::Math::normal_quantile_c((1-confidenceLevel)/2,1)
    / sqrt(data->numEntries()) ;

  std::cout << "----------------" << std::endl;
  std::cout << "--- Mu : " <<  data->mean(*x) << std::endl;
  std::cout << "--- Mu saved : " <<  tmp_mean_mu << std::endl;
  // Use the intervals
  std::cout << "expected interval is [" <<
    expectedLL << ", " <<
    expectedUL << "]" << endl;

  cout << "plc interval is [" <<
    plInt->LowerLimit(*mu) << ", " <<
    plInt->UpperLimit(*mu) << "]" << endl;

  std::cout << "fc interval is ["<<
    interval->LowerLimit(*mu) << " , "  <<
    interval->UpperLimit(*mu) << "]" << endl;

  cout << "bc interval is [" <<
    bcInt->LowerLimit() << ", " <<
    bcInt->UpperLimit() << "]" << endl;

  cout << "mc interval is [" <<
    mcInt->LowerLimit(*mu) << ", " <<
    mcInt->UpperLimit(*mu) << "]" << endl;

  mu->setVal(10);
  cout << "is mu=10 in the interval? " <<
    plInt->IsInInterval(RooArgSet(*mu)) << endl;


  // make a reasonable style
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetStatColor(0);

  // plot the profile likeihood
  canvas->cd(2);
  LikelihoodIntervalPlot plot(plInt);
  plot.Draw();

  // plot the MCMC interval
  canvas->cd(3);
  MCMCIntervalPlot* mcPlot = new MCMCIntervalPlot(*mcInt);
  mcPlot->SetLineColor(kGreen);
  mcPlot->SetLineWidth(2);
  mcPlot->Draw();

  canvas->cd(4);
  RooPlot * bcPlot = bc.GetPosteriorPlot();
  bcPlot->Draw();

  canvas->Update();

  t.Stop();
  t.Print();

}
