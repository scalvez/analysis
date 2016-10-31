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

  w.factory("SUM:model(nsig[0,1000]*sig_pdf, nbkg[0,1000000]*bkg_pdf)");  // for extended model

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
  plot->Draw();

  pdf->fitTo(*data);

  // RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  // r->Print();

  pdf->plotOn(plot);
  pdf->plotOn(plot, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed) );
  pdf->plotOn(plot, RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  // add also the fit parameters
  pdf->paramOn(plot);
  plot->Draw();

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

  // // write the workspace in the file
  // TString fileName = "SPlusBExpoModel.root";
  // w.writeToFile(fileName,true);
  // cout << "model written to file " << fileName << endl;

  ProfileLikelihoodCalculator pl(*data,mc);
  pl.SetConfidenceLevel(0.90); // 68.3% interval
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
  lim += upperLimit;
  }
  std::cout << std::endl << "--- Average limit is " << lim/100 << std::endl;
  return;
}
