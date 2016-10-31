#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooProfileLL.h"
#include "RooAbsPdf.h"
#include "RooStats/HypoTestResult.h"
#include "RooRealVar.h"
#include "RooPlot.h"
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

using namespace RooFit;
using namespace RooStats;

void limit_simple_counting( int nobs = 3,           // number of observed events
                            double b = 1,           // number of background events
                            double sigmab = 0.2 )   // relative uncertainty in b
{
  nobs = 0;
  b = 0;
  sigmab = 0.1;

  RooWorkspace w("w");

// make Poisson model * Gaussian constraint
   w.factory("sum:nexp(s[3,0,15],b[1,0,10])");
// Poisson of (n | s+b)
   w.factory("Poisson:pdf(nobs[0,50],nexp)");
   w.factory("Gaussian:constraint(b0[0,10],b,sigmab[1])");
   w.factory("PROD:model(pdf,constraint)");


   w.var("b0")->setVal(b);
   w.var("b0")->setConstant(true); // needed for being treated as global observables
   w.var("sigmab")->setVal(sigmab*b);


   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("s"));
   mc.SetObservables(*w.var("nobs"));
   mc.SetNuisanceParameters(*w.var("b"));

   // these are needed for the hypothesis tests
   mc.SetSnapshot(*w.var("s"));
   mc.SetGlobalObservables(*w.var("b0"));

   mc.Print();
   // import model in the workspace
   w.import(mc);

   // make data set with the number of observed events
   RooDataSet data("data","", *w.var("nobs"));
   w.var("nobs")->setVal(3);
   data.add(*w.var("nobs") );
   // import data set in workspace and save it in a file
   w.import(data);

   w.Print();

   TString fileName = "CountingModel.root";

   // write workspace in the file (recreate file if already existing)
   w.writeToFile(fileName, true);

   ProfileLikelihoodCalculator plc(data, mc);
  plc.SetTestSize(.10);
  ConfInterval* lrint = plc.GetInterval();  // that was easy.

  // Let's make a plot
  TCanvas* dataCanvas = new TCanvas("dataCanvas");
  // dataCanvas->Divide(2,1);

  // dataCanvas->cd(1);
  LikelihoodIntervalPlot plotInt((LikelihoodInterval*)lrint);
  plotInt.SetTitle("Profile Likelihood Ratio and Posterior for S");
  plotInt.Draw();


}
