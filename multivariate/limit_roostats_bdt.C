#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TFile.h"
#include "RooProfileLL.h"
#include "RooAbsPdf.h"
#include "RooStats/HypoTestResult.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
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

#include "config_sensitivity.h"

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;

void limit_roostats_bdt() {

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TFile *infile = TFile::Open("./bdt_scores/source/A/bdt_scores.root");
  TH1D* sig = (TH1D*)infile->Get("0nu");
  TH1D* bkg1 = (TH1D*)infile->Get("2nu");
  TH1D* bkg2 = (TH1D*)infile->Get("tl208");
  TH1D* bkg3 = (TH1D*)infile->Get("bi214");
  TH1D* bkg4 = (TH1D*)infile->Get("radon");

  // TFile *infile_sig = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/0nu_pdf_trunc.root");
  // TH1D* sig = (TH1D*)infile_sig->Get("2e_electrons_energy_sum");

  // TFile *infile_bkg1 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/2nu_pdf_trunc.root");
  // TH1D* bkg1 = (TH1D*)infile_bkg1->Get("2e_electrons_energy_sum");

  // TFile *infile_bkg2 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/tl208_pdf_trunc.root");
  // TH1D* bkg2 = (TH1D*)infile_bkg2->Get("2e_electrons_energy_sum");

  // TFile *infile_bkg3 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/bi214_pdf_trunc.root");
  // TH1D* bkg3 = (TH1D*)infile_bkg3->Get("2e_electrons_energy_sum");

  // TFile *infile_bkg4 = TFile::Open("$SW_WORK_DIR/analysis/sensitivity/data/pdf/radon_pdf_trunc.root");
  // TH1D* bkg4 = (TH1D*)infile_bkg4->Get("2e_electrons_energy_sum");

  // TString calculator = "simple_fc";
  // TString calculator = "bayesian";
  // TString calculator = "plc";
  TString calculator = "mcmc";
  // TString calculator = "fc";
  // TString calculator = "simple_fc_bis";
  // TString calculator = "my_fc";

  std::vector<double> limit_vector;
  std::vector<double> nwindow;

  double lim = 0;
  double count_lim = 0;

  double n_pseudo = 100
  for(auto i = 1; i<=n_pseudo; ++i) {
    if(i%1==0)
      std::cout << std::endl << " -------Pseudo experiment n° " << i << std::endl;

    // Recreate all variables, model and workspace, otherwise initialization issue
    RooRealVar score("score","score",-1,1);
    score.setBins(100);

    RooDataHist sig_hist("sig_hist","0nu histogram ",RooArgList(score),sig);
    RooDataHist bkg1_hist("bkg1_hist","2nu histogram ",RooArgList(score),bkg1);
    RooDataHist bkg2_hist("bkg2_hist","tl208 histogram ",RooArgList(score),bkg2);
    RooDataHist bkg3_hist("bkg3_hist","bi214 histogram ",RooArgList(score),bkg3);
    RooDataHist bkg4_hist("bkg4_hist","radon histogram ",RooArgList(score),bkg4);

    RooHistPdf hist_pdf_sig("hist_pdf_sig","hist_pdf_sig",score,sig_hist,0);
    RooHistPdf hist_pdf_bkg1("hist_pdf_bkg1","hist_pdf_bkg1",score,bkg1_hist,0);
    RooHistPdf hist_pdf_bkg2("hist_pdf_bkg2","hist_pdf_bkg2",score,bkg2_hist,0);
    RooHistPdf hist_pdf_bkg3("hist_pdf_bkg3","hist_pdf_bkg3",score,bkg3_hist,0);
    RooHistPdf hist_pdf_bkg4("hist_pdf_bkg4","hist_pdf_bkg4",score,bkg4_hist,0);

    // RooRealVar nsig("nsig","fitted number of 0nu events",0,0,50);
    // RooRealVar nbkg1("nbkg1","fitted number of 2nu events",34,0,200);
    // RooRealVar nbkg2("nbkg2","fitted number of tl208 events",0.06,0,200);
    // RooRealVar nbkg3("nbkg3","fitted number of bi214 events",0.34,0,200);
    // RooRealVar nbkg4("nbkg4","fitted number of radon events",0.64,0,200);

    RooRealVar nsig("nsig","fitted number of 0nu events",0,0,50);
    RooRealVar nbkg1("nbkg1","fitted number of 2nu events",conf_sens::N_source_2nu,0,2*conf_sens::N_source_2nu);
    RooRealVar nbkg2("nbkg2","fitted number of tl208 events",conf_sens::N_source_tl208,0,2*conf_sens::N_source_tl208);
    RooRealVar nbkg3("nbkg3","fitted number of bi214 events",conf_sens::N_source_bi214,0,2*conf_sens::N_source_bi214);
    RooRealVar nbkg4("nbkg4","fitted number of radon events",conf_sens::N_source_radon,0,2*conf_sens::N_source_radon);

    RooAddPdf model("model","model",RooArgList(hist_pdf_sig,hist_pdf_bkg1,hist_pdf_bkg2,hist_pdf_bkg3,hist_pdf_bkg4),RooArgList(nsig,nbkg1,nbkg2,nbkg3,nbkg4));

    RooRandom::randomGenerator()->SetSeed(i);
    // RooRandom::randomGenerator()->SetSeed(50);
    RooDataSet* data = model.generate(score);

    RooPlot* frame = score.frame(Title("BDT score"),Bins(100));

    // sig_hist.plotOn(frame);
    // bkg1_hist.plotOn(frame);
    data->plotOn(frame);
    model.plotOn(frame, RooFit::LineColor(kGray));
    model.plotOn(frame, RooFit::Components("hist_pdf_sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    model.plotOn(frame, RooFit::Components("hist_pdf_bkg1"), RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed) );
    model.plotOn(frame, RooFit::Components("hist_pdf_bkg2"), RooFit::LineColor(kGreen+1), RooFit::LineStyle(kDashed) );
    model.plotOn(frame, RooFit::Components("hist_pdf_bkg3"), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed) );
    model.plotOn(frame, RooFit::Components("hist_pdf_bkg4"), RooFit::LineColor(kMagenta), RooFit::LineStyle(kDashed) );

    RooFitResult * fit_res = model.fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit","Migrad"));

    RooWorkspace *w = new RooWorkspace("w","workspace") ;
    w->import(model) ;
    w->import(*data) ;
    // w->Print() ;

    RooStats::ModelConfig mc("ModelConfig",w);
    mc.SetPdf(model);
    mc.SetParametersOfInterest(*w->var("nsig"));
    mc.SetObservables(*w->var("score"));
    RooArgSet nuisance_pars(nbkg1,nbkg2,nbkg3,nbkg4);
    w->defineSet("nuisParams",nuisance_pars);

    mc.SetNuisanceParameters(*w->set("nuisParams"));

    w->import(mc);

    std::cout << " debug " << std::endl;
    double confidenceLevel = 0.9;

    double lowerLimit = 0;
    double upperLimit = 0;
    if(calculator.EqualTo("simple_fc")) {
      RooDataHist *h_data = data->binnedClone();

      // double best_i = -1;
      // double best_upper_limit = 99;
      // for(double i = 0.3; i<=1;i+=0.7) {
      //   TString cut = "score>"+ std::to_string(i);
      //   double n_window = h_data->sumEntries(cut);
      //   double n_exp = h_data->sumEntries(cut);
      //   TFeldmanCousins *simple_fc = new TFeldmanCousins(0.9);
      //   upperLimit = simple_fc->CalculateUpperLimit(n_window,n_exp);
      //   if(upperLimit<best_upper_limit) {
      //     best_upper_limit = upperLimit;
      //     best_i = i;
      //   }
      // }
      // upperLimit=best_upper_limit;

      double n_window = h_data->sumEntries("score>0.4");
      double n_exp = h_data->sumEntries("score>0.4");
      TFeldmanCousins *simple_fc = new TFeldmanCousins(0.9);
      lowerLimit = simple_fc->CalculateLowerLimit(n_window,n_exp);
      upperLimit = simple_fc->CalculateUpperLimit(n_window,n_exp);

    }
    else if(calculator.EqualTo("bayesian")) {
      // Bayesian calculator
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

      lowerLimit = interval->LowerLimit();
      upperLimit = interval->UpperLimit();

      cout << "\n90% interval on " <<firstPOI->GetName()<<" is : ["<<
        lowerLimit << ", "<<
        upperLimit <<"] "<<endl;

      if(i==1) {
        // draw plot of posterior function

        TCanvas* c_bc = new TCanvas("IntervalPlot");

        RooPlot * plot_bc = bayesianCalc.GetPosteriorPlot();
        if (plot_bc) plot_bc->Draw();
      }

    }        //  --- End Bayesian Calculator
    else if(calculator.EqualTo("plc")) {
      // ---- Profile Likelihood Calculator
      ProfileLikelihoodCalculator pl(*data,mc);
      pl.SetConfidenceLevel(confidenceLevel); // 68.3% interval
      LikelihoodInterval* interval = pl.GetInterval();

      // find the interval on the first Parameter of Interest (nsig)
      RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

      lowerLimit = interval->LowerLimit(*firstPOI);
      upperLimit = interval->UpperLimit(*firstPOI);

      cout << "\n90% interval on " <<firstPOI->GetName()<<" is : ["<<
        lowerLimit << ", "<<
        upperLimit <<"] "<<endl;

    } //  --- End Profile Likelihood Calculator
    else if (calculator.EqualTo("mcmc")) {
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

      upperLimit = interval->UpperLimit(*firstPOI);
    }// End MCMC
    else if (calculator.EqualTo("fc")) {
      // Feldman-Cousins

      FeldmanCousins fc(*data, mc);
      fc.SetConfidenceLevel( confidenceLevel);
      fc.SetNBins(100); // number of points to test per parameter
      fc.UseAdaptiveSampling(true); // make it go faster

      // Here, we consider only ensembles with 100 events
      // The PDF could be extended and this could be removed
      fc.FluctuateNumDataEntries(false);

      PointSetInterval* interval = (PointSetInterval*) fc.GetInterval();
      RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

      std::cout << "fc interval is ["<<
        interval->LowerLimit(*firstPOI) << " , "  <<
        interval->UpperLimit(*firstPOI) << "]" << endl;

      upperLimit = interval->UpperLimit(*firstPOI);

    }    // End of Feldman-Cousins
    else if(calculator.EqualTo("simple_fc_bis")) {
      RooDataHist *h_data = data->binnedClone();
      double n_window = h_data->sumEntries("score>0.4");
      // double n_exp = 0;
      double n_exp = h_data->sumEntries("score>0.4");
      TFeldmanCousins *simple_fc = new TFeldmanCousins(0.9);
      lowerLimit = simple_fc->CalculateLowerLimit(n_window,n_exp);
      upperLimit = simple_fc->CalculateUpperLimit(n_window,n_exp);
      nwindow.push_back(n_window);
    }  // End of Simple FC bis
    else {
      std::cout << "No Calculator specified " << std::endl;
      return;
    }

    if(upperLimit == upperLimit) {
      count_lim++;
      lim += upperLimit;
      limit_vector.push_back(upperLimit);
    }

    frame->Draw();
  }

  std::cout << " average limit is " << lim/count_lim << std::endl;
  std::cout << " count limit " << count_lim << std::endl;

  std::sort(limit_vector.begin(),limit_vector.end());
  std::cout << " median is " << limit_vector[int(limit_vector.size()/2)] << std::endl;

  // for(auto i=0; i<limit_vector.size();++i)
  //   std::cout << limit_vector[i] << std::endl;
  // for(auto i=0; i<nwindow.size();++i)
  //   std::cout << nwindow[i] << std::endl;
  return;
}
