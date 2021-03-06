#include "CollieIOFile.hh"
#include "TRandom.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////
  ///Create IO file with input parameters
  /////////////////////////////////////////
  CollieIOFile* cfile = new CollieIOFile();
  // Specify outputfile and channel name
  cfile->initFile("./collie_io/spectra/io_file_spectra_CD_Ecut.root", "0nu");

  // Define your input histograms
  double Xmin = 2;
  double Xmax = 5;
  int Nbins = 60;

  cfile->setInputHist(Xmin,Xmax,Nbins);

  // cfile->setInputHist2D(Xmin,Xmax,Nbins);

  // Option to define physical cutoffs where events should not exist (in terms of your historam range)
  //
  //  cfile->setCutoffs(0.0,1.0);

  cfile->setRebin(1);

  // cfile->setSmooth(false);
  // cfile->setHistNorm(50000);

  //Define backgrounds
  vector<string> bkgdNames;
  bkgdNames.push_back("2nu");
  bkgdNames.push_back("tl208");
  bkgdNames.push_back("bi214");
  bkgdNames.push_back("radon");
  cfile->createChannel(bkgdNames);

  TFile infile("./spectra/source/CD/spectra_counts_Ecut.root");
  TH1D* sig = (TH1D*)infile.Get("0nu");
  TH1D* bkgd1 = (TH1D*)infile.Get("2nu");
  TH1D* bkgd2 = (TH1D*)infile.Get("tl208");
  TH1D* bkgd3 = (TH1D*)infile.Get("bi214");
  TH1D* bkgd4 = (TH1D*)infile.Get("radon");
  TH1D* data = (TH1D*)infile.Get("h_data");

  // Make sure you keep track of statistical uncertainties in histograms correctly
  // bkgd1->Sumw2();
  // bkgd2->Sumw2();
  // sig->Sumw2();
  // data->Sumw2();

  TRandom r(1234);
  double niter = 5e5;

  //We'll make three mass points...
  // for(int m=0; m<=10; m+=1){
  for(int m=0; m<=0; m+=1){

    //full range
    // bkgd1->Scale(7655);
    // bkgd2->Scale(1);
    // bkgd3->Scale(7);
    // bkgd4->Scale(15);
    // sig->Scale(1);

    // sig->Scale(1/sig->Integral(1,sig->GetNbinsX()));
    // bkgd1->Scale(34.3/bkgd1->Integral(1,bkgd1->GetNbinsX()));
    // bkgd2->Scale(0.06/bkgd2->Integral(1,bkgd2->GetNbinsX()));
    // bkgd3->Scale(0.34/bkgd3->Integral(1,bkgd3->GetNbinsX()));
    // bkgd4->Scale(0.64/bkgd4->Integral(1,bkgd4->GetNbinsX()));

    //Backgrounds are passed in via vector
    vector<TH1D*> vbkgd;
    vbkgd.push_back(bkgd1);
    vbkgd.push_back(bkgd2);
    vbkgd.push_back(bkgd3);
    vbkgd.push_back(bkgd4);

    //Alpha parameters only matter when smoothing is utilized
    //  Input values don't matter if you're not smoothing.
    //  Don't smooth unless you know what you're doing.
    vector<double> valpha;
    valpha.push_back(-1);
    valpha.push_back(-1);
    valpha.push_back(-1);
    valpha.push_back(-1);

    ///Use this tool to allow collie to generate a low-stats safe
    // binning for your histograms.  Histogram binning should be
    // larger than your desired final number of bins.  For more
    // details see CollieIOfile.h
    // Should be done for each mass point in the file, but with
    // the same number of output bins (specified above).  IE, if
    // you have a mass point loop, put this inside the loop.
    //
    //         TH1D* btotal = (TH1D*)bkgd1->Clone("btotal");
    //         btotal->Add(bkgd2);
    //         cfile->generateBinMap(btotal,sig,"MVA");
    //         delete btotal; btotal = NULL;


    //Each parameter point has a signal histo, data histo, and an array of backgrounds...
    //  Smoothing parameters are also passed in.
    cfile->createMassPoint(m, data, sig, -1, vbkgd,valpha);


    // If you have more than one mass point, you may choose to interpolate on some parameter grid
    //cfile->interpolateMassGrid(5,100,110);


    // Add systematics...either flat or by shape (ie, function of final variable)
    //   if by shape, must supply a histogram of the values in percent(%) fluctuations...
    //   Signal requires no index, but backgrounds must be specifically indexed (0->N bkgds)
    //   Read the instructions in collie/io/include/CollieIOFile.hh if you're in doubt
    // cfile->createFlatSigSystematic("Lumi",0.01,0.01,m);

    double eff_error = 0.05;

    //test
    cfile->createFlatSigSystematic("Eff",eff_error,eff_error,m);

    // // cfile->createFlatBkgdSystematic(0,"Lumi",0.01,0.01,m);
    cfile->createFlatBkgdSystematic(0,"Eff",eff_error,eff_error,m);
    cfile->createFlatBkgdSystematic(1,"Eff",eff_error,eff_error,m);
    cfile->createFlatBkgdSystematic(2,"Eff",eff_error,eff_error,m);
    cfile->createFlatBkgdSystematic(3,"Eff",eff_error,eff_error,m);

    // Example of systematics input as histograms, can be flat or function of final variable
    //==>Use this method if you're inputing fractional shape systematics
    //    TH1D* systP = (TH1D*)infile.Get("signal_Systematic_positive");
    //    TH1D* systN = (TH1D*)infile.Get("signal_Systematic_negative");
    //    cfile->createSigSystematic("ShapeSyst",systP,systN,m);

    //==>Use this method if you're inputing a different shape template
    //    systP = (TH1D*)infile.Get("bkgd_BkgdShape_positive");
    //    systN = (TH1D*)infile.Get("bkgd_BkgdShape_negative");
    //    cfile->createShapeBkgdSystematic(0,"BkgdShape",systP,systN,m);
    //    cfile->createShapeBkgdSystematic(1,"BkgdShape",systP,systN,m);


    //  ==>Option to remove prior constraint on systematic uncertainty PDF.
    //     Floating makes a parameter a free parameter in the fit.

    // The systematics are still taken into account
    cfile->setBkgdFloatFlag(0,"Eff",true,m);
    cfile->setBkgdFloatFlag(1,"Eff",true,m);
    cfile->setBkgdFloatFlag(2,"Eff",true,m);
    cfile->setBkgdFloatFlag(3,"Eff",true,m);

    cfile->setSigFloatFlag("Eff",true,m);


    //  ==>For large uncertainties (eg, >30%) use a log-normal PDF to avoid
    //     problems in Gaussian PDF modeling.
    //
    //  cfile->setLogNormalFlag("Eff",true,m);

  }
  ///store and output channel information
  cfile->storeFile();
}
