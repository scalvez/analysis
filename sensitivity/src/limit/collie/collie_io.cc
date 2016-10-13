#include "CollieIOFile.hh"
#include "TRandom.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////
  ///Create IO file with input parameters
  /////////////////////////////////////////
  CollieIOFile* cfile = new CollieIOFile();
  // Specify outputfile and channel name
  cfile->initFile("exampleCollieIOfile.root", "My Channel");

  // Define your input histograms
  double Xmin = 0.0;
  double Xmax = 100.0;
  int Nbins = 50;
  //

  cfile->setInputHist(Xmin,Xmax,Nbins);

  // Option to define physical cutoffs where events should not exist (in terms of your historam range)
  //
  //  cfile->setCutoffs(0.0,1.0);

  cfile->setRebin(1);

  // cfile->setSmooth(false);
  //  cfile->setHistNorm(50000);


  //Define backgrounds
  vector<string> bkgdNames;
  bkgdNames.push_back("Bkgd1");
  cfile->createChannel(bkgdNames);


  // For this example, generate backgrounds, signal, data, and systematics
  TH1D* bkgd1 = new TH1D("bkgd1","bkgd1",Nbins,Xmin,Xmax);

  TH1D* sig = new TH1D("sig","sig",Nbins,Xmin,Xmax);
  TH1D* data = new TH1D("data","data",Nbins,Xmin,Xmax);

  // Make sure you keep track of statistical uncertainties in histograms correctly
  bkgd1->Sumw2();
  sig->Sumw2();
  data->Sumw2();


  //Otherwise, get your input histograms from an external file
  //
  //  TFile infile("myInputFile.root");
  //  TH1D* data = (TH1D*)infile.Get("data");
  //  TH1D* signal = (TH1D*)infile.Get("signal");
  //  TH1D* bkgd1 = (TH1D*)infile.Get("bkgd1");

  TRandom r(1234);
  double niter = 5e5;

  //We'll make three mass points...
  for(int m=75; m<=125; m+=1){
  // for(int m=100; m<=100; m+=25){

    bkgd1->Scale(0);
    sig->Scale(0);
    data->Scale(0);

    // For this example, randomly fill histograms
    for(int i=0; i<niter; i++){
      //Flat background
      bkgd1->Fill(r.Rndm()*100,1000.0/niter);

      //Fill data with 5% higher bkgd1 rate
      data->Fill(r.Rndm()*100, 1000.0/niter);

      //Peaked signal
      sig->Fill(r.Gaus(50,2),100/niter);
    }

    //Backgrounds are passed in via vector
    vector<TH1D*> vbkgd;
    vbkgd.push_back(bkgd1);

    //Alpha parameters only matter when smoothing is utilized
    //  Input values don't matter if you're not smoothing.
    //  Don't smooth unless you know what you're doing.
    vector<double> valpha;
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
    cfile->createFlatSigSystematic("Lumi",0.1,0.1,m);

    cfile->createFlatBkgdSystematic(0,"Lumi",0.1,0.1,m);

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
    //
    //  cfile->setBkgdFloatFlag(0,"Eff",true,m);
    //  cfile->setBkgdFloatFlag(1,"Eff",true,m);
    //  cfile->setSigFloatFlag("Eff",true,m);


    //  ==>For large uncertainties (eg, >30%) use a log-normal PDF to avoid
    //     problems in Gaussian PDF modeling.
    //
    //  cfile->setLogNormalFlag("Eff",true,m);

  }
  ///store and output channel information
  cfile->storeFile();
}
