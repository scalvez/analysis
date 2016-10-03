#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>

#include "sensitivity_constants.h"
#include "analysis_config.h"
// #include "sensitivity_measurements.h"

extern std::map < TString , double > quantity_efficiency;
extern std::map < TString , TH1F* > quantity_pdf;

void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity, double seed) {

  std::cout << " ---- Isotope : " << isotope << std::endl;

  TString input_file = "../pdf/" + isotope + "_pdf.root";
  TString output_file = "../pseudo/" + isotope + "_pseudo.root";

  TFile *f = TFile::Open(input_file);
  TFile *f_output= new TFile(output_file, "RECREATE");

  TRandom *rdm = new TRandom();
  rdm->SetSeed(seed);

  int n_decays = 0;
  if (isotope.Contains("2nu"))
    n_decays = int(mass*exposure_y/halflife_2nu*const_se);
  else if (isotope.Contains("radon"))
    n_decays = int(activity*exposure_sec*tracker_volume);
  else
    n_decays = int(activity*exposure_sec*mass);

  int n_decays_rdm = 0;

  if(poisson_pseudo)
    n_decays_rdm = rdm->Poisson(n_decays);
  else
    n_decays_rdm = n_decays;

  for(auto it = quantities.begin(); it != quantities.end(); ++it) {
    TString qty = *it;
    TH1 *h = (TH1*)f->Get(qty);
    TH1 *h_cdf = h->GetCumulative();
    TString key = isotope + "_" + qty;

    std::cout << " ---------- Qty : " << qty << std::endl;

    int n_events = int(n_decays_rdm*quantity_efficiency.at(key));

    // std::cout << " Generating pseudo-experiment with " << n_decays << " decays" << std::endl;
    // std::cout << " Generating pseudo-experiment with " << n_decays_rdm << " rdm decays" << std::endl;
    std::cout << " Generating pseudo-experiment with " << n_events << " events" << "(eff=" << quantity_efficiency.at(key) << ")" << std::endl;
    int nbins = h_cdf->GetNbinsX();

    TH1F *h_pseudo = new TH1F(qty,qty,nbins,h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

    for(int i=0; i<n_events;++i) {
      double rdm_number = rdm->Uniform(1);
      for(int j=1; j<=nbins;++j) {
        if(h_cdf->GetBinContent(j)>rdm_number) {
          h_pseudo->Fill(h_cdf->GetXaxis()->GetBinCenter(j));
          break;
        }
      }
    }

    TString pseudo_quantity = "pseudo_" + qty;
    quantity_pdf.insert(std::pair<TString,TH1F*>(pseudo_quantity,h_pseudo));

    h_pseudo->Write();
  }
  f->Close();
  f_output->Close();
  return;
}
