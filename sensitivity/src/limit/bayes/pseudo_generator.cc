#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>

void pseudo_generator(TH1 * h_pdf, TH1 * h_pseudo, int n_events, double seed) {

  TRandom *rdm = new TRandom();
  rdm->SetSeed(seed);

  int n_events_rdm = rdm->Poisson(n_events);

  TH1 *h_cdf = h_pdf->GetCumulative();

  std::cout << " Generating pseudo-experiment with " << n_events_rdm << " events" << std::endl;

  int nbins = h_cdf->GetNbinsX();

  // TH1F *h_pseudo = new TH1F(qty,qty,nbins,h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

  for(int i=0; i<n_events_rdm;++i) {
    double rdm_number = rdm->Uniform(1);
    for(int j=1; j<=nbins;++j) {
      if(h_cdf->GetBinContent(j)>rdm_number) {
        h_pseudo->Fill(h_cdf->GetXaxis()->GetBinCenter(j));
        break;
      }
    }
  }

return;
}
