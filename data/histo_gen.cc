#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include "histo_gen.h"

void histo_gen(TString input_file = "", TString output_file = "") {

  // TString input_file = "./test_tree.root";
  // TString output_file = "./test_histo.root";

  if (input_file == "" || output_file == "") {
    std::cout << "ERROR : No input/output file(s) specified !"  << std::endl;
    return;
  }

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");
  TFile *f_output= new TFile(output_file, "RECREATE");

  int nbins;
  double xmin, xmax;
  //default values
  nbins = 100;
  xmin = 0;
  xmax = 5;

  TString qty = "2e_electrons_energy_sum";
  get_histogram_options(qty,nbins,xmin,xmax);

  TH1F* h = new TH1F(qty,qty,nbins,xmin,xmax);

  // tree->Project(qty,qty,"2e_electrons_internal_probability >= 0.04 && 2e_electrons_vertices_distance_y <= 60 && 2e_electrons_vertices_distance_z <= 70 && 2e_electrons_vertex_location == 0");

  tree->Project(qty,qty);

  h->SetBinContent(0,0);
  h->SetBinContent(nbins+1,0);

  h->Sumw2();
  h->SetDrawOption("APL");
  h->Write();

  f->Close();
  f_output->Close();

  return;
}
