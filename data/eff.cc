#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include "histo_gen.h"

void eff() {

  double norm = 1e8;

  // TString input_file = "./trees_source_2e/0nu.root";
  // TString input_file = "./trees_source_2e/2nu.root";
  // TString input_file = "./trees_source_2e/tl208.root";
  // TString input_file = "./trees_source_2e/bi214.root";
  // TString input_file = "./trees_source_2e/radon.root";
  // TString input_file = "./trees_source_2e/bi214_calo_wrapper_surface.root";
  // TString input_file = "./trees_source_2e/tl208_pmt_glass_bulk.root";
  // TString input_file = "./trees_source_2e/bi214_pmt_glass_bulk.root";
  // TString input_file = "./trees_source_2e/k40_pmt_glass_bulk.root";

  // TString input_file = "./trees_topo_cuts_2e/0nu.root";
  // TString input_file = "./trees_topo_cuts_2e/2nu.root";
  // TString input_file = "./trees_topo_cuts_2e/tl208.root";
  // TString input_file = "./trees_topo_cuts_2e/bi214.root";
  // TString input_file = "./trees_topo_cuts_2e/radon.root";
  // TString input_file = "./trees_topo_cuts_2e/bi214_calo_wrapper_surface.root";
  // TString input_file = "./trees_topo_cuts_2e/tl208_pmt_glass_bulk.root";
  // TString input_file = "./trees_topo_cuts_2e/bi214_pmt_glass_bulk.root";
  TString input_file = "./trees_topo_cuts_2e/k40_pmt_glass_bulk.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");

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

  tree->Project(qty,qty,"2e_electrons_energy_sum >= 2.7 && 2e_electrons_energy_sum <= 3.2");

  h->SetBinContent(0,0);
  h->SetBinContent(nbins+1,0);

  std::cout << " N_events " << h->Integral(1,100) << std::endl;
  std::cout << " N_events " << h->Integral(1,100) /norm << std::endl;
  // h->Sumw2();
  // h->SetDrawOption("APL");

  // f->Close();

  return;
}
