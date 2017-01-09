#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

// void tree_cut(TString input_file ="", TString output_file = "",
//               double cut_vertex_location = 0, double cut_pint = 0.04, double cut_delta_y = 60, double cut_delta_z = 70) {
void tree_cut() {

  // if (input_file == "" || output_file == "") {
  //   std::cout << "ERROR : No input/output file(s) specified !"  << std::endl;
  //   return;
  // }

  // TString input_file = "./raw_trees/0nu/merge.root";
  // // TString output_file = "./trees_source_2e/0nu.root";
  // TString output_file = "./trees_topo_cuts_2e/0nu.root";

  // TString input_file = "./raw_trees/2nu_full/merge.root";
  // // TString output_file = "./trees_source_2e/2nu_full.root";
  // TString output_file = "./trees_topo_cuts_2e/2nu_full.root";

  // TString input_file = "./raw_trees/2nu/merge.root";
  // // TString output_file = "./trees_source_2e/2nu_2MeV.root";
  // TString output_file = "./trees_topo_cuts_2e/2nu_2MeV.root";

  // TString input_file = "./raw_trees/tl208/merge.root";
  // // TString output_file = "./trees_source_2e/tl208.root";
  // TString output_file = "./trees_topo_cuts_2e/tl208.root";

  // TString input_file = "./raw_trees/bi214/merge.root";
  // // TString output_file = "./trees_source_2e/bi214.root";
  // TString output_file = "./trees_topo_cuts_2e/bi214.root";

  TString input_file = "./raw_trees/radon/merge.root";
  TString output_file = "./trees_source_2e/radon.root";
  // TString output_file = "./trees_topo_cuts_2e/radon.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");

  Long64_t nentries = tree->GetEntries();
  double vertex_location = -1;
  double pint = -1;
  double delta_y = 999;
  double delta_z = 999;

  double Esum = -1;

  double cut_vertex_location = 0;
  double cut_pint = 0.04;
  double cut_delta_y = 60;
  double cut_delta_z = 70;

  tree->SetBranchAddress("2e_electrons_vertex_location",&vertex_location);
  tree->SetBranchAddress("2e_electrons_internal_probability",&pint);
  tree->SetBranchAddress("2e_electrons_vertices_distance_y",&delta_y);
  tree->SetBranchAddress("2e_electrons_vertices_distance_z",&delta_z);
  tree->SetBranchAddress("2e_electrons_energy_sum",&Esum);

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("2e_*",1);

  TFile *newfile = new TFile(output_file,"recreate");
  TTree *newtree = tree->CloneTree(0);
  for (Long64_t i=0;i<nentries; i++) {
    if (i%100000==0)
      std::cout << i << std::endl;
    tree->GetEntry(i);

    // if (pint>cut_pint &&
    //     delta_y<cut_delta_y &&
    //     delta_z<cut_delta_z &&
    //     vertex_location==cut_vertex_location)
    //   newtree->Fill();

    if (vertex_location==0) {
      newtree->Fill();
    }
  }

  newtree->Print();
  newfile->Write();

  delete f;
  delete newfile;

  return;
}
