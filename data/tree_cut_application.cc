#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void tree_cut_application() {

  // int n_training =  17933000;
  // int n_application = 17933000;

  int n_training =  37500;
  int n_application = 37500;

  // TString input_file = "./raw_trees/0nu/merge.root";
  // TString output_file = "./trees_source_2e_application/0nu_1M_bis.root";
  // TString output_file = "./trees_topo_cuts_2e_application/0nu.root";

  // TString input_file = "./raw_trees/2nu_full/merge.root";
  // // TString output_file = "./trees_source_2e_application/2nu_full.root";
  // TString output_file = "./trees_source_2e_application/2nu_full_2MeV.root";
  // // TString output_file = "./trees_topo_cuts_2e_application/2nu.root";

  // TString input_file = "./raw_trees/tl208/merge.root";
  // TString output_file = "./trees_source_2e_application/tl208.root";
  // // TString output_file = "./trees_topo_cuts_2e_application/tl208.root";

  // TString input_file = "./raw_trees/bi214/merge.root";
  // TString output_file = "./trees_source_2e_application/bi214.root";
  // // TString output_file = "./trees_topo_cuts_2e_application/bi214.root";

  TString input_file = "./raw_trees/radon/merge.root";
  TString output_file = "./trees_source_2e_application/radon.root";
  // TString output_file = "./trees_topo_cuts_2e_application/radon.root";

  // TString input_file = "./raw_trees/2nu_full/merge.root";
  // TString output_file = "./trees_source_2e_application/2nu_full_1M.root";
  // // TString output_file = "./trees_topo_cuts_2e_application/2nu.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");

  Long64_t nentries = tree->GetEntries();
  double vertex_location = -1;
  double pint = -1;
  double delta_y = -1;
  double delta_z = -1;

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
  int count = 0;
  int count_application = 0;
  for (Long64_t i=0;i<nentries; i++) {
    if (i%100000==0)
      std::cout << i << std::endl;
    tree->GetEntry(i);

    //source
    if (vertex_location==0) {
      ++count;
      if (count >= n_training) {
        ++count_application;
        // if(Esum>2)
          newtree->Fill();
      }
    }

    // //topo
    // if (pint>cut_pint &&
    //     delta_y<cut_delta_y &&
    //     delta_z<cut_delta_z &&
    //     vertex_location==cut_vertex_location) {
    //   ++count;
    //   if (count >= n_training) {
    //     ++count_application;
    //    wtree->Fill();
    //   }
    // }

    if(count_application == n_application)
      break;
  }

  newtree->Print();
  newfile->Write();

  delete f;
  delete newfile;

  return;
}
