#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void tree_cut_subset() {

  // Topo cuts
  // // 0nu : 1/4
  // int n_start =    3*5947500;
  // int n_end = 4*5947500;

  // // 0nu : 1/4
  // int n_start = 3*1000000;
  // int n_end =   4*1000000;

  // // 2nu
  // int n_start = 3*7880000;
  // int n_end =   4*7880000;

  // // tl208
  // int n_start = 3*87750;
  // int n_end =   4*87750;

  // //bi214
  // int n_start = 3*138000;
  // int n_end =   4*138000;

  // // radon
  // int n_start = 3*5500;
  // int n_end =   4*5500;


  // TString input_file = "./trees_topo_cuts_2e/0nu.root";
  // TString output_file = "./topo_cuts/D/0nu.root";

  // TString input_file = "./trees_topo_cuts_2e/0nu.root";
  // TString output_file = "./topo_cuts/D/0nu.root";

  // TString input_file = "./trees_topo_cuts_2e/2nu_full.root";
  // TString output_file = "./topo_cuts/D/2nu_full.root";

  // TString input_file = "./trees_topo_cuts_2e/tl208.root";
  // TString output_file = "./topo_cuts/D/tl208.root";

  // TString input_file = "./trees_topo_cuts_2e/bi214.root";
  // TString output_file = "./topo_cuts/D/bi214.root";

  // TString input_file = "./trees_topo_cuts_2e/radon.root";
  // TString output_file = "./topo_cuts/A/radon.root";


  // Source : divide in 4 samples A,B,C,D

  double i_low =  3;
  double i_high = 4;
  TString sample = "D";

  // 0nu
  int n_start = i_low*1000000;
  int n_end =   i_high*1000000;

  // // 2nu
  // int n_start = 3*8750000;
  // int n_end =   4*8750000;

  // // tl208
  // int n_start = i_low*105750;
  // int n_end =   i_high*105750;

  // //bi214
  // int n_start = i_low*162500;
  // int n_end =   i_high*162500;

  // // radon
  // int n_start = i_low*34000;
  // int n_end =   i_high*34000;

  // TString input_file = "./trees_source_2e/0nu.root";
  // TString output_file = "./source/" + sample + "/0nu.root";

  TString input_file = "./trees_source_2e/0nu_rhc.root";
  TString output_file = "./source/" + sample + "/0nu_rhc.root";

  // TString input_file = "./trees_source_2e/2nu_full.root";
  // TString output_file = "./source/D/2nu_full.root";

  // TString input_file = "./trees_source_2e/tl208.root";
  // TString output_file = "./source/"+ sample + "/tl208.root";

  // TString input_file = "./trees_source_2e/bi214.root";
  // TString output_file = "./source/" + sample + "/bi214.root";

  // TString input_file = "./trees_source_2e/radon.root";
  // TString output_file = "./source/" + sample + "/radon.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");

  Long64_t nentries = tree->GetEntries();

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

    ++count;
    if (count >= n_start) {
      ++count_application;
      newtree->Fill();
    }

    if(count_application == n_end-n_start)
      break;
  }

  newtree->Print();
  newfile->Write();

  delete f;
  delete newfile;

  return;
}
