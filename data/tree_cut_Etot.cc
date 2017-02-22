#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void tree_cut_Etot() {

  TString input_file =  "./source/CD/0nu_rhc.root";
  TString output_file = "./source/CD/0nu_rhc_Ecut.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");

  Long64_t nentries = tree->GetEntries();
  double Esum = -1;

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
    if(Esum>2)
      newtree->Fill();
  }

  newtree->Print();
  newfile->Write();

  delete f;
  delete newfile;

  return;
}
